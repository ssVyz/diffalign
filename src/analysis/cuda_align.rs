//! CUDA-accelerated bitap aligner backend.
//!
//! Mirrors the API surface of `simple_align.rs` and `simd_align.rs` so it
//! plugs into the `WorkerAligner` dispatch in `screener.rs` without further
//! changes:
//!
//! * `CudaAligner` — per-worker handle. The actual GPU state is global; this
//!   struct is essentially a marker carried by Rayon's `map_init`.
//! * `create_cuda_aligner` — constructed once per worker (cheap, no GPU work).
//! * `collect_matches_with_cuda_aligner` / `collect_mismatch_counts_with_cuda_aligner`
//!
//! ## Architecture
//!
//! The GPU is a single shared resource. We initialize one CUDA context per
//! process and serialize all kernel calls through a global `Mutex`. Per-window
//! work scales with the GPU's internal parallelism (one thread per reference),
//! not with the number of CPU workers; the workers contend on the mutex but
//! their CPU-side work (variant grouping, fragment extraction) still overlaps.
//!
//! ## Reference upload
//!
//! References are uploaded to the GPU exactly once per run. The screener
//! (`run_screening` in screener.rs) calls [`register_slot`] up front for
//! the main reference set (slot 0) and, if differential mode is on, the
//! exclusivity set (slot 1). Per-window dispatches identify which slot to
//! use by the base pointer + length of the `&[Vec<u8>]` slice they receive —
//! the same slice the screener passed at registration time. If a slice
//! reaches the dispatch that wasn't pre-registered, we return an error
//! result (treated as "all no-match" so the run doesn't abort silently).

use std::os::raw::c_char;
use std::sync::Mutex;

use once_cell::sync::Lazy;

use super::simplescreen::pattern::{MAX_PATTERN_LEN, PreparedPattern};
use super::simplescreen::screener::Orientation;
use super::types::{AnchorHit, SimpleParams};

// ───── FFI: declarations match cuda/diffalign_cuda.h ──────────────────────

const MAX_K: u32 = 16;

#[repr(C)]
#[derive(Copy, Clone, Default)]
struct DiffalignCudaHit {
    end_pos: u32,
    mismatches: u32,
    found: u8,
    _pad: [u8; 3],
}

#[link(name = "diffalign_cuda", kind = "static")]
unsafe extern "C" {
    fn diffalign_cuda_init(err_buf: *mut c_char, err_buf_len: i32) -> i32;
    fn diffalign_cuda_shutdown();
    fn diffalign_cuda_upload(
        slot: i32,
        concat: *const u8,
        concat_len: u64,
        offsets: *const u32,
        num_refs: u32,
        err_buf: *mut c_char,
        err_buf_len: i32,
    ) -> i32;
    fn diffalign_cuda_slot_count(slot: i32) -> u32;
    fn diffalign_cuda_scan(
        slot: i32,
        masks: *const u64,
        pattern_len: u32,
        max_mismatches: u32,
        out: *mut DiffalignCudaHit,
        err_buf: *mut c_char,
        err_buf_len: i32,
    ) -> i32;
}

// ───── global state ───────────────────────────────────────────────────────

/// Identity of a registered reference set, captured at upload time. Pointer
/// is stored as `usize` so the struct stays `Send` (raw pointers aren't).
#[derive(Copy, Clone, Default)]
struct SlotIdent {
    base_ptr: usize,
    len: usize,
    /// Cached copy of `num_refs` for sanity checks; equals `len`.
    num_refs: u32,
}

struct CudaState {
    initialized: bool,
    slot0: Option<SlotIdent>,
    slot1: Option<SlotIdent>,
}

static STATE: Lazy<Mutex<CudaState>> = Lazy::new(|| {
    Mutex::new(CudaState {
        initialized: false,
        slot0: None,
        slot1: None,
    })
});

fn capture_err(buf: &[c_char]) -> String {
    // Find the null terminator (or take the whole buffer if none).
    let mut end = buf.len();
    for (i, &b) in buf.iter().enumerate() {
        if b == 0 {
            end = i;
            break;
        }
    }
    let bytes: Vec<u8> = buf[..end].iter().map(|&c| c as u8).collect();
    String::from_utf8_lossy(&bytes).into_owned()
}

/// Initialize the CUDA backend if it hasn't been. Idempotent; safe to call
/// many times.
pub fn ensure_initialized() -> Result<(), String> {
    let mut s = STATE.lock().map_err(|e| e.to_string())?;
    if s.initialized {
        return Ok(());
    }
    let mut err_buf = [0 as c_char; 512];
    let rc = unsafe { diffalign_cuda_init(err_buf.as_mut_ptr(), err_buf.len() as i32) };
    if rc != 0 {
        return Err(capture_err(&err_buf));
    }
    s.initialized = true;
    Ok(())
}

/// Upload a reference set into `slot` (0 = main, 1 = exclusivity). The slot's
/// identity is recorded so per-window dispatches can match the slice they
/// receive against the registered slot.
pub fn register_slot(slot: i32, references: &[Vec<u8>]) -> Result<(), String> {
    ensure_initialized()?;

    // Flatten host-side: concat bytes + offsets table.
    let mut offsets: Vec<u32> = Vec::with_capacity(references.len() + 1);
    let mut total: u64 = 0;
    offsets.push(0);
    for r in references {
        total += r.len() as u64;
        if total > u32::MAX as u64 {
            return Err(format!(
                "concatenated reference size exceeds 4 GiB (limit of the u32 offset table); \
                 got {} bytes after {} references",
                total,
                offsets.len()
            ));
        }
        offsets.push(total as u32);
    }
    let mut concat: Vec<u8> = Vec::with_capacity(total as usize);
    for r in references {
        concat.extend_from_slice(r);
    }

    let mut err_buf = [0 as c_char; 512];
    let rc = unsafe {
        diffalign_cuda_upload(
            slot,
            concat.as_ptr(),
            total,
            offsets.as_ptr(),
            references.len() as u32,
            err_buf.as_mut_ptr(),
            err_buf.len() as i32,
        )
    };
    if rc != 0 {
        return Err(capture_err(&err_buf));
    }

    let ident = SlotIdent {
        base_ptr: references.as_ptr() as usize,
        len: references.len(),
        num_refs: references.len() as u32,
    };
    let mut s = STATE.lock().map_err(|e| e.to_string())?;
    match slot {
        0 => s.slot0 = Some(ident),
        1 => s.slot1 = Some(ident),
        other => return Err(format!("invalid slot {}", other)),
    }
    Ok(())
}

/// Identify which (pre-registered) slot a `&[Vec<u8>]` slice corresponds to.
fn slot_for(references: &[Vec<u8>]) -> Option<i32> {
    let p = references.as_ptr() as usize;
    let l = references.len();
    let s = STATE.lock().ok()?;
    if let Some(ident) = s.slot0 {
        if ident.base_ptr == p && ident.len == l {
            return Some(0);
        }
    }
    if let Some(ident) = s.slot1 {
        if ident.base_ptr == p && ident.len == l {
            return Some(1);
        }
    }
    None
}

/// Issue a single-orientation scan kernel call. Returns one hit per ref.
fn scan_one(slot: i32, masks: &[u64; 256], pattern_len: u32, k: u32) -> Result<Vec<DiffalignCudaHit>, String> {
    let num_refs = unsafe { diffalign_cuda_slot_count(slot) } as usize;
    let mut out: Vec<DiffalignCudaHit> = vec![DiffalignCudaHit::default(); num_refs];

    // Serialize GPU access through the global state mutex; held only for the
    // launch + sync window.
    let _guard = STATE.lock().map_err(|e| e.to_string())?;
    let mut err_buf = [0 as c_char; 512];
    let rc = unsafe {
        diffalign_cuda_scan(
            slot,
            masks.as_ptr(),
            pattern_len,
            k,
            out.as_mut_ptr(),
            err_buf.as_mut_ptr(),
            err_buf.len() as i32,
        )
    };
    if rc != 0 {
        return Err(capture_err(&err_buf));
    }
    Ok(out)
}

// ───── public aligner API ─────────────────────────────────────────────────

/// Per-worker scratch state for the CUDA aligner. Carries no GPU resources;
/// the GPU state is global. `Send`, not `Sync` (we let Rayon hand one per
/// worker for API symmetry with the other backends).
pub struct CudaAligner {
    /// Reusable buffer of merged per-orientation results.
    _scratch: (),
}

pub fn create_cuda_aligner(_params: &SimpleParams) -> CudaAligner {
    CudaAligner { _scratch: () }
}

/// Largest oligo length the CUDA backend can handle. Same single-u64 cap as
/// the scalar/SIMD bitap.
pub const CUDA_MAX_OLIGO_LEN: usize = MAX_PATTERN_LEN;

pub const CUDA_MAX_MISMATCHES: u32 = MAX_K;

#[inline]
fn complement_byte(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' | b'u' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

fn reverse_complement(bytes: &[u8]) -> String {
    let mut out = Vec::with_capacity(bytes.len());
    for &b in bytes.iter().rev() {
        out.push(complement_byte(b));
    }
    String::from_utf8(out).unwrap_or_default()
}

/// Per-ref best hit after merging fwd + rc orientations.
struct MergedHit {
    mismatches: u32,
    end_pos: u32,
    orientation: Orientation,
}

fn merge_orientation(
    fwd: DiffalignCudaHit,
    rc: DiffalignCudaHit,
) -> Option<MergedHit> {
    match (fwd.found, rc.found) {
        (0, 0) => None,
        (1, 0) => Some(MergedHit {
            mismatches: fwd.mismatches,
            end_pos: fwd.end_pos,
            orientation: Orientation::Forward,
        }),
        (0, 1) => Some(MergedHit {
            mismatches: rc.mismatches,
            end_pos: rc.end_pos,
            orientation: Orientation::Reverse,
        }),
        _ => {
            // Both fired: same tiebreak as simple_align::best_hit —
            // (min mismatches, earliest start [== earliest end_pos for fixed
            // length], forward preferred over reverse).
            let pick_fwd = if fwd.mismatches != rc.mismatches {
                fwd.mismatches < rc.mismatches
            } else if fwd.end_pos != rc.end_pos {
                fwd.end_pos < rc.end_pos
            } else {
                true
            };
            if pick_fwd {
                Some(MergedHit {
                    mismatches: fwd.mismatches,
                    end_pos: fwd.end_pos,
                    orientation: Orientation::Forward,
                })
            } else {
                Some(MergedHit {
                    mismatches: rc.mismatches,
                    end_pos: rc.end_pos,
                    orientation: Orientation::Reverse,
                })
            }
        }
    }
}

/// Run the GPU scan for one window, returning per-reference merged best hits.
/// Errors out (rather than silently degrading) if the references slice wasn't
/// pre-registered or the kernel call fails.
fn scan_window(
    pattern: &PreparedPattern,
    references: &[Vec<u8>],
) -> Result<Vec<Option<MergedHit>>, String> {
    let slot = slot_for(references)
        .ok_or_else(|| "references slice not registered with CUDA backend".to_string())?;

    let fwd = scan_one(slot, &pattern.fwd_masks, pattern.len, pattern.max_mismatches)?;
    let rc  = scan_one(slot, &pattern.rc_masks,  pattern.len, pattern.max_mismatches)?;

    if fwd.len() != references.len() || rc.len() != references.len() {
        return Err(format!(
            "CUDA scan returned {}/{} hits, expected {} (slot mismatch?)",
            fwd.len(), rc.len(), references.len()
        ));
    }

    let merged: Vec<Option<MergedHit>> = (0..references.len())
        .map(|i| merge_orientation(fwd[i], rc[i]))
        .collect();
    Ok(merged)
}

pub fn collect_matches_with_cuda_aligner(
    _aligner: &mut CudaAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &SimpleParams,
) -> (Vec<String>, usize) {
    let pattern = match PreparedPattern::build(oligo, params.max_mismatches) {
        Ok(p) => p,
        Err(_) => return (Vec::new(), references.len()),
    };
    if !pattern.valid {
        return (Vec::new(), references.len());
    }

    let merged = match scan_window(&pattern, references) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("CUDA scan failed: {}", e);
            return (Vec::new(), references.len());
        }
    };

    let len = pattern.len as usize;
    let mut matched = Vec::new();
    let mut no_match_count = 0usize;
    for (i, hit) in merged.into_iter().enumerate() {
        match hit {
            Some(h) => {
                let end_0 = (h.end_pos as usize) + 1;
                let start_0 = end_0 - len;
                let frag = &references[i][start_0..end_0];
                let s = match h.orientation {
                    Orientation::Forward => String::from_utf8_lossy(frag).into_owned(),
                    Orientation::Reverse => reverse_complement(frag),
                };
                matched.push(s);
            }
            None => no_match_count += 1,
        }
    }
    (matched, no_match_count)
}

/// Per-reference anchor positions from a single GPU scan window. Same accept
/// rules as `collect_matches_with_cuda_aligner`; the returned start is
/// 0-based on the reference's forward strand.
pub fn collect_anchors_with_cuda_aligner(
    _aligner: &mut CudaAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &SimpleParams,
) -> Vec<Option<AnchorHit>> {
    let pattern = match PreparedPattern::build(oligo, params.max_mismatches) {
        Ok(p) => p,
        Err(_) => return references.iter().map(|_| None).collect(),
    };
    if !pattern.valid {
        return references.iter().map(|_| None).collect();
    }

    let merged = match scan_window(&pattern, references) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("CUDA scan failed: {}", e);
            return references.iter().map(|_| None).collect();
        }
    };

    let len = pattern.len as usize;
    merged
        .into_iter()
        .map(|h| {
            h.map(|m| {
                let end_0 = (m.end_pos as usize) + 1;
                let start_0 = end_0 - len;
                AnchorHit {
                    start: start_0,
                    orientation: m.orientation,
                    mismatches: m.mismatches,
                }
            })
        })
        .collect()
}

pub fn collect_mismatch_counts_with_cuda_aligner(
    _aligner: &mut CudaAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &SimpleParams,
) -> Vec<Option<u32>> {
    let pattern = match PreparedPattern::build(oligo, params.max_mismatches) {
        Ok(p) => p,
        Err(_) => return references.iter().map(|_| None).collect(),
    };
    if !pattern.valid {
        return references.iter().map(|_| None).collect();
    }

    let merged = match scan_window(&pattern, references) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("CUDA scan failed: {}", e);
            return references.iter().map(|_| None).collect();
        }
    };

    merged.into_iter().map(|h| h.map(|m| m.mismatches)).collect()
}

/// Tear down the GPU context. Optional but tidy at program exit.
#[allow(dead_code)]
pub fn shutdown() {
    let _ = STATE.lock().map(|mut s| {
        if s.initialized {
            unsafe { diffalign_cuda_shutdown() };
            s.initialized = false;
            s.slot0 = None;
            s.slot1 = None;
        }
    });
}

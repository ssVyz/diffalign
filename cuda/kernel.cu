// Bitap (Wu-Manber substitutions-only) kernel for the diffalign CUDA backend.
//
// One thread = one reference. Each thread runs the same recurrence as the
// scalar bitap in `src/analysis/simplescreen/bitap.rs`, producing exactly the
// same per-position registers `R[d]` and reporting the smallest `d` with
// `R[d] & end_bit != 0` at each text position. The per-reference best hit is
// kept by (min mismatches, earliest end position) — same tiebreak as the
// scalar/SIMD `LaneBest`.
//
// The mask table is staged into shared memory at block start so the inner
// loop hits ~zero global-memory latency for mask lookups. The reference data
// itself is read from global memory once per byte per thread; no caching is
// done since the next byte is consumed immediately.
//
// Ragged-length handling: each thread iterates its own reference length only.
// Within a warp this means short-reference threads finish early and idle
// while longer-reference threads continue. The host can mitigate this by
// sorting references by length before upload; this kernel makes no
// assumption about input order.

#include <cuda_runtime.h>
#include <stdint.h>

#include "diffalign_cuda.h"

#ifndef DIFFALIGN_CUDA_THREADS_PER_BLOCK
#define DIFFALIGN_CUDA_THREADS_PER_BLOCK 128
#endif

// Fixed register-array size for the per-thread bitap state. Compile-time so
// the array can live in registers (or local memory on overflow). Capped at
// DIFFALIGN_CUDA_MAX_K + 1.
#define DIFFALIGN_CUDA_REGS (DIFFALIGN_CUDA_MAX_K + 1)

extern "C" __global__ void diffalign_bitap_kernel(
    const uint8_t* __restrict__  concat,
    const uint32_t* __restrict__ offsets,
    uint32_t                      num_refs,
    const uint64_t* __restrict__ masks_global,
    uint32_t                      pattern_len,
    uint32_t                      k,           // max_mismatches
    uint64_t                      end_bit,
    DiffalignCudaHit*             out)
{
    // Stage the 256-entry mask table into shared memory.
    __shared__ uint64_t masks_shared[256];

    // 128 threads/block, 256 mask slots → each thread copies 2 entries.
    for (uint32_t i = threadIdx.x; i < 256; i += blockDim.x) {
        masks_shared[i] = masks_global[i];
    }
    __syncthreads();

    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_refs) return;

    uint32_t start = offsets[tid];
    uint32_t end   = offsets[tid + 1];
    uint32_t len   = end - start;

    // Per-thread bitap state. We allocate the maximum and only iterate up to
    // `k`; the extra slots cost a few registers but keep the loop bounds
    // static, which lets the compiler unroll the d-loop on small k.
    uint64_t regs[DIFFALIGN_CUDA_REGS];
    #pragma unroll
    for (uint32_t d = 0; d < DIFFALIGN_CUDA_REGS; d++) regs[d] = 0ULL;

    DiffalignCudaHit best;
    best.end_pos    = 0u;
    best.mismatches = 0u;
    best.found      = 0u;
    best._pad[0] = best._pad[1] = best._pad[2] = 0u;

    const uint8_t* ref = concat + start;

    for (uint32_t j = 0; j < len; j++) {
        uint8_t c = ref[j];
        uint64_t cm = masks_shared[c];

        // Update R[d] from d = k down to d = 1, then R[0] last. Same order
        // as the scalar bitap; ensures each R[d] still sees the previous
        // step's R[d-1] when computing its new value.
        // (Signed loop variable so d=1, d--, d=0, condition fails — no
        // underflow concerns.)
        for (int32_t d = (int32_t)k; d >= 1; d--) {
            uint64_t shifted_d   = regs[d]     << 1;
            uint64_t shifted_dm1 = regs[d - 1] << 1;
            regs[d] = ((shifted_d & cm) | shifted_dm1) | 1ULL;
        }
        regs[0] = ((regs[0] << 1) | 1ULL) & cm;

        // Match check: regs[k] is the union of all R[0..=k] (since
        // R[d] ⊆ R[d+1] in this recurrence).
        if ((regs[k] & end_bit) != 0ULL) {
            // Find the smallest d with the end bit set — that's the actual
            // mismatch count of this hit.
            uint32_t mm = 0;
            for (uint32_t d = 0; d <= k; d++) {
                if ((regs[d] & end_bit) != 0ULL) {
                    mm = d;
                    break;
                }
            }

            // Update best by (min mismatches, earliest end position).
            bool better = false;
            if (!best.found) {
                better = true;
            } else if (mm < best.mismatches) {
                better = true;
            } else if (mm == best.mismatches && j < best.end_pos) {
                better = true;
            }
            if (better) {
                best.found      = 1u;
                best.mismatches = mm;
                best.end_pos    = j;
            }
        }
    }

    out[tid] = best;
}

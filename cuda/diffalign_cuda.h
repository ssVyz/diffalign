/* C ABI for the diffalign CUDA bitap backend.
 *
 * The whole API is process-global: a single CUDA context manages two
 * reference-set slots (0 = main references, 1 = exclusivity references).
 * The host (Rust) side serializes calls through a mutex; this header makes
 * no internal threading guarantees beyond that.
 *
 * Error reporting: every fallible function returns 0 on success, non-zero on
 * error. On non-zero return, a human-readable message is written to `err_buf`
 * (null-terminated, truncated to `err_buf_len - 1` bytes). Pass NULL/0 to
 * skip the message.
 */
#ifndef DIFFALIGN_CUDA_H
#define DIFFALIGN_CUDA_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Hard upper bound on the supported max_mismatches value for the CUDA kernel.
 * Matches the kernel's fixed-size register array. Anything beyond this is
 * rejected up front. Realistic DNA bitap uses values well below this. */
#define DIFFALIGN_CUDA_MAX_K 16

/* Per-reference best hit returned by `diffalign_cuda_scan`.
 *
 * `end_pos` is the 0-based index of the last reference byte covered by the
 * best hit. The start (0-based) is `end_pos + 1 - pattern_len`. If `found`
 * is 0, the other fields are undefined.
 */
typedef struct DiffalignCudaHit {
    uint32_t end_pos;
    uint32_t mismatches;
    uint8_t  found;       /* 0 or 1 */
    uint8_t  _pad[3];
} DiffalignCudaHit;

/* Initialize CUDA on device 0.
 *
 * Idempotent: subsequent calls return 0 without re-initializing. Returns
 * non-zero if no CUDA device is available or the context cannot be created.
 */
int32_t diffalign_cuda_init(char* err_buf, int32_t err_buf_len);

/* Tear down the global context, freeing all uploaded buffers. Safe to call
 * even if init was never called. */
void diffalign_cuda_shutdown(void);

/* Upload one reference set into slot `slot` (0 or 1).
 *
 * `concat` is the concatenation of all `num_refs` reference byte strings;
 * its length is `concat_len`. `offsets` has `num_refs + 1` entries —
 * reference `i` occupies `concat[offsets[i]..offsets[i+1]]`. Replaces
 * whatever was previously in `slot`. Empty references are allowed.
 */
int32_t diffalign_cuda_upload(
    int32_t        slot,
    const uint8_t* concat,
    uint64_t       concat_len,
    const uint32_t* offsets,
    uint32_t       num_refs,
    char*          err_buf,
    int32_t        err_buf_len);

/* Run the bitap kernel for one (pattern, orientation) against the
 * references in `slot`. `masks` is the 256-entry mask table built by the
 * caller (the same table the scalar/SIMD bitap uses). `out` must have space
 * for `num_refs` `DiffalignCudaHit` entries — caller knows num_refs from
 * upload time, so we don't repeat it. `pattern_len` must be in [1, 64];
 * `max_mismatches` in [0, DIFFALIGN_CUDA_MAX_K] and strictly less than
 * `pattern_len`.
 *
 * The kernel computes, per reference, the best hit by (min mismatches,
 * earliest end position) — the same per-orientation rule the scalar bitap
 * and SIMD bitap apply. Cross-orientation tie-breaking is done host-side by
 * the caller, who issues two scans (fwd, rc) and merges the per-ref results.
 */
int32_t diffalign_cuda_scan(
    int32_t            slot,
    const uint64_t*    masks,
    uint32_t           pattern_len,
    uint32_t           max_mismatches,
    DiffalignCudaHit*  out,
    char*              err_buf,
    int32_t            err_buf_len);

/* Number of references currently uploaded in `slot`. Returns 0 if the slot
 * is empty or `slot` is out of range. */
uint32_t diffalign_cuda_slot_count(int32_t slot);

#ifdef __cplusplus
}
#endif

#endif /* DIFFALIGN_CUDA_H */

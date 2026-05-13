// Host-side glue for the diffalign CUDA backend.
//
// Manages a process-global CUDA context with two reference-set slots and
// dispatches the bitap kernel. Caller (Rust) is responsible for serializing
// calls — no internal mutex.

#include <cuda_runtime.h>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <new>

#include "diffalign_cuda.h"

// Declared in kernel.cu.
extern "C" __global__ void diffalign_bitap_kernel(
    const uint8_t*   concat,
    const uint32_t*  offsets,
    uint32_t         num_refs,
    const uint64_t*  masks_global,
    uint32_t         pattern_len,
    uint32_t         k,
    uint64_t         end_bit,
    DiffalignCudaHit* out);

namespace {

constexpr int kSlotCount = 2;
constexpr int kThreadsPerBlock = 128;

struct Slot {
    uint8_t*  d_concat   = nullptr;
    uint64_t  concat_len = 0;
    uint32_t* d_offsets  = nullptr;  // (num_refs + 1) entries
    uint32_t  num_refs   = 0;
    // Result buffer sized to `num_refs`. Reused across calls; reallocated if
    // num_refs grows.
    DiffalignCudaHit* d_out      = nullptr;
    uint32_t          d_out_cap  = 0;
};

struct Context {
    bool initialized = false;
    int  device_id   = 0;
    Slot slots[kSlotCount];
    // Per-call mask buffer. 256 u64. Reused.
    uint64_t* d_masks = nullptr;
};

Context g_ctx;

void write_err(char* buf, int32_t buf_len, const char* msg) {
    if (!buf || buf_len <= 0) return;
    std::snprintf(buf, (size_t)buf_len, "%s", msg);
}

void write_err_fmt(char* buf, int32_t buf_len, const char* fmt, ...) {
    if (!buf || buf_len <= 0) return;
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, (size_t)buf_len, fmt, args);
    va_end(args);
}

bool check_cuda(cudaError_t err, char* err_buf, int32_t err_buf_len, const char* what) {
    if (err == cudaSuccess) return true;
    write_err_fmt(err_buf, err_buf_len, "CUDA error during %s: %s", what, cudaGetErrorString(err));
    return false;
}

void free_slot(Slot& s) {
    if (s.d_concat)  { cudaFree(s.d_concat);  s.d_concat  = nullptr; }
    if (s.d_offsets) { cudaFree(s.d_offsets); s.d_offsets = nullptr; }
    if (s.d_out)     { cudaFree(s.d_out);     s.d_out     = nullptr; }
    s.concat_len = 0;
    s.num_refs   = 0;
    s.d_out_cap  = 0;
}

}  // namespace

extern "C" {

int32_t diffalign_cuda_init(char* err_buf, int32_t err_buf_len) {
    if (g_ctx.initialized) return 0;

    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess) {
        write_err_fmt(err_buf, err_buf_len,
            "no CUDA driver / runtime available: %s",
            cudaGetErrorString(err));
        return 1;
    }
    if (device_count <= 0) {
        write_err(err_buf, err_buf_len,
            "no CUDA-capable GPU detected on this system");
        return 1;
    }

    err = cudaSetDevice(0);
    if (!check_cuda(err, err_buf, err_buf_len, "cudaSetDevice(0)")) return 1;

    // Lazy-allocate the per-call mask buffer.
    err = cudaMalloc(&g_ctx.d_masks, 256 * sizeof(uint64_t));
    if (!check_cuda(err, err_buf, err_buf_len, "cudaMalloc(masks)")) return 1;

    g_ctx.device_id   = 0;
    g_ctx.initialized = true;
    return 0;
}

void diffalign_cuda_shutdown(void) {
    if (!g_ctx.initialized) return;
    for (int i = 0; i < kSlotCount; i++) free_slot(g_ctx.slots[i]);
    if (g_ctx.d_masks) { cudaFree(g_ctx.d_masks); g_ctx.d_masks = nullptr; }
    g_ctx.initialized = false;
}

int32_t diffalign_cuda_upload(
    int32_t        slot,
    const uint8_t* concat,
    uint64_t       concat_len,
    const uint32_t* offsets,
    uint32_t       num_refs,
    char*          err_buf,
    int32_t        err_buf_len)
{
    if (!g_ctx.initialized) {
        write_err(err_buf, err_buf_len, "CUDA backend not initialized");
        return 1;
    }
    if (slot < 0 || slot >= kSlotCount) {
        write_err_fmt(err_buf, err_buf_len, "invalid slot %d (expected 0 or 1)", slot);
        return 1;
    }

    Slot& s = g_ctx.slots[slot];
    free_slot(s);

    // Allocate ref bytes (handle empty slot gracefully).
    if (concat_len > 0) {
        cudaError_t err = cudaMalloc(&s.d_concat, (size_t)concat_len);
        if (!check_cuda(err, err_buf, err_buf_len, "cudaMalloc(refs)")) return 1;
        err = cudaMemcpy(s.d_concat, concat, (size_t)concat_len, cudaMemcpyHostToDevice);
        if (!check_cuda(err, err_buf, err_buf_len, "cudaMemcpy(refs)")) {
            free_slot(s);
            return 1;
        }
    }
    s.concat_len = concat_len;

    // Allocate offsets (always num_refs + 1 entries, even if num_refs = 0).
    {
        size_t bytes = (size_t)(num_refs + 1) * sizeof(uint32_t);
        cudaError_t err = cudaMalloc(&s.d_offsets, bytes);
        if (!check_cuda(err, err_buf, err_buf_len, "cudaMalloc(offsets)")) {
            free_slot(s);
            return 1;
        }
        err = cudaMemcpy(s.d_offsets, offsets, bytes, cudaMemcpyHostToDevice);
        if (!check_cuda(err, err_buf, err_buf_len, "cudaMemcpy(offsets)")) {
            free_slot(s);
            return 1;
        }
    }
    s.num_refs = num_refs;

    // Pre-allocate the output buffer to num_refs.
    if (num_refs > 0) {
        cudaError_t err = cudaMalloc(&s.d_out, (size_t)num_refs * sizeof(DiffalignCudaHit));
        if (!check_cuda(err, err_buf, err_buf_len, "cudaMalloc(out)")) {
            free_slot(s);
            return 1;
        }
        s.d_out_cap = num_refs;
    }

    return 0;
}

uint32_t diffalign_cuda_slot_count(int32_t slot) {
    if (!g_ctx.initialized) return 0;
    if (slot < 0 || slot >= kSlotCount) return 0;
    return g_ctx.slots[slot].num_refs;
}

int32_t diffalign_cuda_scan(
    int32_t            slot,
    const uint64_t*    masks,
    uint32_t           pattern_len,
    uint32_t           max_mismatches,
    DiffalignCudaHit*  out,
    char*              err_buf,
    int32_t            err_buf_len)
{
    if (!g_ctx.initialized) {
        write_err(err_buf, err_buf_len, "CUDA backend not initialized");
        return 1;
    }
    if (slot < 0 || slot >= kSlotCount) {
        write_err_fmt(err_buf, err_buf_len, "invalid slot %d", slot);
        return 1;
    }
    if (pattern_len == 0 || pattern_len > 64) {
        write_err_fmt(err_buf, err_buf_len,
            "pattern_len %u out of range [1, 64]", pattern_len);
        return 1;
    }
    if (max_mismatches > (uint32_t)DIFFALIGN_CUDA_MAX_K) {
        write_err_fmt(err_buf, err_buf_len,
            "max_mismatches %u exceeds CUDA build cap of %d",
            max_mismatches, DIFFALIGN_CUDA_MAX_K);
        return 1;
    }
    if (max_mismatches >= pattern_len) {
        write_err_fmt(err_buf, err_buf_len,
            "max_mismatches (%u) must be < pattern_len (%u)",
            max_mismatches, pattern_len);
        return 1;
    }

    Slot& s = g_ctx.slots[slot];
    if (s.num_refs == 0) return 0;  // nothing to do; out buffer untouched.

    // Upload mask table for this call.
    cudaError_t err = cudaMemcpy(g_ctx.d_masks, masks, 256 * sizeof(uint64_t),
                                  cudaMemcpyHostToDevice);
    if (!check_cuda(err, err_buf, err_buf_len, "cudaMemcpy(masks)")) return 1;

    // Launch the kernel.
    uint64_t end_bit = 1ULL << (pattern_len - 1);
    uint32_t blocks  = (s.num_refs + kThreadsPerBlock - 1) / kThreadsPerBlock;

    diffalign_bitap_kernel<<<blocks, kThreadsPerBlock>>>(
        s.d_concat,
        s.d_offsets,
        s.num_refs,
        g_ctx.d_masks,
        pattern_len,
        max_mismatches,
        end_bit,
        s.d_out);

    err = cudaGetLastError();
    if (!check_cuda(err, err_buf, err_buf_len, "kernel launch")) return 1;

    // Copy results back. (cudaMemcpy is synchronous; no extra sync needed.)
    err = cudaMemcpy(out, s.d_out,
                     (size_t)s.num_refs * sizeof(DiffalignCudaHit),
                     cudaMemcpyDeviceToHost);
    if (!check_cuda(err, err_buf, err_buf_len, "cudaMemcpy(hits)")) return 1;

    return 0;
}

}  // extern "C"

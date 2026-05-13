//! Build script for the diffalign crate.
//!
//! When the `cuda` feature is enabled, this compiles the C++/CUDA sources in
//! `cuda/` into a static library and tells Cargo to link it (plus the CUDA
//! runtime). When the feature is disabled, build.rs is a no-op.

use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=cuda/kernel.cu");
    println!("cargo:rerun-if-changed=cuda/host.cu");
    println!("cargo:rerun-if-changed=cuda/diffalign_cuda.h");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=CUDA_PATH");

    if env::var_os("CARGO_FEATURE_CUDA").is_none() {
        return;
    }

    build_cuda();
}

fn build_cuda() {
    let target = env::var("TARGET").expect("TARGET");
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR"));
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR"));

    let is_windows = target.contains("windows");
    let is_msvc = target.contains("msvc");

    let cuda_path = find_cuda_path().unwrap_or_else(|| {
        panic!(
            "could not locate CUDA toolkit. Set the CUDA_PATH environment \
             variable (e.g. C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v13.2 \
             on Windows) or install the CUDA toolkit."
        )
    });

    let nvcc_exe = cuda_path
        .join("bin")
        .join(if is_windows { "nvcc.exe" } else { "nvcc" });
    if !nvcc_exe.exists() {
        panic!(
            "nvcc not found at {} — check CUDA_PATH",
            nvcc_exe.display()
        );
    }

    let cuda_src = manifest_dir.join("cuda");
    let kernel = cuda_src.join("kernel.cu");
    let host = cuda_src.join("host.cu");

    let lib_name = "diffalign_cuda";
    let lib_filename = if is_msvc {
        format!("{}.lib", lib_name)
    } else {
        format!("lib{}.a", lib_name)
    };
    let lib_out = out_dir.join(&lib_filename);

    // Locate cl.exe via the `cc` crate on Windows MSVC targets and put its
    // directory on PATH for the nvcc invocation. nvcc requires a host compiler
    // (cl.exe on MSVC) and won't find it otherwise since VS doesn't ship cl
    // on the global PATH.
    let mut extra_path: Option<PathBuf> = None;
    if is_msvc {
        let compiler = cc::Build::new()
            .target(&target)
            .host(&env::var("HOST").unwrap_or_else(|_| target.clone()))
            .opt_level(3)
            .cargo_metadata(false)
            .cargo_warnings(false)
            .get_compiler();
        let cl_path = compiler.path();
        if let Some(parent) = cl_path.parent() {
            extra_path = Some(parent.to_path_buf());
        }
    }

    // CUDA 13.x dropped Pascal (sm_61) support; minimum is sm_75 (Turing).
    // If you need to run on Pascal hardware, install CUDA 12.x and add
    // `-gencode=arch=compute_61,code=sm_61` back to this list.
    let gencode_flags: &[&str] = &[
        // Turing: RTX 20-series.
        "-gencode=arch=compute_75,code=sm_75",
        // Ampere: RTX 30-series including 3050 Ti, A100.
        "-gencode=arch=compute_86,code=sm_86",
        // Ada Lovelace: RTX 40-series.
        "-gencode=arch=compute_89,code=sm_89",
        // Hopper: H100.
        "-gencode=arch=compute_90,code=sm_90",
        // Blackwell (consumer): RTX 50-series.
        "-gencode=arch=compute_120,code=sm_120",
        // Forward-compat: PTX for the latest virtual arch; JIT-compiled at
        // runtime for any newer real GPU.
        "-gencode=arch=compute_120,code=compute_120",
    ];

    let mut cmd = Command::new(&nvcc_exe);
    cmd.arg("-lib")
        .arg("-o")
        .arg(&lib_out)
        .arg("-std=c++17")
        .arg("-O3")
        .arg("-I")
        .arg(&cuda_src);

    for flag in gencode_flags {
        cmd.arg(flag);
    }

    if is_msvc {
        // Match MSVC's CRT used by the Rust msvc toolchain (dynamic /MD).
        cmd.arg("--compiler-options").arg("/MD /EHsc /O2 /nologo");
    } else {
        cmd.arg("--compiler-options").arg("-fPIC -O3");
    }

    cmd.arg(&kernel).arg(&host);

    if let Some(p) = extra_path.as_ref() {
        let cur = env::var_os("PATH").unwrap_or_default();
        let mut paths: Vec<PathBuf> = env::split_paths(&cur).collect();
        paths.insert(0, p.clone());
        let joined = env::join_paths(paths).expect("join PATH");
        cmd.env("PATH", joined);
    }

    eprintln!("running: {:?}", cmd);
    let status = cmd.status().expect("failed to invoke nvcc");
    if !status.success() {
        panic!("nvcc failed with status {}", status);
    }
    if !lib_out.exists() {
        panic!("expected nvcc output {} not found", lib_out.display());
    }

    // Tell Cargo to link the static lib and the CUDA runtime.
    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static={}", lib_name);

    let cuda_libdir = if is_windows {
        cuda_path.join("lib").join("x64")
    } else {
        cuda_path.join("lib64")
    };
    println!(
        "cargo:rustc-link-search=native={}",
        cuda_libdir.display()
    );
    println!("cargo:rustc-link-lib=dylib=cudart");

    // MSVC needs the C++ runtime when the static lib references it. The
    // compiler driver normally handles this, but with link.exe driven by
    // rustc, we hint at the standard libs nvcc-compiled MSVC objects need.
    if !is_msvc {
        println!("cargo:rustc-link-lib=dylib=stdc++");
    }
}

fn find_cuda_path() -> Option<PathBuf> {
    if let Ok(p) = env::var("CUDA_PATH") {
        let path = PathBuf::from(p);
        if path.exists() {
            return Some(path);
        }
    }
    if let Ok(p) = env::var("CUDA_HOME") {
        let path = PathBuf::from(p);
        if path.exists() {
            return Some(path);
        }
    }
    // Last resort: probe a couple of common locations.
    #[cfg(windows)]
    {
        let prog = env::var("ProgramFiles")
            .unwrap_or_else(|_| String::from("C:\\Program Files"));
        let root = PathBuf::from(prog).join("NVIDIA GPU Computing Toolkit").join("CUDA");
        if root.exists() {
            // Pick the highest-version subdir.
            if let Ok(entries) = std::fs::read_dir(&root) {
                let mut versions: Vec<PathBuf> =
                    entries.flatten().map(|e| e.path()).filter(|p| p.is_dir()).collect();
                versions.sort();
                if let Some(latest) = versions.into_iter().last() {
                    return Some(latest);
                }
            }
        }
    }
    None
}

use std::fs;
use std::io::{self, IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, anyhow, bail};
use clap::Parser;
use mimalloc::MiMalloc;

use diffalign::analysis::{
    ReferenceData, ScreeningResults, parse_reference_fasta, parse_template_fasta, run_screening,
};
use diffalign::cli::Cli;
use diffalign::config::{self, Config};
use diffalign::key_listener::KeyListener;
use diffalign::pause::PauseFlag;
use diffalign::progress::{Reporter, build_length_plan};

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() {
    if let Err(err) = run() {
        eprintln!("error: {:#}", err);
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let cli = Cli::parse();
    let started = Instant::now();
    let timer_enabled = cli.timer;

    let ini_path = match &cli.config {
        Some(path) => path.clone(),
        None => config::ini_path()?,
    };

    if cli.mkini {
        return cmd_mkini(&ini_path);
    }
    if cli.check_config {
        return cmd_check_config(&ini_path);
    }

    let cfg = match load_or_prompt(&ini_path)? {
        Some(cfg) => cfg,
        None => return Ok(()),
    };

    let template_path = cli
        .template
        .as_ref()
        .ok_or_else(|| anyhow!("template FASTA is required"))?;
    let reference_path = cli
        .reference
        .as_ref()
        .ok_or_else(|| anyhow!("reference FASTA is required"))?;
    let output_path = cli
        .output
        .as_ref()
        .ok_or_else(|| anyhow!("--output is required"))?;

    if output_path.exists() {
        bail!(
            "output file already exists: {}\nrefusing to overwrite — choose a different path",
            output_path.display()
        );
    }

    let template = read_template(template_path)?;
    let references = read_references(reference_path, "reference")?;
    let exclusivity = match cli.diff.as_ref() {
        Some(path) => Some(read_references(path, "exclusivity")?),
        None => None,
    };

    let params = cli.into_analysis_params(&cfg)?;

    let plan = build_length_plan(
        template.sequence.len(),
        params.min_oligo_length,
        params.max_oligo_length,
        params.length_skip,
        params.resolution,
    );

    if plan.iter().all(|(_, n)| *n == 0) {
        bail!(
            "template ({} bp) is too short for any requested oligo length ({}..={})",
            template.sequence.len(),
            params.min_oligo_length,
            params.max_oligo_length
        );
    }

    eprintln!(
        "diffalign: template={} ({} bp), references={}, {}{}using {}",
        template.name,
        template.sequence.len(),
        references.len(),
        match &exclusivity {
            Some(e) => format!("exclusivity={} sequences, ", e.len()),
            None => String::new(),
        },
        if exclusivity.is_some() {
            "differential mode, "
        } else {
            ""
        },
        format_threads(&params.thread_count),
    );

    let render_progress = !cli.quiet && io::stderr().is_terminal();
    let reporter = if render_progress {
        Reporter::new(plan)
    } else {
        Reporter::quiet()
    };

    let pause_flag = if render_progress {
        Some(PauseFlag::new())
    } else {
        None
    };

    let _key_listener = match (pause_flag.clone(), reporter.multi()) {
        (Some(pf), Some(multi)) => KeyListener::try_spawn(pf, multi),
        _ => None,
    };

    let results = run_screening(
        &template,
        &references,
        &params,
        exclusivity.as_ref(),
        reporter.sender(),
        pause_flag,
    );
    drop(_key_listener);
    reporter.finish();

    write_results(&results, output_path)?;

    eprintln!("✓ wrote {}", output_path.display());
    if timer_enabled {
        let elapsed = started.elapsed();
        eprintln!("elapsed: {}", format_duration(elapsed));
    }
    Ok(())
}

fn format_duration(d: std::time::Duration) -> String {
    let total_secs = d.as_secs();
    let millis = d.subsec_millis();
    let hours = total_secs / 3600;
    let minutes = (total_secs % 3600) / 60;
    let seconds = total_secs % 60;
    if hours > 0 {
        format!("{}h {:02}m {:02}.{:03}s", hours, minutes, seconds, millis)
    } else if minutes > 0 {
        format!("{}m {:02}.{:03}s", minutes, seconds, millis)
    } else {
        format!("{}.{:03}s", seconds, millis)
    }
}

fn cmd_mkini(path: &Path) -> Result<()> {
    config::write_default_ini(path)?;
    println!("wrote default config to {}", path.display());
    Ok(())
}

fn cmd_check_config(path: &Path) -> Result<()> {
    println!("config path: {}", path.display());
    if !path.exists() {
        println!("status: missing");
        std::process::exit(2);
    }
    match config::load(path) {
        Ok(_) => {
            println!("status: ok");
            Ok(())
        }
        Err(e) => {
            println!("status: invalid ({:#})", e);
            std::process::exit(2);
        }
    }
}

/// Load the INI, or — if it's missing and we're on a TTY — ask whether to
/// create it. Returns `Ok(None)` when the user declined and the program
/// should exit cleanly.
fn load_or_prompt(path: &Path) -> Result<Option<Config>> {
    if path.exists() {
        return config::load(path).map(Some);
    }

    eprintln!("No config file found at {}", path.display());

    if !io::stdin().is_terminal() {
        bail!(
            "no diffalign.ini found at {}. Run `diffalign --mkini` to create one.",
            path.display()
        );
    }

    eprint!("Create one with default values? [y/N]: ");
    io::stderr().flush().ok();

    let mut answer = String::new();
    io::stdin()
        .read_line(&mut answer)
        .context("reading confirmation")?;
    let answer = answer.trim().to_ascii_lowercase();
    if !(answer == "y" || answer == "yes") {
        eprintln!("aborting; run `diffalign --mkini` to create one later.");
        return Ok(None);
    }

    config::write_default_ini(path)?;
    eprintln!("created default config at {}", path.display());
    config::load(path).map(Some)
}

fn read_template(path: &PathBuf) -> Result<diffalign::analysis::TemplateData> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("reading template file {}", path.display()))?;
    parse_template_fasta(&text)
        .map_err(|e| anyhow!("parsing template {}: {}", path.display(), e))
}

fn read_references(path: &PathBuf, kind: &str) -> Result<ReferenceData> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("reading {} file {}", kind, path.display()))?;
    parse_reference_fasta(&text)
        .map_err(|e| anyhow!("parsing {} {}: {}", kind, path.display(), e))
}

fn write_results(results: &ScreeningResults, path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() && !parent.exists() {
            fs::create_dir_all(parent)
                .with_context(|| format!("creating output directory {}", parent.display()))?;
        }
    }
    let json = serde_json::to_string_pretty(results).context("serializing results to JSON")?;
    fs::write(path, json).with_context(|| format!("writing results to {}", path.display()))?;
    Ok(())
}

fn format_threads(tc: &diffalign::analysis::ThreadCount) -> String {
    match tc {
        diffalign::analysis::ThreadCount::Auto => "auto threads".to_string(),
        diffalign::analysis::ThreadCount::Fixed(n) => format!("{} thread(s)", n),
    }
}

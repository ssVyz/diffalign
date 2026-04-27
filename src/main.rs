use std::fs;
use std::io::{self, IsTerminal, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};
use clap::Parser;
use mimalloc::MiMalloc;

use diffalign::analysis::{
    ReferenceData, ScreeningResults, parse_reference_fasta, parse_template_fasta, run_screening,
};
use diffalign::cli::Cli;
use diffalign::config::{self, Config};
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

    let render_progress = !cli.quiet && io::stderr().is_terminal();
    let reporter = if render_progress {
        Reporter::new(plan)
    } else {
        Reporter::quiet()
    };

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

    let results = run_screening(
        &template,
        &references,
        &params,
        exclusivity.as_ref(),
        reporter.sender(),
    );
    reporter.finish();

    write_results(&results, output_path)?;

    eprintln!("✓ wrote {}", output_path.display());
    Ok(())
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

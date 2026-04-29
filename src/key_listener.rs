//! Background thread that watches for hotkeys (pause/resume, abort) while
//! the screening run is in progress.
//!
//! Puts the terminal into raw mode so single keypresses are delivered without
//! waiting for a newline. A Drop guard on the listener restores cooked mode.
//!
//! Set `DIFFALIGN_DEBUG_KEYS=1` to print every received key event for
//! diagnosing terminal/raw-mode issues.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread::{self, JoinHandle};
use std::time::Duration;

use crossterm::event::{self, Event, KeyCode, KeyEvent, KeyEventKind, KeyModifiers};
use crossterm::terminal::{disable_raw_mode, enable_raw_mode};
use indicatif::MultiProgress;

use crate::pause::PauseFlag;

pub struct KeyListener {
    shutdown: Arc<AtomicBool>,
    handle: Option<JoinHandle<()>>,
}

impl KeyListener {
    /// Enable raw mode and spawn the listener thread. Returns `None` if raw
    /// mode could not be enabled — a warning is printed so the user knows
    /// pause/resume is unavailable rather than silently broken.
    pub fn try_spawn(pause: PauseFlag, multi: MultiProgress) -> Option<Self> {
        if let Err(e) = enable_raw_mode() {
            let _ = multi.println(format!(
                "[pause/resume unavailable: could not enable terminal raw mode: {}]",
                e
            ));
            return None;
        }

        let shutdown = Arc::new(AtomicBool::new(false));
        let shutdown_clone = Arc::clone(&shutdown);

        let _ = multi.println("[press P or Space to pause, Ctrl+C to abort]");

        let handle = thread::spawn(move || {
            run(pause, multi, shutdown_clone);
        });

        Some(Self {
            shutdown,
            handle: Some(handle),
        })
    }
}

impl Drop for KeyListener {
    fn drop(&mut self) {
        self.shutdown.store(true, Ordering::Relaxed);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
        let _ = disable_raw_mode();
    }
}

fn debug_keys_enabled() -> bool {
    std::env::var_os("DIFFALIGN_DEBUG_KEYS")
        .map(|v| v != "0" && !v.is_empty())
        .unwrap_or(false)
}

fn run(pause: PauseFlag, multi: MultiProgress, shutdown: Arc<AtomicBool>) {
    let debug = debug_keys_enabled();
    if debug {
        let _ = multi.println("[debug-keys: listener thread started]");
    }

    while !shutdown.load(Ordering::Relaxed) {
        match event::poll(Duration::from_millis(200)) {
            Ok(true) => match event::read() {
                Ok(Event::Key(k)) => {
                    if debug {
                        let _ = multi.println(format!("[debug-keys: {:?}]", k));
                    }
                    handle_key(&pause, &multi, k);
                }
                Ok(other) => {
                    if debug {
                        let _ = multi.println(format!("[debug-keys: non-key event {:?}]", other));
                    }
                }
                Err(e) => {
                    let _ = multi.println(format!("[key listener: read error: {}]", e));
                }
            },
            Ok(false) => {}
            Err(e) => {
                let _ = multi.println(format!("[key listener: poll error: {}]", e));
                // Avoid a hot error loop if poll consistently fails.
                thread::sleep(Duration::from_millis(500));
            }
        }
    }
}

fn handle_key(pause: &PauseFlag, multi: &MultiProgress, k: KeyEvent) {
    // crossterm fires Press and (where supported) Release/Repeat. Only act on
    // key-down events so a single tap doesn't toggle twice.
    if !matches!(k.kind, KeyEventKind::Press | KeyEventKind::Repeat) {
        return;
    }

    if k.modifiers.contains(KeyModifiers::CONTROL)
        && matches!(k.code, KeyCode::Char('c') | KeyCode::Char('C'))
    {
        let _ = disable_raw_mode();
        eprintln!("\n[aborted]");
        std::process::exit(130);
    }

    match k.code {
        KeyCode::Char('p') | KeyCode::Char('P') | KeyCode::Char(' ') => {
            let now_paused = pause.toggle();
            let _ = multi.println(if now_paused {
                "[paused — press P or Space to resume]"
            } else {
                "[resumed]"
            });
        }
        _ => {}
    }
}

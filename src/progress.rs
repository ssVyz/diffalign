//! Progress reporting backed by `indicatif`.
//!
//! Renders two stacked bars on stderr (lengths + positions within the current
//! length). Auto-degrades to nothing when stderr is not a TTY (indicatif
//! handles that), and can be disabled outright via `--quiet`.

use std::sync::mpsc::{Receiver, Sender, channel};
use std::thread::{self, JoinHandle};
use std::time::Duration;

use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

use crate::analysis::ProgressUpdate;

/// (oligo_length, total_positions_for_that_length)
pub type LengthPlan = Vec<(u32, usize)>;

/// Compute how many positions each oligo length will analyze, given
/// the template length and resolution. Skipped lengths are pre-filtered out.
pub fn build_length_plan(
    template_len: usize,
    min_len: u32,
    max_len: u32,
    length_skip: u32,
    resolution: u32,
) -> LengthPlan {
    let step = length_skip.saturating_add(1).max(1) as usize;
    let res = resolution.max(1) as usize;
    let mut plan = Vec::new();
    if min_len > max_len {
        return plan;
    }
    for length in (min_len..=max_len).step_by(step) {
        let length_usize = length as usize;
        if template_len < length_usize {
            plan.push((length, 0));
            continue;
        }
        let max_start = template_len - length_usize;
        let total_positions = (0..=max_start).step_by(res).count();
        plan.push((length, total_positions));
    }
    plan
}

/// Owns the indicatif bars and the worker thread that drains
/// `ProgressUpdate`s into them.
pub struct Reporter {
    sender: Option<Sender<ProgressUpdate>>,
    handle: Option<JoinHandle<()>>,
    multi: Option<MultiProgress>,
}

impl Reporter {
    pub fn quiet() -> Self {
        Self {
            sender: None,
            handle: None,
            multi: None,
        }
    }

    pub fn new(plan: LengthPlan) -> Self {
        let (tx, rx) = channel::<ProgressUpdate>();
        let multi = MultiProgress::new();

        let length_bar = multi.add(ProgressBar::new(plan.len() as u64));
        length_bar.set_style(
            ProgressStyle::with_template(
                "{prefix:>9} [{bar:30.cyan/blue}] {pos}/{len}  {msg}",
            )
            .unwrap_or(ProgressStyle::default_bar())
            .progress_chars("=>-"),
        );
        length_bar.set_prefix("Lengths");

        let initial_total = plan
            .iter()
            .find(|(_, n)| *n > 0)
            .map(|(_, n)| *n as u64)
            .unwrap_or(0);
        let position_bar = multi.add(ProgressBar::new(initial_total));
        position_bar.set_style(
            ProgressStyle::with_template(
                "{prefix:>9} [{bar:30.green/blue}] {pos}/{len}  {msg}",
            )
            .unwrap_or(ProgressStyle::default_bar())
            .progress_chars("=>-"),
        );
        position_bar.set_prefix("Position");
        position_bar.enable_steady_tick(Duration::from_millis(120));

        let handle = thread::spawn(move || {
            run(rx, plan, length_bar, position_bar);
        });

        Self {
            sender: Some(tx),
            handle: Some(handle),
            multi: Some(multi),
        }
    }

    /// Hand out a sender for the screening engine. Returns `None` in quiet
    /// mode so the screener skips its progress emission entirely.
    pub fn sender(&self) -> Option<Sender<ProgressUpdate>> {
        self.sender.clone()
    }

    /// Borrow a clone of the underlying `MultiProgress` so other components
    /// (e.g. the key listener) can emit status lines that play nicely with
    /// the live bars.
    pub fn multi(&self) -> Option<MultiProgress> {
        self.multi.clone()
    }

    /// Drop the sender and wait for the worker thread to finalize the bars.
    pub fn finish(mut self) {
        self.sender.take();
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
        if let Some(multi) = self.multi.take() {
            let _ = multi.clear();
        }
    }
}

fn run(
    rx: Receiver<ProgressUpdate>,
    plan: LengthPlan,
    length_bar: ProgressBar,
    position_bar: ProgressBar,
) {
    let mut current_length: Option<u32> = None;

    while let Ok(update) = rx.recv() {
        // Re-anchor the position bar when the length changes.
        if current_length != Some(update.current_length) {
            current_length = Some(update.current_length);
            position_bar.reset();
            let total = plan
                .iter()
                .find(|(l, _)| *l == update.current_length)
                .map(|(_, t)| *t as u64)
                .unwrap_or(update.total_positions as u64);
            position_bar.set_length(total);
            position_bar.set_message(format!("length {}", update.current_length));
            length_bar.set_position(update.lengths_completed as u64);
            length_bar.set_message(format!(
                "current: {} bp",
                update.current_length
            ));
        }

        // The screener emits "completed positions" via the progress message;
        // the position field is the most recent processed position.
        // We approximate completed by (lengths_completed * 0) + emitted index;
        // the message-derived count already reflects this — read from the
        // engine's own progress field.
        let completed_for_length =
            extract_completed(&update.message).unwrap_or(update.current_position + 1);
        position_bar.set_position(completed_for_length as u64);

        // If a length finished, bump the length bar.
        let target_total = plan
            .iter()
            .find(|(l, _)| *l == update.current_length)
            .map(|(_, t)| *t)
            .unwrap_or(0);
        if target_total > 0 && completed_for_length >= target_total {
            length_bar.inc(1);
            current_length = None; // force re-anchor on the next update
        }
    }

    position_bar.finish_and_clear();
    length_bar.finish_and_clear();
}

/// The screening engine encodes the completed position count in
/// `ProgressUpdate.message` (e.g. "Length 1/3: Position 340/2912"). Pull it
/// out so we can drive the position bar accurately even though the field
/// `current_position` is the *index* of the most recent position rather than
/// a count.
fn extract_completed(message: &str) -> Option<usize> {
    let after_position = message.split("Position ").nth(1)?;
    let head = after_position.split('/').next()?;
    head.trim().parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plan_handles_skip() {
        let plan = build_length_plan(50, 18, 25, 1, 1);
        assert_eq!(plan.iter().map(|(l, _)| *l).collect::<Vec<_>>(), vec![18, 20, 22, 24]);
        // template_len 50, length 18 → positions 0..=32 = 33
        assert_eq!(plan[0].1, 33);
    }

    #[test]
    fn plan_handles_resolution() {
        let plan = build_length_plan(50, 20, 20, 0, 5);
        assert_eq!(plan.len(), 1);
        // positions 0,5,10,15,20,25,30 = 7
        assert_eq!(plan[0].1, 7);
    }

    #[test]
    fn plan_handles_template_too_short() {
        let plan = build_length_plan(10, 18, 18, 0, 1);
        assert_eq!(plan, vec![(18, 0)]);
    }

    #[test]
    fn extract_completed_parses_message() {
        assert_eq!(
            extract_completed("Length 1/3: Position 340/2912"),
            Some(340)
        );
        assert_eq!(extract_completed("garbage"), None);
    }
}

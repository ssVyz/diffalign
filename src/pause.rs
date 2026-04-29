//! Cooperative pause/resume primitive for parallel workers.
//!
//! Workers call `wait_if_paused` at task boundaries (in our case, the top of
//! the per-position rayon closure). The fast path is a single relaxed atomic
//! load — no mutex, no contention. When toggled, parked workers sleep on a
//! Condvar until resumed, so paused work consumes zero CPU.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex};

#[derive(Clone)]
pub struct PauseFlag {
    inner: Arc<Inner>,
}

struct Inner {
    paused: AtomicBool,
    lock: Mutex<()>,
    condvar: Condvar,
}

impl PauseFlag {
    pub fn new() -> Self {
        Self {
            inner: Arc::new(Inner {
                paused: AtomicBool::new(false),
                lock: Mutex::new(()),
                condvar: Condvar::new(),
            }),
        }
    }

    /// Block while the pause flag is set. Cheap when not paused.
    pub fn wait_if_paused(&self) {
        if !self.inner.paused.load(Ordering::Acquire) {
            return;
        }
        let mut guard = self.inner.lock.lock().unwrap();
        while self.inner.paused.load(Ordering::Acquire) {
            guard = self.inner.condvar.wait(guard).unwrap();
        }
    }

    /// Flip pause state. Returns the new state (true = now paused).
    pub fn toggle(&self) -> bool {
        let _g = self.inner.lock.lock().unwrap();
        let new_state = !self.inner.paused.load(Ordering::Acquire);
        self.inner.paused.store(new_state, Ordering::Release);
        if !new_state {
            self.inner.condvar.notify_all();
        }
        new_state
    }

    pub fn is_paused(&self) -> bool {
        self.inner.paused.load(Ordering::Acquire)
    }
}

impl Default for PauseFlag {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;
    use std::thread;
    use std::time::Duration;

    #[test]
    fn fast_path_when_not_paused() {
        let p = PauseFlag::new();
        for _ in 0..1000 {
            p.wait_if_paused();
        }
    }

    #[test]
    fn toggle_pauses_and_resumes() {
        let p = PauseFlag::new();
        assert!(!p.is_paused());
        assert!(p.toggle());
        assert!(p.is_paused());
        assert!(!p.toggle());
        assert!(!p.is_paused());
    }

    #[test]
    fn worker_blocks_then_resumes() {
        let p = PauseFlag::new();
        p.toggle();
        let counter = Arc::new(AtomicUsize::new(0));

        let pc = p.clone();
        let cc = Arc::clone(&counter);
        let h = thread::spawn(move || {
            pc.wait_if_paused();
            cc.fetch_add(1, Ordering::Relaxed);
        });

        thread::sleep(Duration::from_millis(50));
        assert_eq!(counter.load(Ordering::Relaxed), 0);

        p.toggle();
        h.join().unwrap();
        assert_eq!(counter.load(Ordering::Relaxed), 1);
    }
}

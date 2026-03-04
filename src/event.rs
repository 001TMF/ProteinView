use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;
use std::sync::Arc;
use std::thread;
use std::time::Duration;

use crossterm::event::{self, Event, KeyEvent};

/// Spawns a dedicated input thread that sends key events through a channel.
/// The thread checks `quit_flag` each iteration and stops when it's set.
/// Returns (receiver, quit_flag).
pub fn spawn_input_thread() -> (mpsc::Receiver<KeyEvent>, Arc<AtomicBool>) {
    let (tx, rx) = mpsc::channel();
    let quit_flag = Arc::new(AtomicBool::new(false));
    let quit = quit_flag.clone();

    thread::spawn(move || {
        loop {
            if quit.load(Ordering::Relaxed) {
                break;
            }
            // Poll with short timeout so we can check quit_flag regularly
            if event::poll(Duration::from_millis(10)).unwrap_or(false) {
                if let Ok(Event::Key(key)) = event::read() {
                    if tx.send(key).is_err() {
                        break; // receiver dropped
                    }
                }
            }
        }
    });

    (rx, quit_flag)
}

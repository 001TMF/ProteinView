use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

use crossterm::event::{self, Event, KeyEvent};

/// Application-level event that wraps terminal events we care about.
pub enum AppEvent {
    /// A keyboard input event.
    Key(KeyEvent),
    /// The terminal was resized to (columns, rows).
    Resize(u16, u16),
}

/// Spawns a dedicated input thread that sends key and resize events through a channel.
/// The thread checks `quit_flag` each iteration and stops when it's set.
/// Returns (receiver, quit_flag).
pub fn spawn_input_thread() -> (mpsc::Receiver<AppEvent>, Arc<AtomicBool>) {
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
                match event::read() {
                    Ok(Event::Key(key)) => {
                        if tx.send(AppEvent::Key(key)).is_err() {
                            break; // receiver dropped
                        }
                    }
                    Ok(Event::Resize(w, h)) => {
                        if tx.send(AppEvent::Resize(w, h)).is_err() {
                            break; // receiver dropped
                        }
                    }
                    _ => {}
                }
            }
        }
    });

    (rx, quit_flag)
}

use std::time::Duration;
use crossterm::event::{self, Event, KeyEvent};
use anyhow::Result;

/// Poll for keyboard events with a timeout
pub fn poll_event(timeout: Duration) -> Result<Option<KeyEvent>> {
    if event::poll(timeout)? {
        if let Event::Key(key) = event::read()? {
            return Ok(Some(key));
        }
    }
    Ok(None)
}

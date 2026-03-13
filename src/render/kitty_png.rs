//! Compressed Kitty graphics protocol transmitter using PNG encoding.
//!
//! ratatui-image sends Kitty images as raw RGBA32 base64-encoded (`f=32`),
//! which is ~1.3MB per frame for a 640x384 render.  This module replaces
//! that path with PNG encoding (`f=100`) which compresses protein renders
//! (mostly black/transparent background) to roughly 5-15% of raw size.

use std::fmt::Write;
use std::io::Cursor;
use std::sync::atomic::{AtomicU32, Ordering};

use base64::Engine;
use image::DynamicImage;
use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::widgets::Widget;

/// Monotonically increasing image ID so Kitty can track/clean up images.
static NEXT_ID: AtomicU32 = AtomicU32::new(1);

/// A ratatui `Widget` that renders a `DynamicImage` via the Kitty graphics
/// protocol using PNG compression (`f=100`).
pub struct KittyPngImage {
    transmit: String,
    area: Rect,
}

impl KittyPngImage {
    /// Create a new compressed Kitty image widget.
    ///
    /// The image is immediately PNG-encoded and base64-chunked into the
    /// Kitty escape sequence.  Call this outside the draw closure if you
    /// want to time encoding separately.
    pub fn new(img: &DynamicImage, area: Rect) -> Self {
        let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
        let (w, h) = (img.width(), img.height());

        // Encode as PNG.
        let mut png_bytes: Vec<u8> = Vec::new();
        img.write_to(&mut Cursor::new(&mut png_bytes), image::ImageFormat::Png)
            .expect("PNG encoding failed");

        // Base64.
        let b64 = base64::engine::general_purpose::STANDARD.encode(&png_bytes);

        // Build Kitty escape, chunked to <=4096 base64 chars.
        const CHARS_PER_CHUNK: usize = 4096;
        let chunks: Vec<&str> = b64
            .as_bytes()
            .chunks(CHARS_PER_CHUNK)
            .map(|c| std::str::from_utf8(c).unwrap())
            .collect();
        let chunk_count = chunks.len();

        let mut transmit = String::with_capacity(b64.len() + chunk_count * 80);

        for (i, chunk) in chunks.iter().enumerate() {
            write!(transmit, "\x1b_Gq=2,").unwrap();
            if i == 0 {
                // a=T = transmit+display, f=100 = PNG, t=d = direct (inline)
                // C=1 = do not move cursor, z=-1 = put behind text
                write!(transmit, "i={id},a=T,f=100,t=d,C=1,s={w},v={h},").unwrap();
            }
            let more = u8::from(chunk_count > (i + 1));
            write!(transmit, "m={more};{chunk}\x1b\\").unwrap();
        }

        Self { transmit, area }
    }
}

impl Widget for KittyPngImage {
    fn render(self, area: Rect, buf: &mut Buffer) {
        // Write the entire Kitty escape into the first cell.
        // Subsequent cells are marked as skip so ratatui doesn't
        // overwrite the image with blanks.
        if let Some(cell) = buf.cell_mut((area.left(), area.top())) {
            cell.set_symbol(&self.transmit);
        }

        // Mark all cells in the area as skip so ratatui's diff engine
        // doesn't emit escape sequences that would overwrite the image.
        for y in 0..area.height.min(self.area.height) {
            for x in 0..area.width.min(self.area.width) {
                // Skip the first cell (0,0) since we wrote our symbol there.
                if y == 0 && x == 0 {
                    continue;
                }
                if let Some(cell) = buf.cell_mut((area.left() + x, area.top() + y)) {
                    cell.set_skip(true);
                }
            }
        }
    }
}

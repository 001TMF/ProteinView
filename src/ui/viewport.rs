use image::DynamicImage;
use ratatui::layout::Rect;
use ratatui::Frame;
use ratatui_image::picker::ProtocolType;
use ratatui_image::{Image, Resize};

use crate::app::App;
use crate::render::braille;
use crate::render::framebuffer::framebuffer_to_braille_widget;
use crate::render::hd;

/// Render the main 3D viewport
pub fn render_viewport(frame: &mut Frame, area: Rect, app: &App) {
    if app.hd_mode {
        render_hd_viewport(frame, area, app);
    } else {
        // Normal mode uses Braille (2x4 dots per cell, higher resolution but monochrome per cell)
        let width = area.width as f64 * 2.0;
        let height = area.height as f64 * 4.0;

        let canvas = braille::render_protein(
            &app.protein,
            &app.camera,
            &app.color_scheme,
            app.viz_mode,
            width,
            height,
        );

        frame.render_widget(canvas, area);
    }
}

/// Render the HD viewport using graphics protocol (Sixel/Kitty/iTerm2) when
/// available, falling back to half-block characters otherwise.
fn render_hd_viewport(frame: &mut Frame, area: Rect, app: &App) {
    let proto = app.picker.protocol_type();
    let (font_w, font_h) = app.picker.font_size();

    // Determine framebuffer pixel dimensions.
    // With a true graphics protocol we render at full pixel resolution
    // (cols * font_width, rows * font_height).  For the colored braille
    // fallback we render at braille resolution: cols*2 wide, rows*4 tall
    // (2 dot columns and 4 dot rows per terminal cell).
    let (px_w, px_h) = if proto != ProtocolType::Halfblocks && font_w > 0 && font_h > 0 {
        (
            area.width as f64 * font_w as f64,
            area.height as f64 * font_h as f64,
        )
    } else {
        (
            area.width as f64 * 2.0,
            area.height as f64 * 4.0,
        )
    };

    // Rasterize the 3D scene into our software framebuffer.
    let fb = hd::render_hd_framebuffer(
        &app.protein,
        &app.camera,
        &app.color_scheme,
        app.viz_mode,
        px_w,
        px_h,
        &app.mesh_cache,
    );

    // If the terminal supports a real graphics protocol, convert the
    // framebuffer to an image and render it through ratatui-image.
    if proto != ProtocolType::Halfblocks {
        let rgb_img = fb.to_rgb_image();
        let dyn_img = DynamicImage::ImageRgb8(rgb_img);
        match app.picker.new_protocol(dyn_img, area, Resize::Fit(None)) {
            Ok(protocol) => {
                let widget = Image::new(&protocol);
                frame.render_widget(widget, area);
                return;
            }
            Err(_) => {
                // Fall through to half-block rendering on error.
            }
        }
    }

    // Fallback: colored braille character rendering (always works).
    // Each terminal cell maps to a 2x4 block of framebuffer pixels,
    // giving 4x the spatial resolution of half-blocks at the cost of
    // per-cell (rather than per-pixel) coloring.
    let widget = framebuffer_to_braille_widget(&fb);
    frame.render_widget(widget, area);
}

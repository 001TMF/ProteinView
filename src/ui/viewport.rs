use image::DynamicImage;
use ratatui::layout::Rect;
use ratatui::Frame;
use ratatui_image::picker::ProtocolType;
use ratatui_image::{Image, Resize};

use crate::app::{App, RenderMode};
use crate::render::braille;
use crate::render::framebuffer::framebuffer_to_braille_widget;
use crate::render::hd;
use crate::render::kitty_png::KittyPngImage;

/// Render the main 3D viewport
pub fn render_viewport(frame: &mut Frame, area: Rect, app: &App) {
    match app.render_mode {
        RenderMode::Braille => {
            // Braille mode: 2x4 dots per cell, higher resolution but monochrome per cell
            let width = area.width as f64 * 2.0;
            let height = area.height as f64 * 4.0;

            let canvas = braille::render_protein(
                &app.protein,
                &app.camera,
                &app.color_scheme,
                app.viz_mode,
                width,
                height,
                app.show_ligands,
            );

            frame.render_widget(canvas, area);
        }
        RenderMode::HalfBlock => {
            // HalfBlock mode: render at braille resolution (2x4 per cell) through
            // the HD rasterizer (Lambert shading, z-buffer, depth fog) and convert
            // to colored braille characters.  This gives the same spatial resolution
            // as the basic Braille renderer but with much higher quality shading.
            let width = area.width as f64 * 2.0;
            let height = area.height as f64 * 4.0;

            let fb = hd::render_hd_framebuffer(
                &app.protein,
                &app.camera,
                &app.color_scheme,
                app.viz_mode,
                width,
                height,
                &app.mesh_cache,
                app.show_ligands,
            );

            let widget = framebuffer_to_braille_widget(&fb);
            frame.render_widget(widget, area);
        }
        RenderMode::FullHD => {
            render_fullhd_viewport(frame, area, app);
        }
    }
}

/// Render the FullHD viewport using graphics protocol (Sixel/Kitty/iTerm2) when
/// available, falling back to colored braille characters otherwise.
fn render_fullhd_viewport(frame: &mut Frame, area: Rect, app: &App) {
    let proto = app.picker.protocol_type();
    let (font_w, font_h) = app.picker.font_size();

    // Determine framebuffer pixel dimensions.
    // With a true graphics protocol we render at full pixel resolution
    // (cols * font_width, rows * font_height).  For the colored braille
    // fallback we render at braille resolution: cols*2 wide, rows*4 tall.
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

    // Rasterize the 3D scene into a software framebuffer.
    let fb = hd::render_hd_framebuffer(
        &app.protein,
        &app.camera,
        &app.color_scheme,
        app.viz_mode,
        px_w,
        px_h,
        &app.mesh_cache,
        app.show_ligands,
    );

    // If the terminal supports a real graphics protocol, convert the
    // framebuffer to an image and send it.
    if proto != ProtocolType::Halfblocks {
        if proto == ProtocolType::Kitty {
            // Use our custom PNG-compressed Kitty transmitter.
            // This is ~10-20x smaller than ratatui-image's raw RGBA path,
            // making FullHD viable over SSH.
            let dyn_img = DynamicImage::ImageRgba8(fb.to_rgba_image());
            let widget = KittyPngImage::new(&dyn_img, area);
            frame.render_widget(widget, area);
            return;
        }

        // Sixel/iTerm2: use ratatui-image (no PNG option for Sixel).
        let dyn_img = DynamicImage::ImageRgb8(fb.to_rgb_image());
        match app.picker.new_protocol(dyn_img, area, Resize::Fit(None)) {
            Ok(protocol) => {
                let widget = Image::new(&protocol);
                frame.render_widget(widget, area);
                return;
            }
            Err(_) => {
                // Fall through to braille rendering on error.
            }
        }
    }

    // Fallback: colored braille character rendering (always works).
    let widget = framebuffer_to_braille_widget(&fb);
    frame.render_widget(widget, area);
}

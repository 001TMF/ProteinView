use image::DynamicImage;
use ratatui::layout::Rect;
use ratatui::Frame;
use ratatui_image::picker::ProtocolType;
use ratatui_image::{Image, Resize};

use crate::app::{App, RenderMode};
use crate::render::braille;
use crate::render::framebuffer::{framebuffer_to_braille_widget, framebuffer_to_widget};
use crate::render::hd;

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
            // HalfBlock mode: 1 pixel per column, 2 pixels per row
            // Rasterize inline and convert to half-block Paragraph widget
            let width = area.width as f64 * 1.0;
            let height = area.height as f64 * 2.0;

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

            let widget = framebuffer_to_widget(&fb);
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
    // framebuffer to an image and render it through ratatui-image.
    // Kitty supports RGBA transparency; Sixel/iTerm2 do not, so use
    // opaque RGB to avoid rendering artifacts on those protocols.
    if proto != ProtocolType::Halfblocks {
        let dyn_img = if proto == ProtocolType::Kitty {
            DynamicImage::ImageRgba8(fb.to_rgba_image())
        } else {
            DynamicImage::ImageRgb8(fb.to_rgb_image())
        };
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

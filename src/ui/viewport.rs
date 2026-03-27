use image::DynamicImage;
use ratatui::layout::Rect;
use ratatui::Frame;
use ratatui_image::picker::ProtocolType;
use ratatui_image::{Image, Resize};

use crate::app::{App, RenderMode};
use crate::model::interface::Interaction;
use crate::render::braille;
use crate::render::framebuffer::framebuffer_to_braille_widget;
use crate::render::hd;
use crate::render::kitty_png::KittyPngImage;

/// Render the main 3D viewport
pub fn render_viewport(frame: &mut Frame, area: Rect, app: &App) {
    let interactions: &[Interaction] = if app.show_interface && app.show_interactions {
        &app.interface_analysis.interactions
    } else {
        &[]
    };

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
                interactions,
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
                interactions,
            );

            let widget = framebuffer_to_braille_widget(&fb);
            frame.render_widget(widget, area);
        }
        RenderMode::FullHD => {
            render_fullhd_viewport(frame, area, app, interactions);
        }
    }
}

/// Render the FullHD viewport using graphics protocol (Sixel/Kitty/iTerm2) when
/// available, falling back to colored braille characters otherwise.
fn render_fullhd_viewport(frame: &mut Frame, area: Rect, app: &App, interactions: &[Interaction]) {
    let proto = app.picker.protocol_type();
    let (font_w, font_h) = app.picker.font_size();

    // Determine framebuffer pixel dimensions.
    // With a true graphics protocol we render at full pixel resolution
    // (cols * font_width, rows * font_height).  For the colored braille
    // fallback we render at braille resolution: cols*2 wide, rows*4 tall.
    //
    // During interaction (auto-rotate), render at half resolution for the
    // graphics-protocol path. The terminal upscales via Kitty `c=/r=` params.
    // Even with parallel rasterization, half-res keeps frame rates smooth
    // on large structures.
    let is_graphics = proto != ProtocolType::Halfblocks && font_w > 0 && font_h > 0;
    let is_large = app.protein.residue_count() > crate::app::LARGE_STRUCTURE_THRESHOLD;
    let scale = if is_graphics && is_large && app.is_interacting() { 0.5 } else { 1.0 };
    let (px_w, px_h) = if is_graphics {
        (
            area.width as f64 * font_w as f64 * scale,
            area.height as f64 * font_h as f64 * scale,
        )
    } else {
        (
            area.width as f64 * 2.0,
            area.height as f64 * 4.0,
        )
    };

    // Rasterize the 3D scene into a software framebuffer.
    // When rendering at reduced resolution, scale camera zoom to match so the
    // protein fills the same relative area of the smaller buffer. Kitty's
    // c=/r= params then upscale the result to fill the full viewport.
    let mut cam = app.camera.clone();
    if scale < 1.0 {
        cam.zoom *= scale;
        cam.pan_x *= scale;
        cam.pan_y *= scale;
    }
    let fb = hd::render_hd_framebuffer(
        &app.protein,
        &cam,
        &app.color_scheme,
        app.viz_mode,
        px_w,
        px_h,
        &app.mesh_cache,
        app.show_ligands,
        interactions,
    );

    // If the terminal supports a real graphics protocol, convert the
    // framebuffer to an image and send it.
    if proto != ProtocolType::Halfblocks {
        if proto == ProtocolType::Kitty {
            // Use our custom PNG-compressed Kitty transmitter.
            // This is ~10-20x smaller than ratatui-image's raw RGBA path,
            // making FullHD viable over SSH.
            let dyn_img = DynamicImage::ImageRgba8(fb.to_rgba_image());
            if let Some(widget) = KittyPngImage::new(&dyn_img, area) {
                frame.render_widget(widget, area);
                return;
            }
            // PNG encoding failed — fall through to braille.
        }

        // Sixel/iTerm2: use ratatui-image (no PNG option for Sixel).
        if proto != ProtocolType::Kitty {
            let dyn_img = DynamicImage::ImageRgb8(fb.to_rgb_image());
            if let Ok(protocol) = app.picker.new_protocol(dyn_img, area, Resize::Fit(None)) {
                let widget = Image::new(&protocol);
                frame.render_widget(widget, area);
                return;
            }
            // Protocol error — fall through to braille.
        }
    }

    // Fallback: colored braille character rendering (always works).
    let widget = framebuffer_to_braille_widget(&fb);
    frame.render_widget(widget, area);
}

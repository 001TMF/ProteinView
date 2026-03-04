use ratatui::layout::Rect;
use ratatui::Frame;

use crate::app::App;
use crate::render::braille;
use crate::render::cartoon;

/// Render the main 3D viewport
pub fn render_viewport(frame: &mut Frame, area: Rect, app: &App) {
    if app.hd_mode {
        // HD / Cartoon mode uses HalfBlock (1x2 colored pixels per cell)
        let width = area.width as f64 * 1.0;
        let height = area.height as f64 * 2.0;

        let canvas = cartoon::render_cartoon(
            &app.protein,
            &app.camera,
            &app.color_scheme,
            app.viz_mode,
            width,
            height,
        );

        frame.render_widget(canvas, area);
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

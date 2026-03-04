use ratatui::layout::Rect;
use ratatui::Frame;

use crate::app::App;
use crate::render::braille;

/// Render the main 3D viewport
pub fn render_viewport(frame: &mut Frame, area: Rect, app: &App) {
    let width = area.width as f64 * 2.0;  // braille gives 2x horizontal resolution
    let height = area.height as f64 * 4.0; // braille gives 4x vertical resolution

    let canvas = braille::render_protein(
        &app.protein,
        &app.camera,
        &app.color_scheme,
        width,
        height,
    );

    frame.render_widget(canvas, area);
}

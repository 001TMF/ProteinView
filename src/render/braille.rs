use ratatui::widgets::canvas::{Canvas, Context, Line};

use crate::model::protein::Protein;
use crate::render::camera::Camera;
use crate::render::color::ColorScheme;

/// Render protein as braille characters on a ratatui Canvas
pub fn render_protein<'a>(
    protein: &'a Protein,
    camera: &'a Camera,
    color_scheme: &'a ColorScheme,
    width: f64,
    height: f64,
) -> Canvas<'a, impl Fn(&mut Context<'_>) + 'a> {
    Canvas::default()
        .marker(ratatui::symbols::Marker::Braille)
        .x_bounds([-width / 2.0, width / 2.0])
        .y_bounds([-height / 2.0, height / 2.0])
        .paint(move |ctx| {
            let backbone = protein.backbone_atoms();
            if backbone.is_empty() { return; }

            // Draw lines between consecutive C-alpha atoms in each chain
            let mut prev: Option<(f64, f64, &str)> = None;
            let mut prev_chain_id = "";

            for (atom, residue, chain) in &backbone {
                let proj = camera.project(atom.x, atom.y, atom.z);

                if chain.id == prev_chain_id {
                    if let Some((px, py, _)) = prev {
                        let color = color_scheme.residue_color(residue, chain);
                        ctx.draw(&Line {
                            x1: px,
                            y1: py,
                            x2: proj.x,
                            y2: proj.y,
                            color,
                        });
                    }
                }

                prev = Some((proj.x, proj.y, &chain.id));
                prev_chain_id = &chain.id;
            }
        })
}

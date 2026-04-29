pub mod bond;
pub mod braille;
pub mod camera;
pub mod color;
pub mod framebuffer;
pub mod hd;
pub mod kitty_png;
pub mod ribbon;

use crate::app::VizMode;
use crate::model::protein::Protein;
use crate::render::camera::Camera;
use crate::render::color::ColorScheme;
use crate::render::framebuffer::Framebuffer;
use crate::render::ribbon::generate_ribbon_mesh;

// NOTE: Refactor decision — FRESH PATH (not extract).
//
// `hd::render_hd_framebuffer` already implements the full rasterization
// pipeline (cartoon tiled, backbone, wireframe, ligands, depth-fog).
// Rather than extracting orchestration from App (which is interleaved with
// TUI state: picker, Sixel/Kitty encoding, frame-skip logic, auto-rotate),
// we provide `draw_protein` as a thin wrapper that calls
// `render_hd_framebuffer` directly.  The TUI path remains unchanged;
// `batch::render_frame` uses this wrapper.  The only shared state is the
// ribbon mesh, which we build inline here (no App needed).

/// Render a protein into a [`Framebuffer`] at the requested pixel dimensions.
///
/// This is the canonical headless rasterization entry-point shared by both
/// the TUI (via `hd::render_hd_framebuffer`) and the batch pipeline.
///
/// Hotspot highlights are not yet supported in this path; pass the protein
/// with pre-applied color adjustments if needed (wired up in Task 6).
pub fn draw_protein(
    protein: &Protein,
    camera: &Camera,
    width: usize,
    height: usize,
    viz_mode: VizMode,
    color_scheme: &ColorScheme,
) -> anyhow::Result<Framebuffer> {
    let mesh = if viz_mode == VizMode::Cartoon {
        generate_ribbon_mesh(protein, color_scheme)
    } else {
        Vec::new()
    };

    let fb = hd::render_hd_framebuffer(
        protein,
        camera,
        color_scheme,
        viz_mode,
        width as f64,
        height as f64,
        &mesh,
        /* show_ligands */ true,
        /* interactions */ &[],
    );
    Ok(fb)
}

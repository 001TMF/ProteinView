use ratatui::widgets::Paragraph;

use crate::app::VizMode;
use crate::model::protein::Protein;
use crate::render::camera::Camera;
use crate::render::color::ColorScheme;
use crate::render::framebuffer::{default_light_dir, framebuffer_to_widget, Framebuffer, Triangle};
use crate::render::ribbon::generate_ribbon_mesh;

fn color_to_rgb(color: ratatui::style::Color) -> [u8; 3] {
    match color {
        ratatui::style::Color::Rgb(r, g, b) => [r, g, b],
        _ => [180, 180, 180],
    }
}

/// Apply the camera's rotation to a direction vector (no zoom/pan).
fn rotate_normal(camera: &Camera, nx: f64, ny: f64, nz: f64) -> [f64; 3] {
    let (sin_x, cos_x) = camera.rot_x.sin_cos();
    let y1 = ny * cos_x - nz * sin_x;
    let z1 = ny * sin_x + nz * cos_x;

    let (sin_y, cos_y) = camera.rot_y.sin_cos();
    let x2 = nx * cos_y + z1 * sin_y;
    let z2 = -nx * sin_y + z1 * cos_y;

    let (sin_z, cos_z) = camera.rot_z.sin_cos();
    let x3 = x2 * cos_z - y1 * sin_z;
    let y3 = x2 * sin_z + y1 * cos_z;

    [x3, y3, z2]
}

/// Convert projected coords (centered at origin) to pixel coords (top-left origin).
#[inline]
fn to_pixel(proj_x: f64, proj_y: f64, proj_z: f64, half_w: f64, half_h: f64) -> [f64; 3] {
    [proj_x + half_w, half_h - proj_y, proj_z]
}

/// Render protein in HD mode using software rasterization.
///
/// Returns a `Paragraph` widget built from half-block characters with per-pixel
/// color and Lambert shading.
///
/// - `width`: pixel width (= terminal columns)
/// - `height`: pixel height (= terminal rows * 2)
pub fn render_hd(
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
    viz_mode: VizMode,
    width: f64,
    height: f64,
) -> Paragraph<'static> {
    let px_w = width as usize;
    let px_h = height as usize;
    if px_w == 0 || px_h == 0 {
        return Paragraph::new("");
    }

    let mut fb = Framebuffer::new(px_w, px_h);
    let light_dir = default_light_dir();
    let half_w = px_w as f64 / 2.0;
    let half_h = px_h as f64 / 2.0;

    match viz_mode {
        VizMode::Backbone => {
            let mesh = generate_ribbon_mesh(protein, color_scheme);
            for tri in &mesh {
                let v0 = camera.project(tri.verts[0][0], tri.verts[0][1], tri.verts[0][2]);
                let v1 = camera.project(tri.verts[1][0], tri.verts[1][1], tri.verts[1][2]);
                let v2 = camera.project(tri.verts[2][0], tri.verts[2][1], tri.verts[2][2]);

                let rotated_normal =
                    rotate_normal(camera, tri.normal[0], tri.normal[1], tri.normal[2]);

                let screen_tri = Triangle {
                    verts: [
                        to_pixel(v0.x, v0.y, v0.z, half_w, half_h),
                        to_pixel(v1.x, v1.y, v1.z, half_w, half_h),
                        to_pixel(v2.x, v2.y, v2.z, half_w, half_h),
                    ],
                    color: tri.color,
                    normal: rotated_normal,
                };
                fb.rasterize_triangle(&screen_tri, light_dir);
            }
        }
        VizMode::Wireframe => {
            render_wireframe_fb(&mut fb, protein, camera, color_scheme, false, half_w, half_h);
        }
        VizMode::BallAndStick => {
            render_wireframe_fb(&mut fb, protein, camera, color_scheme, true, half_w, half_h);
        }
    }

    framebuffer_to_widget(&fb)
}

fn atoms_bonded_3d(
    a1_x: f64,
    a1_y: f64,
    a1_z: f64,
    a2_x: f64,
    a2_y: f64,
    a2_z: f64,
) -> bool {
    let dx = a2_x - a1_x;
    let dy = a2_y - a1_y;
    let dz = a2_z - a1_z;
    dx * dx + dy * dy + dz * dz <= 1.9 * 1.9
}

/// Render wireframe or ball-and-stick to framebuffer.
fn render_wireframe_fb(
    fb: &mut Framebuffer,
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
    ball_and_stick: bool,
    half_w: f64,
    half_h: f64,
) {
    let atom_radius = if ball_and_stick { 2.0 } else { 0.0 };

    for chain in &protein.chains {
        for residue in &chain.residues {
            let projected: Vec<_> = residue
                .atoms
                .iter()
                .map(|a| {
                    let p = camera.project(a.x, a.y, a.z);
                    let px = to_pixel(p.x, p.y, p.z, half_w, half_h);
                    let color = color_to_rgb(color_scheme.atom_color(a, residue, chain));
                    (a, px, color)
                })
                .collect();

            // Draw atom dots for ball+stick
            if ball_and_stick {
                for (_, px, color) in &projected {
                    fb.draw_circle_z(px[0], px[1], px[2], atom_radius, *color);
                }
            }

            // Intra-residue bonds
            for i in 0..projected.len() {
                for j in (i + 1)..projected.len() {
                    let (a1, p1, c1) = &projected[i];
                    let (a2, p2, _) = &projected[j];
                    if atoms_bonded_3d(a1.x, a1.y, a1.z, a2.x, a2.y, a2.z) {
                        fb.draw_line_3d(*p1, *p2, *c1);
                    }
                }
            }
        }

        // Inter-residue peptide bonds
        for i in 0..chain.residues.len().saturating_sub(1) {
            let res_curr = &chain.residues[i];
            let res_next = &chain.residues[i + 1];
            let c_atom = res_curr.atoms.iter().find(|a| a.name.trim() == "C");
            let n_atom = res_next.atoms.iter().find(|a| a.name.trim() == "N");

            if let (Some(c), Some(n)) = (c_atom, n_atom) {
                let p1 = camera.project(c.x, c.y, c.z);
                let p2 = camera.project(n.x, n.y, n.z);
                let px1 = to_pixel(p1.x, p1.y, p1.z, half_w, half_h);
                let px2 = to_pixel(p2.x, p2.y, p2.z, half_w, half_h);
                let color = color_to_rgb(color_scheme.atom_color(c, res_curr, chain));
                fb.draw_line_3d(px1, px2, color);
            }
        }
    }
}

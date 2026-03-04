use ratatui::symbols::Marker;
use ratatui::widgets::canvas::{Canvas, Context, Line};

use crate::app::VizMode;
use crate::model::protein::Protein;
use crate::render::camera::Camera;
use crate::render::color::ColorScheme;

/// Draw a thick line by rendering parallel offset lines along the perpendicular direction.
fn draw_thick_line(
    ctx: &mut Context<'_>,
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
    color: ratatui::style::Color,
    offsets: &[f64],
) {
    let dx = x2 - x1;
    let dy = y2 - y1;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 0.001 {
        return;
    }

    // Perpendicular direction: (-dy, dx) normalized
    let nx = -dy / len;
    let ny = dx / len;

    for &off in offsets {
        let ox = nx * off;
        let oy = ny * off;
        ctx.draw(&Line {
            x1: x1 + ox,
            y1: y1 + oy,
            x2: x2 + ox,
            y2: y2 + oy,
            color,
        });
    }
}

/// Check whether two atoms are bonded in 3D space.
/// Returns true if they are in the same residue and within 1.9 Angstroms,
/// or if they form a peptide bond (C of residue i to N of residue i+1 in same chain).
fn atoms_bonded_3d(
    a1_x: f64, a1_y: f64, a1_z: f64,
    a2_x: f64, a2_y: f64, a2_z: f64,
) -> bool {
    let dx = a2_x - a1_x;
    let dy = a2_y - a1_y;
    let dz = a2_z - a1_z;
    let dist_sq = dx * dx + dy * dy + dz * dz;
    dist_sq <= 1.9 * 1.9
}

/// Render protein on a ratatui Canvas with the Braille marker.
/// Behavior depends on VizMode:
/// - Backbone/Cartoon: Connect C-alpha atoms with thick lines (5 offsets)
/// - Wireframe: All atoms, thin single bond lines
pub fn render_protein<'a>(
    protein: &'a Protein,
    camera: &'a Camera,
    color_scheme: &'a ColorScheme,
    viz_mode: VizMode,
    width: f64,
    height: f64,
) -> Canvas<'a, impl Fn(&mut Context<'_>) + 'a> {
    Canvas::default()
        .marker(Marker::Braille)
        .x_bounds([-width / 2.0, width / 2.0])
        .y_bounds([-height / 2.0, height / 2.0])
        .paint(move |ctx| {
            match viz_mode {
                VizMode::Backbone | VizMode::Cartoon => {
                    render_backbone(ctx, protein, camera, color_scheme);
                }
                VizMode::Wireframe => {
                    render_wireframe(ctx, protein, camera, color_scheme);
                }
            }
        })
}

/// Backbone rendering: thick lines between consecutive C-alpha atoms.
fn render_backbone(
    ctx: &mut Context<'_>,
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
) {
    let backbone = protein.backbone_atoms();
    if backbone.is_empty() { return; }

    // Perpendicular offsets: centre line + 2 offsets on each side
    let offsets: [f64; 5] = [0.0, 0.3, -0.3, 0.6, -0.6];

    let mut prev: Option<(f64, f64, &str)> = None;
    let mut prev_chain_id = "";

    for (atom, residue, chain) in &backbone {
        let proj = camera.project(atom.x, atom.y, atom.z);

        if chain.id == prev_chain_id {
            if let Some((px, py, _)) = prev {
                let color = color_scheme.residue_color(residue, chain);
                draw_thick_line(ctx, px, py, proj.x, proj.y, color, &offsets);
            }
        }

        prev = Some((proj.x, proj.y, &chain.id));
        prev_chain_id = &chain.id;
    }
}

/// Wireframe rendering: all atoms with thin single bond lines.
fn render_wireframe(
    ctx: &mut Context<'_>,
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
) {
    let offsets: [f64; 1] = [0.0];

    for chain in &protein.chains {
        // Process each residue: intra-residue bonds
        for residue in &chain.residues {
            let atom_count = residue.atoms.len();
            // Project all atoms in this residue
            let projected: Vec<_> = residue.atoms.iter().map(|a| {
                let proj = camera.project(a.x, a.y, a.z);
                (a, proj)
            }).collect();

            // Intra-residue bonds: check all atom pairs within the residue
            for i in 0..atom_count {
                for j in (i + 1)..atom_count {
                    let (a1, p1) = &projected[i];
                    let (a2, p2) = &projected[j];
                    if atoms_bonded_3d(a1.x, a1.y, a1.z, a2.x, a2.y, a2.z) {
                        let color = color_scheme.atom_color(a1, residue, chain);
                        draw_thick_line(ctx, p1.x, p1.y, p2.x, p2.y, color, &offsets);
                    }
                }
            }
        }

        // Inter-residue peptide bonds: C of residue i to N of residue i+1
        for i in 0..(chain.residues.len().saturating_sub(1)) {
            let res_curr = &chain.residues[i];
            let res_next = &chain.residues[i + 1];

            // Find C atom in current residue and N atom in next residue
            let c_atom = res_curr.atoms.iter().find(|a| a.name.trim() == "C");
            let n_atom = res_next.atoms.iter().find(|a| a.name.trim() == "N");

            if let (Some(c), Some(n)) = (c_atom, n_atom) {
                let p1 = camera.project(c.x, c.y, c.z);
                let p2 = camera.project(n.x, n.y, n.z);
                let color = color_scheme.atom_color(c, res_curr, chain);
                draw_thick_line(ctx, p1.x, p1.y, p2.x, p2.y, color, &offsets);
            }
        }
    }
}

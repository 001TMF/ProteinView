use ratatui::symbols::Marker;
use ratatui::widgets::canvas::{Canvas, Context, Line};

use crate::app::VizMode;
use crate::model::protein::{Protein, SecondaryStructure};
use crate::render::camera::Camera;
use crate::render::color::ColorScheme;

/// Draw a thick line by rendering `num_lines` parallel lines spread across
/// `thickness` units in the perpendicular direction.
fn draw_thick_line(
    ctx: &mut Context<'_>,
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
    color: ratatui::style::Color,
    thickness: f64,
    num_lines: usize,
) {
    let dx = x2 - x1;
    let dy = y2 - y1;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 0.001 {
        return;
    }

    // Perpendicular direction
    let nx = -dy / len;
    let ny = dx / len;

    let half = num_lines / 2;
    for i in 0..num_lines {
        let offset = (i as f64 - half as f64) * (thickness / half.max(1) as f64);
        let ox = nx * offset;
        let oy = ny * offset;
        ctx.draw(&Line {
            x1: x1 + ox,
            y1: y1 + oy,
            x2: x2 + ox,
            y2: y2 + oy,
            color,
        });
    }
}

/// Draw an arrowhead by fanning out the line at the endpoint.
/// `base_thickness` is the normal band width; the arrow widens to
/// `arrow_thickness` at the tip.
fn draw_arrow_segment(
    ctx: &mut Context<'_>,
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
    color: ratatui::style::Color,
    base_thickness: f64,
    arrow_thickness: f64,
    num_lines: usize,
) {
    let dx = x2 - x1;
    let dy = y2 - y1;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 0.001 {
        return;
    }

    let nx = -dy / len;
    let ny = dx / len;

    let half = num_lines / 2;
    for i in 0..num_lines {
        let frac = (i as f64 - half as f64) / half.max(1) as f64;
        // Start offset uses base thickness, end offset uses arrow thickness
        let start_off = frac * base_thickness;
        let end_off = frac * arrow_thickness;
        ctx.draw(&Line {
            x1: x1 + nx * start_off,
            y1: y1 + ny * start_off,
            x2: x2 + nx * end_off,
            y2: y2 + ny * end_off,
            color,
        });
    }
}

/// Draw a small cross/dot marker at a position to indicate an atom.
fn draw_atom_dot(ctx: &mut Context<'_>, x: f64, y: f64, color: ratatui::style::Color) {
    let d = 0.4;
    ctx.draw(&Line { x1: x - d, y1: y, x2: x + d, y2: y, color });
    ctx.draw(&Line { x1: x, y1: y - d, x2: x, y2: y + d, color });
}

/// Check whether two atoms are bonded (within 1.9 Angstroms in 3D).
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

/// Render protein as a cartoon / ribbon view on a ratatui Canvas with
/// HalfBlock marker for coloured pixel output.
///
/// Behavior depends on VizMode:
/// - Backbone: Cartoon rendering (thick ribbons for helices, arrows for sheets, thin coils)
/// - Wireframe: All atoms, thin lines, no thick ribbons
/// - BallAndStick: All atoms with dots + thin bond lines
pub fn render_cartoon<'a>(
    protein: &'a Protein,
    camera: &'a Camera,
    color_scheme: &'a ColorScheme,
    viz_mode: VizMode,
    width: f64,
    height: f64,
) -> Canvas<'a, impl Fn(&mut Context<'_>) + 'a> {
    Canvas::default()
        .marker(Marker::HalfBlock)
        .x_bounds([-width / 2.0, width / 2.0])
        .y_bounds([-height / 2.0, height / 2.0])
        .paint(move |ctx| {
            match viz_mode {
                VizMode::Backbone => {
                    render_cartoon_backbone(ctx, protein, camera, color_scheme);
                }
                VizMode::Wireframe => {
                    render_cartoon_wireframe(ctx, protein, camera, color_scheme, false);
                }
                VizMode::BallAndStick => {
                    render_cartoon_wireframe(ctx, protein, camera, color_scheme, true);
                }
            }
        })
}

/// Original cartoon backbone rendering with secondary structure features.
fn render_cartoon_backbone(
    ctx: &mut Context<'_>,
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
) {
    let backbone = protein.backbone_atoms();
    if backbone.is_empty() {
        return;
    }

    // We need a small lookahead to detect the last residue of a sheet
    // run so we can draw the arrowhead. Collect projected points first.
    struct SegInfo {
        px: f64,
        py: f64,
        x: f64,
        y: f64,
        color: ratatui::style::Color,
        ss: SecondaryStructure,
        next_ss: Option<SecondaryStructure>,
    }

    // Build segment list
    let mut segments: Vec<SegInfo> = Vec::new();
    {
        let mut prev_pt: Option<(f64, f64, String)> = None;

        for (i, (atom, residue, chain)) in backbone.iter().enumerate() {
            let proj = camera.project(atom.x, atom.y, atom.z);

            let same_chain = chain.id == prev_pt.as_ref().map(|p| p.2.as_str()).unwrap_or("");

            if same_chain {
                if let Some((px, py, _)) = &prev_pt {
                    let color = color_scheme.residue_color(residue, chain);
                    // Peek at next residue's SS for arrowhead detection
                    let next_ss = backbone.get(i + 1).and_then(|(_, next_res, next_chain)| {
                        if next_chain.id == chain.id {
                            Some(next_res.secondary_structure)
                        } else {
                            None
                        }
                    });
                    segments.push(SegInfo {
                        px: *px,
                        py: *py,
                        x: proj.x,
                        y: proj.y,
                        color,
                        ss: residue.secondary_structure,
                        next_ss,
                    });
                }
            }

            prev_pt = Some((proj.x, proj.y, chain.id.clone()));
        }
    }

    // Draw segments
    for seg in &segments {
        match seg.ss {
            SecondaryStructure::Helix => {
                // Wide ribbon: 11 parallel lines, 1.5 unit thickness
                draw_thick_line(ctx, seg.px, seg.py, seg.x, seg.y, seg.color, 1.5, 11);
            }
            SecondaryStructure::Sheet => {
                // Check if this is the last segment of a sheet run
                // (next residue is not a sheet, or end of chain)
                let is_last = match seg.next_ss {
                    Some(SecondaryStructure::Sheet) => false,
                    _ => true,
                };

                if is_last {
                    // Arrowhead: widen from 1.8 to 3.0 at the tip
                    draw_arrow_segment(
                        ctx, seg.px, seg.py, seg.x, seg.y, seg.color, 1.8, 3.0, 15,
                    );
                } else {
                    // Normal sheet band
                    draw_thick_line(ctx, seg.px, seg.py, seg.x, seg.y, seg.color, 1.8, 13);
                }
            }
            SecondaryStructure::Turn | SecondaryStructure::Coil => {
                // Thin coil/turn: 3 parallel lines, 0.3 unit thickness
                draw_thick_line(ctx, seg.px, seg.py, seg.x, seg.y, seg.color, 0.3, 3);
            }
        }
    }
}

/// Wireframe / Ball+Stick rendering for cartoon (HD) mode.
/// All atoms shown with bonds. If `ball_and_stick` is true, atom dots are drawn
/// and bonds use slightly thicker lines.
fn render_cartoon_wireframe(
    ctx: &mut Context<'_>,
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
    ball_and_stick: bool,
) {
    let bond_thickness = if ball_and_stick { 0.15 } else { 0.0 };
    let bond_lines = if ball_and_stick { 3 } else { 1 };

    for chain in &protein.chains {
        // Process each residue: intra-residue bonds
        for residue in &chain.residues {
            let atom_count = residue.atoms.len();
            // Project all atoms in this residue
            let projected: Vec<_> = residue.atoms.iter().map(|a| {
                let proj = camera.project(a.x, a.y, a.z);
                (a, proj)
            }).collect();

            // Draw atom dots for ball+stick mode
            if ball_and_stick {
                for (atom, proj) in &projected {
                    let color = color_scheme.atom_color(atom, residue, chain);
                    draw_atom_dot(ctx, proj.x, proj.y, color);
                }
            }

            // Intra-residue bonds: check all atom pairs within the residue
            for i in 0..atom_count {
                for j in (i + 1)..atom_count {
                    let (a1, p1) = &projected[i];
                    let (a2, p2) = &projected[j];
                    if atoms_bonded_3d(a1.x, a1.y, a1.z, a2.x, a2.y, a2.z) {
                        let color = color_scheme.atom_color(a1, residue, chain);
                        draw_thick_line(ctx, p1.x, p1.y, p2.x, p2.y, color, bond_thickness, bond_lines);
                    }
                }
            }
        }

        // Inter-residue peptide bonds: C of residue i to N of residue i+1
        for i in 0..(chain.residues.len().saturating_sub(1)) {
            let res_curr = &chain.residues[i];
            let res_next = &chain.residues[i + 1];

            let c_atom = res_curr.atoms.iter().find(|a| a.name.trim() == "C");
            let n_atom = res_next.atoms.iter().find(|a| a.name.trim() == "N");

            if let (Some(c), Some(n)) = (c_atom, n_atom) {
                let p1 = camera.project(c.x, c.y, c.z);
                let p2 = camera.project(n.x, n.y, n.z);
                let color = color_scheme.atom_color(c, res_curr, chain);
                draw_thick_line(ctx, p1.x, p1.y, p2.x, p2.y, color, bond_thickness, bond_lines);
            }
        }
    }
}

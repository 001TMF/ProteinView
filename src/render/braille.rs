use ratatui::symbols::Marker;
use ratatui::widgets::canvas::{Canvas, Context, Line};

use crate::app::VizMode;
use crate::model::interface::{Interaction, InteractionType};
use crate::model::protein::{LigandType, MoleculeType, Protein};
use crate::render::bond::atoms_bonded;
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

/// Render protein on a ratatui Canvas with the Braille marker.
/// Behavior depends on VizMode:
/// - Backbone/Cartoon: Connect C-alpha atoms with thick lines (5 offsets)
/// - Wireframe: All atoms, thin single bond lines
#[allow(clippy::too_many_arguments)]
pub fn render_protein<'a>(
    protein: &'a Protein,
    camera: &'a Camera,
    color_scheme: &'a ColorScheme,
    viz_mode: VizMode,
    width: f64,
    height: f64,
    show_ligands: bool,
    interactions: &'a [Interaction],
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

            if show_ligands {
                render_ligands(ctx, protein, camera, color_scheme);
            }

            if !interactions.is_empty() {
                render_interactions(ctx, interactions, camera);
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
    if backbone.is_empty() {
        return;
    }

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
            let projected: Vec<_> = residue
                .atoms
                .iter()
                .map(|a| {
                    let proj = camera.project(a.x, a.y, a.z);
                    (a, proj)
                })
                .collect();

            // Intra-residue bonds: check all atom pairs within the residue
            for i in 0..atom_count {
                for j in (i + 1)..atom_count {
                    let (a1, p1) = &projected[i];
                    let (a2, p2) = &projected[j];
                    if atoms_bonded(&a1.element, a1.x, a1.y, a1.z, &a2.element, a2.x, a2.y, a2.z) {
                        let color = color_scheme.atom_color(a1, residue, chain);
                        draw_thick_line(ctx, p1.x, p1.y, p2.x, p2.y, color, &offsets);
                    }
                }
            }
        }

        // Inter-residue bonds: peptide (C->N) for proteins,
        // phosphodiester (O3'->P) for nucleic acids
        for i in 0..(chain.residues.len().saturating_sub(1)) {
            let res_curr = &chain.residues[i];
            let res_next = &chain.residues[i + 1];

            let (from_atom, to_atom) = match chain.molecule_type {
                MoleculeType::RNA | MoleculeType::DNA => {
                    let o3 = res_curr.atoms.iter().find(|a| a.name.trim() == "O3'");
                    let p = res_next.atoms.iter().find(|a| a.name.trim() == "P");
                    (o3, p)
                }
                MoleculeType::Protein => {
                    let c = res_curr.atoms.iter().find(|a| a.name.trim() == "C");
                    let n = res_next.atoms.iter().find(|a| a.name.trim() == "N");
                    (c, n)
                }
                MoleculeType::SmallMolecule => (None, None),
            };

            if let (Some(a1), Some(a2)) = (from_atom, to_atom) {
                let p1 = camera.project(a1.x, a1.y, a1.z);
                let p2 = camera.project(a2.x, a2.y, a2.z);
                let color = color_scheme.atom_color(a1, res_curr, chain);
                draw_thick_line(ctx, p1.x, p1.y, p2.x, p2.y, color, &offsets);
            }
        }
    }
}

/// Render small molecules: dots and lines for ligands, single dots for ions.
fn render_ligands(
    ctx: &mut Context<'_>,
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
) {
    let offsets: [f64; 1] = [0.0];

    for ligand in &protein.ligands {
        match ligand.ligand_type {
            LigandType::Ion => {
                // Single dot for ions — draw as a small cross for visibility
                if let Some(atom) = ligand.atoms.first() {
                    let proj = camera.project(atom.x, atom.y, atom.z);
                    let color = color_scheme.ligand_atom_color(atom, ligand);
                    let sz = 0.5;
                    ctx.draw(&Line {
                        x1: proj.x - sz,
                        y1: proj.y,
                        x2: proj.x + sz,
                        y2: proj.y,
                        color,
                    });
                    ctx.draw(&Line {
                        x1: proj.x,
                        y1: proj.y - sz,
                        x2: proj.x,
                        y2: proj.y + sz,
                        color,
                    });
                }
            }
            LigandType::Ligand => {
                // Project all atoms
                let projected: Vec<_> = ligand
                    .atoms
                    .iter()
                    .map(|a| {
                        let proj = camera.project(a.x, a.y, a.z);
                        let color = color_scheme.ligand_atom_color(a, ligand);
                        (a, proj, color)
                    })
                    .collect();

                // Draw intra-ligand bonds
                for i in 0..projected.len() {
                    for j in (i + 1)..projected.len() {
                        let (a1, p1, color) = &projected[i];
                        let (a2, p2, _) = &projected[j];
                        if atoms_bonded(
                            &a1.element,
                            a1.x,
                            a1.y,
                            a1.z,
                            &a2.element,
                            a2.x,
                            a2.y,
                            a2.z,
                        ) {
                            draw_thick_line(ctx, p1.x, p1.y, p2.x, p2.y, *color, &offsets);
                        }
                    }
                }
            }
        }
    }
}

/// Map interaction type to a ratatui color for braille rendering.
fn braille_interaction_color(t: InteractionType) -> ratatui::style::Color {
    match t {
        InteractionType::HydrogenBond => ratatui::style::Color::Rgb(0, 220, 255), // cyan
        InteractionType::SaltBridge => ratatui::style::Color::Rgb(255, 80, 80),   // red
        InteractionType::HydrophobicContact => ratatui::style::Color::Rgb(220, 200, 60), // yellow
        InteractionType::Other => ratatui::style::Color::Rgb(160, 160, 160),      // gray
    }
}

/// Render non-covalent interaction lines in braille mode.
///
/// Uses solid lines (not dashed) because braille resolution is too low for
/// dashes to be visually meaningful.
fn render_interactions(ctx: &mut Context<'_>, interactions: &[Interaction], camera: &Camera) {
    for interaction in interactions {
        let p1 = camera.project(
            interaction.atom_a[0],
            interaction.atom_a[1],
            interaction.atom_a[2],
        );
        let p2 = camera.project(
            interaction.atom_b[0],
            interaction.atom_b[1],
            interaction.atom_b[2],
        );
        let color = braille_interaction_color(interaction.interaction_type);
        ctx.draw(&Line {
            x1: p1.x,
            y1: p1.y,
            x2: p2.x,
            y2: p2.y,
            color,
        });
    }
}

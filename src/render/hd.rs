use crate::app::VizMode;
use crate::model::interface::{Interaction, InteractionType};
use crate::model::protein::{LigandType, MoleculeType, Protein};
use crate::render::bond::atoms_bonded;
use crate::render::camera::Camera;
use crate::render::color::{color_to_rgb, ColorScheme};
use crate::render::framebuffer::{default_light_dir, Framebuffer, Triangle};
use crate::render::ribbon::RibbonTriangle;

/// Render the protein into a raw [`Framebuffer`] at the given pixel dimensions.
///
/// This is the core rasterization entry-point.  Callers decide how to present
/// the result -- either via braille characters or via a graphics-protocol
/// image (Sixel / Kitty) through ratatui-image.
#[allow(clippy::too_many_arguments)]
pub fn render_hd_framebuffer(
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
    viz_mode: VizMode,
    width: f64,
    height: f64,
    mesh: &[RibbonTriangle],
    show_ligands: bool,
    interactions: &[Interaction],
) -> Framebuffer {
    let px_w = width as usize;
    let px_h = height as usize;
    if px_w == 0 || px_h == 0 {
        return Framebuffer::new(1, 1);
    }

    let mut fb = Framebuffer::new(px_w, px_h);
    let light_dir = default_light_dir();
    let half_w = px_w as f64 / 2.0;
    let half_h = px_h as f64 / 2.0;

    // Scale line thickness and circle radii relative to framebuffer size.
    // Values were tuned at ~160px wide (braille resolution) where 1.5px
    // lines and circles look correct.  At FullHD (~640px+) we scale up
    // proportionally.  Floor of 1.0 preserves the original look at low
    // resolutions; ceiling of 3.0 caps growth on 4K terminals.
    let ts = (px_w as f64 / 500.0).clamp(1.0, 3.0);

    match viz_mode {
        VizMode::Cartoon => {
            for tri in mesh {
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
                fb.rasterize_triangle_depth(&screen_tri, light_dir);
            }
        }
        VizMode::Backbone => {
            render_backbone_fb(&mut fb, protein, camera, color_scheme, half_w, half_h, ts);
        }
        VizMode::Wireframe => {
            render_wireframe_fb(&mut fb, protein, camera, color_scheme, half_w, half_h, ts);
        }
    }

    // Render small molecules as ball-and-stick overlay
    if show_ligands {
        render_ligands_fb(&mut fb, protein, camera, color_scheme, half_w, half_h, ts);
    }

    // Post-pass: blend all rasterized pixels toward a cool blue-gray fog color
    // based on their z-buffer depth.  This gives uniform depth cues across all
    // rendering modes (triangles, lines, circles).
    fb.apply_depth_tint([40, 50, 70], 0.35);

    // Render interaction lines AFTER depth tint so their color coding stays vivid.
    if !interactions.is_empty() {
        render_interactions_fb(&mut fb, interactions, camera, half_w, half_h);
    }

    fb
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

/// Render backbone CA trace to framebuffer.
fn render_backbone_fb(
    fb: &mut Framebuffer,
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
    half_w: f64,
    half_h: f64,
    ts: f64,
) {
    for chain in &protein.chains {
        let mut prev: Option<([f64; 3], [u8; 3])> = None;
        for residue in &chain.residues {
            if let Some(ca) = residue.atoms.iter().find(|a| a.is_backbone) {
                let p = camera.project(ca.x, ca.y, ca.z);
                let px = to_pixel(p.x, p.y, p.z, half_w, half_h);
                let color = color_to_rgb(color_scheme.residue_color(residue, chain));
                fb.draw_circle_z(px[0], px[1], px[2], 2.5 * ts, color);
                if let Some((prev_px, prev_color)) = prev {
                    fb.draw_thick_line_3d(prev_px, px, prev_color, 2.0 * ts);
                }
                prev = Some((px, color));
            }
        }
    }
}

/// Render wireframe mode to framebuffer.
///
/// All atoms are always rendered (the integer underflow fix in `draw_circle_z`
/// prevents the freeze that previously required skipping atoms for large
/// proteins).  Small dots are drawn at every atom position so that atoms are
/// visible at bond intersections.
fn render_wireframe_fb(
    fb: &mut Framebuffer,
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
    half_w: f64,
    half_h: f64,
    ts: f64,
) {
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

            // Draw small dots at atom positions so atoms are visible at bond
            // intersections.
            for (_, px, color) in &projected {
                fb.draw_circle_z(px[0], px[1], px[2], 1.5 * ts, *color);
            }

            // Intra-residue bonds (thick lines)
            for i in 0..projected.len() {
                for j in (i + 1)..projected.len() {
                    let (a1, p1, c1) = &projected[i];
                    let (a2, p2, _) = &projected[j];
                    if atoms_bonded(&a1.element, a1.x, a1.y, a1.z, &a2.element, a2.x, a2.y, a2.z) {
                        fb.draw_thick_line_3d(*p1, *p2, *c1, 1.5 * ts);
                    }
                }
            }
        }

        // Inter-residue bonds: peptide (C->N) for proteins,
        // phosphodiester (O3'->P) for nucleic acids
        for i in 0..chain.residues.len().saturating_sub(1) {
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
                let px1 = to_pixel(p1.x, p1.y, p1.z, half_w, half_h);
                let px2 = to_pixel(p2.x, p2.y, p2.z, half_w, half_h);
                let color = color_to_rgb(color_scheme.atom_color(a1, res_curr, chain));
                fb.draw_thick_line_3d(px1, px2, color, 1.5 * ts);
            }
        }
    }
}

/// Render small molecules: ball-and-stick for ligands, spheres for ions.
fn render_ligands_fb(
    fb: &mut Framebuffer,
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
    half_w: f64,
    half_h: f64,
    ts: f64,
) {
    for ligand in &protein.ligands {
        match ligand.ligand_type {
            LigandType::Ion => {
                // Single sphere for ions (larger radius)
                if let Some(atom) = ligand.atoms.first() {
                    let p = camera.project(atom.x, atom.y, atom.z);
                    let px = to_pixel(p.x, p.y, p.z, half_w, half_h);
                    let color = color_to_rgb(color_scheme.ligand_atom_color(atom, ligand));
                    fb.draw_circle_z(px[0], px[1], px[2], 4.5 * ts, color);
                }
            }
            LigandType::Ligand => {
                // Ball-and-stick: atom spheres + bond sticks
                let projected: Vec<_> = ligand
                    .atoms
                    .iter()
                    .map(|a| {
                        let p = camera.project(a.x, a.y, a.z);
                        let px = to_pixel(p.x, p.y, p.z, half_w, half_h);
                        let color = color_to_rgb(color_scheme.ligand_atom_color(a, ligand));
                        (a, px, color)
                    })
                    .collect();

                // Draw atom spheres (radius varies by element)
                for (atom, px, color) in &projected {
                    let radius = match atom.element.trim() {
                        "H" => 1.5,
                        "C" => 2.5,
                        "N" | "O" | "S" => 2.8,
                        "P" => 3.0,
                        "FE" | "Fe" | "ZN" | "Zn" | "MG" | "Mg" => 3.5,
                        _ => 2.5,
                    } * ts;
                    fb.draw_circle_z(px[0], px[1], px[2], radius, *color);
                }

                // Draw bonds between nearby atoms
                for i in 0..projected.len() {
                    for j in (i + 1)..projected.len() {
                        let (a1, p1, c1) = &projected[i];
                        let (a2, p2, _) = &projected[j];
                        if atoms_bonded(&a1.element, a1.x, a1.y, a1.z, &a2.element, a2.x, a2.y, a2.z) {
                            fb.draw_thick_line_3d(*p1, *p2, *c1, 1.5 * ts);
                        }
                    }
                }
            }
        }
    }
}

/// Render non-covalent interaction lines as dashed segments in the framebuffer.
fn render_interactions_fb(
    fb: &mut Framebuffer,
    interactions: &[Interaction],
    camera: &Camera,
    half_w: f64,
    half_h: f64,
) {
    for interaction in interactions {
        let p1 = camera.project(interaction.atom_a[0], interaction.atom_a[1], interaction.atom_a[2]);
        let p2 = camera.project(interaction.atom_b[0], interaction.atom_b[1], interaction.atom_b[2]);
        let px1 = to_pixel(p1.x, p1.y, p1.z, half_w, half_h);
        let px2 = to_pixel(p2.x, p2.y, p2.z, half_w, half_h);
        let color = interaction_color(interaction.interaction_type);
        fb.draw_dashed_line_3d(px1, px2, color, 4.0, 3.0);
    }
}

/// Map interaction type to an RGB color for rendering.
fn interaction_color(t: InteractionType) -> [u8; 3] {
    match t {
        InteractionType::HydrogenBond => [0, 220, 255],       // cyan
        InteractionType::SaltBridge => [255, 80, 80],          // red
        InteractionType::HydrophobicContact => [220, 200, 60], // yellow
        InteractionType::Other => [160, 160, 160],             // gray
    }
}

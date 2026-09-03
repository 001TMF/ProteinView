use crate::app::VizMode;
use crate::model::interface::{Interaction, InteractionType};
use crate::model::protein::{LigandType, MoleculeType, Protein};
use crate::render::bond::atoms_bonded;
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, color_to_rgb};
use crate::render::framebuffer::{Framebuffer, default_light_dir};
use crate::render::ribbon::RibbonTriangle;
use rayon::prelude::*;

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
    render_hd_framebuffer_ssaa(
        protein,
        camera,
        color_scheme,
        viz_mode,
        width,
        height,
        mesh,
        show_ligands,
        interactions,
        1.0,
    )
}

/// Like [`render_hd_framebuffer`], but aware that the caller intends to
/// downsample the result by a factor of `ssaa` before display.
///
/// `width` / `height` are the *supersampled* framebuffer dimensions, i.e. the
/// output resolution already multiplied by `ssaa`.  The factor is needed
/// separately because line thickness and circle radii must be derived from the
/// **output** resolution and then scaled up, so that features downsample to the
/// same apparent size rather than becoming proportionally thinner.
#[allow(clippy::too_many_arguments)]
pub fn render_hd_framebuffer_ssaa(
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
    viz_mode: VizMode,
    width: f64,
    height: f64,
    mesh: &[RibbonTriangle],
    show_ligands: bool,
    interactions: &[Interaction],
    ssaa: f64,
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

    // Scale line thickness and circle radii relative to the *output* size.
    // Values were tuned at ~160px wide (braille resolution) where 1.5px
    // lines and circles look correct.  At FullHD (~640px+) we scale up
    // proportionally.  Floor of 1.0 preserves the original look at low
    // resolutions; ceiling of 3.0 caps growth on 4K terminals.
    //
    // When supersampling, the clamp must be evaluated against the resolution
    // the user actually sees and only then multiplied by `ssaa`.  Deriving it
    // from the supersampled width instead would let the clamp floor absorb the
    // factor and render features too thin once downsampled.
    let ssaa = if ssaa.is_finite() && ssaa >= 1.0 {
        ssaa
    } else {
        1.0
    };
    let output_px_w = px_w as f64 / ssaa;
    let ts = (output_px_w / 500.0).clamp(1.0, 3.0) * ssaa;

    // Pre-compute sin/cos once for the entire frame instead of per-vertex.
    let cache = camera.projection_cache();

    match viz_mode {
        VizMode::Cartoon => {
            let ctx = TiledRenderCtx {
                half_w,
                half_h,
                px_w,
                px_h,
                light_dir,
            };
            render_cartoon_tiled(&mut fb, mesh, &cache, &ctx);
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

/// Convert projected coords (centered at origin) to pixel coords (top-left origin).
#[inline]
fn to_pixel(proj_x: f64, proj_y: f64, proj_z: f64, half_w: f64, half_h: f64) -> [f64; 3] {
    [proj_x + half_w, half_h - proj_y, proj_z]
}

// ---------------------------------------------------------------------------
// Tile-based parallel cartoon rasterization
// ---------------------------------------------------------------------------

/// Tile size in pixels.  64x64 is a good balance between parallelism (many
/// tiles) and per-tile overhead (triangle binning, allocation).
const TILE_SIZE: usize = 64;

/// A projected, shaded triangle ready for rasterization into tiles.
struct ProjectedTriangle {
    /// Screen-space vertices `[x, y, z]`.
    verts: [[f64; 3]; 3],
    /// Pre-computed flat-shaded color (Lambert applied).
    shaded: [u8; 3],
    /// Screen-space bounding box (clamped to framebuffer).
    min_x: usize,
    max_x: usize,
    min_y: usize,
    max_y: usize,
}

/// Context for tile-based cartoon rasterization, reducing parameter count.
struct TiledRenderCtx {
    half_w: f64,
    half_h: f64,
    px_w: usize,
    px_h: usize,
    light_dir: [f64; 3],
}

/// A rasterized tile with its position, dimensions, and pixel data.
struct RenderedTile {
    x: usize,
    w: usize,
    h: usize,
    color: Vec<[u8; 3]>,
    depth: Vec<f32>,
}

/// Render the cartoon mesh using tile-based parallel rasterization.
///
/// 1. Project all triangles serially (Lambert shade).
/// 2. Bin projected triangles into screen-space tiles.
/// 3. Rasterize each tile in parallel via rayon -- each tile owns its own
///    color/depth arrays so no synchronization is needed.
/// 4. Merge tile results back into the main framebuffer.
fn render_cartoon_tiled(
    fb: &mut Framebuffer,
    mesh: &[RibbonTriangle],
    cache: &crate::render::camera::ProjectionCache,
    ctx: &TiledRenderCtx,
) {
    const AMBIENT: f64 = 0.55;
    let half_w = ctx.half_w;
    let half_h = ctx.half_h;
    let px_w = ctx.px_w;
    let px_h = ctx.px_h;
    let light_dir = ctx.light_dir;

    // ------------------------------------------------------------------
    // Step 1: Project and shade all triangles (serial).
    // ------------------------------------------------------------------
    let projected: Vec<ProjectedTriangle> = mesh
        .par_iter()
        .filter_map(|tri| {
            let v0 = cache.project(tri.verts[0][0], tri.verts[0][1], tri.verts[0][2]);
            let v1 = cache.project(tri.verts[1][0], tri.verts[1][1], tri.verts[1][2]);
            let v2 = cache.project(tri.verts[2][0], tri.verts[2][1], tri.verts[2][2]);

            let sv0 = to_pixel(v0.x, v0.y, v0.z, half_w, half_h);
            let sv1 = to_pixel(v1.x, v1.y, v1.z, half_w, half_h);
            let sv2 = to_pixel(v2.x, v2.y, v2.z, half_w, half_h);

            // Screen-space bounding box clamped to framebuffer.
            let fmin_x = sv0[0].min(sv1[0]).min(sv2[0]).floor() as isize;
            let fmax_x = sv0[0].max(sv1[0]).max(sv2[0]).ceil() as isize;
            let fmin_y = sv0[1].min(sv1[1]).min(sv2[1]).floor() as isize;
            let fmax_y = sv0[1].max(sv1[1]).max(sv2[1]).ceil() as isize;

            let min_x = fmin_x.max(0) as usize;
            let max_x = (fmax_x.max(0) as usize).min(px_w.saturating_sub(1));
            let min_y = fmin_y.max(0) as usize;
            let max_y = (fmax_y.max(0) as usize).min(px_h.saturating_sub(1));

            if min_x > max_x || min_y > max_y {
                return None;
            }

            // Two-sided half-Lambert shading (identical to `rasterize_triangle_depth`).
            let rn = cache.rotate_normal(tri.normal[0], tri.normal[1], tri.normal[2]);
            let dot = rn[0] * light_dir[0] + rn[1] * light_dir[1] + rn[2] * light_dir[2];
            let half_lambert = dot.abs() * 0.4 + 0.6;
            let intensity = AMBIENT + (1.0 - AMBIENT) * half_lambert;
            let shaded: [u8; 3] = [
                (tri.color[0] as f64 * intensity).min(255.0) as u8,
                (tri.color[1] as f64 * intensity).min(255.0) as u8,
                (tri.color[2] as f64 * intensity).min(255.0) as u8,
            ];

            Some(ProjectedTriangle {
                verts: [sv0, sv1, sv2],
                shaded,
                min_x,
                max_x,
                min_y,
                max_y,
            })
        })
        .collect();

    if projected.is_empty() {
        return;
    }

    // ------------------------------------------------------------------
    // Step 2: Create tile grid and bin triangles into tiles.
    // ------------------------------------------------------------------
    let cols = px_w.div_ceil(TILE_SIZE);
    let rows = px_h.div_ceil(TILE_SIZE);
    let num_tiles = cols * rows;
    let mut tile_triangles: Vec<Vec<usize>> = vec![Vec::new(); num_tiles];

    for (tri_idx, tri) in projected.iter().enumerate() {
        let tc0 = tri.min_x / TILE_SIZE;
        let tc1 = tri.max_x / TILE_SIZE;
        let tr0 = tri.min_y / TILE_SIZE;
        let tr1 = tri.max_y / TILE_SIZE;
        for tr in tr0..=tr1 {
            for tc in tc0..=tc1 {
                tile_triangles[tr * cols + tc].push(tri_idx);
            }
        }
    }

    // ------------------------------------------------------------------
    // Step 3: Rasterize tiles in parallel.
    // ------------------------------------------------------------------
    // Each tile produces its own local color and depth buffers.
    let tiles: Vec<RenderedTile> = tile_triangles
        .into_par_iter()
        .enumerate()
        .map(|(tile_idx, tri_indices)| {
            let tx = (tile_idx % cols) * TILE_SIZE;
            let ty = (tile_idx / cols) * TILE_SIZE;
            let tw = TILE_SIZE.min(px_w.saturating_sub(tx));
            let th = TILE_SIZE.min(px_h.saturating_sub(ty));

            let mut color = vec![[0u8; 3]; tw * th];
            let mut depth = vec![f32::INFINITY; tw * th];

            for &tri_idx in &tri_indices {
                let tri = &projected[tri_idx];
                rasterize_into_tile(&mut color, &mut depth, tw, th, tx, ty, tri);
            }

            RenderedTile {
                x: tx,
                w: tw,
                h: th,
                color,
                depth,
            }
        })
        .collect();

    // ------------------------------------------------------------------
    // Step 4: Merge tiles back into the main framebuffer.
    // ------------------------------------------------------------------
    // Tiles cover disjoint screen rectangles, so the merge parallelizes cleanly
    // over framebuffer rows: row `y` is covered by exactly the `cols` tiles in
    // tile-row `y / TILE_SIZE`.
    let tiles = &tiles;
    fb.color
        .par_chunks_mut(px_w)
        .zip(fb.depth.par_chunks_mut(px_w))
        .enumerate()
        .for_each(|(y, (color_row, depth_row))| {
            let tr = y / TILE_SIZE;
            let ly = y % TILE_SIZE;
            for tc in 0..cols {
                let tile = &tiles[tr * cols + tc];
                if ly >= tile.h {
                    continue;
                }
                let src_base = ly * tile.w;
                for lx in 0..tile.w {
                    let ti = src_base + lx;
                    let fi = tile.x + lx;
                    if tile.depth[ti] < depth_row[fi] {
                        color_row[fi] = tile.color[ti];
                        depth_row[fi] = tile.depth[ti];
                    }
                }
            }
        });
}

/// Rasterize a single projected triangle into a tile's local buffers.
///
/// The algorithm is identical to `Framebuffer::rasterize_triangle_depth` but
/// operates on tile-local coordinate arrays.  `tx`/`ty` are the pixel
/// coordinates of the tile's top-left corner in the full framebuffer.
#[inline]
fn rasterize_into_tile(
    color: &mut [[u8; 3]],
    depth: &mut [f32],
    tw: usize,
    th: usize,
    tx: usize,
    ty: usize,
    tri: &ProjectedTriangle,
) {
    let [v0, v1, v2] = tri.verts;

    // Clamp the triangle's bounding box to this tile.
    let min_x = tri.min_x.max(tx);
    let max_x = tri.max_x.min((tx + tw).saturating_sub(1));
    let min_y = tri.min_y.max(ty);
    let max_y = tri.max_y.min((ty + th).saturating_sub(1));

    if min_x > max_x || min_y > max_y {
        return;
    }

    // Barycentric denominator (same math as `rasterize_triangle_depth`).
    let denom = (v1[1] - v2[1]) * (v0[0] - v2[0]) + (v2[0] - v1[0]) * (v0[1] - v2[1]);
    if denom.abs() < 1e-12 {
        return; // degenerate triangle
    }
    let inv_denom = 1.0 / denom;

    let u_x_step = (v1[1] - v2[1]) * inv_denom;
    let v_x_step = (v2[1] - v0[1]) * inv_denom;
    let u_y_coeff = (v2[0] - v1[0]) * inv_denom;
    let v_y_coeff = (v0[0] - v2[0]) * inv_denom;

    for py in min_y..=max_y {
        let pf_y = py as f64 + 0.5;
        let dy = pf_y - v2[1];
        let u_y = u_y_coeff * dy;
        let v_y = v_y_coeff * dy;

        for px in min_x..=max_x {
            let pf_x = px as f64 + 0.5;
            let dx = pf_x - v2[0];

            let u = u_x_step * dx + u_y;
            let v = v_x_step * dx + v_y;
            let w = 1.0 - u - v;

            if u >= -1e-6 && v >= -1e-6 && w >= -1e-6 {
                let z = (u * v0[2] + v * v1[2] + w * v2[2]) as f32;
                let lx = px - tx;
                let ly = py - ty;
                let ti = ly * tw + lx;
                if z < depth[ti] {
                    depth[ti] = z;
                    color[ti] = tri.shaded;
                }
            }
        }
    }
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
        let px1 = to_pixel(p1.x, p1.y, p1.z, half_w, half_h);
        let px2 = to_pixel(p2.x, p2.y, p2.z, half_w, half_h);
        let color = interaction_color(interaction.interaction_type);
        fb.draw_dashed_line_3d(px1, px2, color, 4.0, 3.0);
    }
}

/// Map interaction type to an RGB color for rendering.
fn interaction_color(t: InteractionType) -> [u8; 3] {
    match t {
        InteractionType::HydrogenBond => [0, 220, 255], // cyan
        InteractionType::SaltBridge => [255, 80, 80],   // red
        InteractionType::HydrophobicContact => [220, 200, 60], // yellow
        InteractionType::Other => [160, 160, 160],      // gray
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::protein::{Atom, Chain, MoleculeType, Residue, SecondaryStructure};
    use crate::render::color::ColorSchemeType;

    /// A short CA trace running diagonally across the view.
    fn trace_protein() -> Protein {
        let atom = |x: f64, y: f64| Atom {
            name: "CA".to_string(),
            element: "C".to_string(),
            x,
            y,
            z: 0.0,
            b_factor: 20.0,
            is_backbone: true,
            is_hetero: false,
        };
        let residues = (0..12)
            .map(|i| Residue {
                name: "ALA".to_string(),
                seq_num: i + 1,
                insertion_code: None,
                atoms: vec![atom(i as f64 * 6.0 - 33.0, i as f64 * 3.0 - 16.0)],
                secondary_structure: SecondaryStructure::Coil,
            })
            .collect();
        Protein {
            name: "trace".to_string(),
            chains: vec![Chain {
                id: "A".to_string(),
                residues,
                molecule_type: MoleculeType::Protein,
            }],
            ligands: Vec::new(),
        }
    }

    /// Fraction of framebuffer pixels carrying geometry.
    fn ink_fraction(fb: &Framebuffer) -> f64 {
        let lit = fb.color.iter().filter(|c| **c != [0, 0, 0]).count();
        lit as f64 / fb.color.len() as f64
    }

    /// Bounding box of lit pixels, expressed as fractions of the framebuffer
    /// dimensions.  This is the protein's *apparent* size: what the viewer sees
    /// once the buffer is downsampled onto the terminal cell grid.
    fn normalized_bbox(fb: &Framebuffer) -> (f64, f64) {
        let (mut min_x, mut max_x) = (usize::MAX, 0usize);
        let (mut min_y, mut max_y) = (usize::MAX, 0usize);
        for y in 0..fb.height {
            for x in 0..fb.width {
                if fb.color[y * fb.width + x] != [0, 0, 0] {
                    min_x = min_x.min(x);
                    max_x = max_x.max(x);
                    min_y = min_y.min(y);
                    max_y = max_y.max(y);
                }
            }
        }
        assert!(min_x != usize::MAX, "fixture should draw something");
        (
            (max_x - min_x) as f64 / fb.width as f64,
            (max_y - min_y) as f64 / fb.height as f64,
        )
    }

    /// Render the fixture the way the viewport does: the framebuffer is scaled
    /// by `ssaa`, and the camera is scaled to match because zoom and pan are
    /// both in framebuffer pixel units.
    fn render_trace(out_w: f64, out_h: f64, ssaa: f64) -> Framebuffer {
        let protein = trace_protein();
        let mut camera = Camera::default();
        camera.zoom = 4.0 * ssaa;
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 12);
        render_hd_framebuffer_ssaa(
            &protein,
            &camera,
            &scheme,
            VizMode::Backbone,
            out_w * ssaa,
            out_h * ssaa,
            &[],
            false,
            &[],
            ssaa,
        )
    }

    #[test]
    fn supersampling_preserves_apparent_size() {
        // HDplus must render the protein at the same apparent size as HD.
        // Scaling the framebuffer by `ssaa` without scaling the camera would
        // leave the protein covering the same pixel count in a buffer twice as
        // wide, i.e. half the size once downsampled.
        let (bw, bh) = normalized_bbox(&render_trace(400.0, 184.0, 1.0));
        let (pw, ph) = normalized_bbox(&render_trace(400.0, 184.0, 2.0));

        assert!(
            (pw - bw).abs() < 0.02 && (ph - bh).abs() < 0.02,
            "supersampled extent ({pw:.3}, {ph:.3}) should match base ({bw:.3}, {bh:.3})"
        );
    }

    #[test]
    fn supersampling_preserves_apparent_line_thickness() {
        // Line thickness must be derived from the *output* resolution and then
        // scaled by the supersampling factor, so strokes occupy the same
        // fraction of the frame.  Deriving it from the supersampled width
        // instead lets the clamp floor in `ts` absorb the factor and renders
        // strokes too thin once downsampled.
        let base = ink_fraction(&render_trace(400.0, 184.0, 1.0));
        let supersampled = ink_fraction(&render_trace(400.0, 184.0, 2.0));

        let ratio = supersampled / base;
        assert!(
            (0.9..=1.1).contains(&ratio),
            "supersampled ink fraction {supersampled:.4} should match base {base:.4} within 10% (ratio {ratio:.3})"
        );
    }

    #[test]
    fn ssaa_factor_is_ignored_when_invalid() {
        // Guards against a NaN or sub-1 factor silently collapsing thickness.
        let protein = trace_protein();
        let mut camera = Camera::default();
        camera.zoom = 4.0;
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 12);
        let render = |ssaa: f64| {
            render_hd_framebuffer_ssaa(
                &protein,
                &camera,
                &scheme,
                VizMode::Backbone,
                400.0,
                184.0,
                &[],
                false,
                &[],
                ssaa,
            )
        };
        let expected = ink_fraction(&render(1.0));
        for bad in [f64::NAN, 0.0, -3.0] {
            assert!(
                (ink_fraction(&render(bad)) - expected).abs() < 1e-9,
                "ssaa {bad} should fall back to 1.0"
            );
        }
    }
}

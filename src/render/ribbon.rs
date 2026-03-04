//! Ribbon/cartoon geometry generation for HD protein rendering.
//!
//! Generates triangle meshes from protein backbone data by:
//! 1. Extracting C-alpha positions per chain
//! 2. Fitting a Catmull-Rom spline through backbone positions
//! 3. Computing Frenet-Serret local coordinate frames along the spline
//! 4. Extruding secondary-structure-dependent cross-sections
//! 5. Connecting consecutive cross-sections into triangle strips
//!
//! The output mesh is in world space; the caller projects through the camera.

use crate::model::protein::{Protein, SecondaryStructure};
use crate::render::color::ColorScheme;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Number of spline subdivisions between each pair of C-alpha atoms.
const SPLINE_SUBDIVISIONS: usize = 14;

/// Number of vertices around the coil/turn tube cross-section.
const COIL_SEGMENTS: usize = 12;

/// Cross-section dimensions (in Angstroms).
const HELIX_HALF_WIDTH: f64 = 1.30;
const HELIX_HALF_HEIGHT: f64 = 0.40;

const SHEET_HALF_WIDTH: f64 = 1.50;
const SHEET_HALF_HEIGHT: f64 = 0.20;

const SHEET_ARROW_HALF_WIDTH: f64 = 2.20;

const COIL_RADIUS: f64 = 0.40;

// ---------------------------------------------------------------------------
// Output type
// ---------------------------------------------------------------------------

/// A single triangle in the ribbon mesh, ready for rasterization.
#[derive(Debug, Clone)]
pub struct RibbonTriangle {
    /// Three vertices in 3D world space, each `[x, y, z]`.
    pub verts: [[f64; 3]; 3],
    /// Base RGB color of this face.
    pub color: [u8; 3],
    /// Outward-facing unit normal of the triangle.
    pub normal: [f64; 3],
}

// ---------------------------------------------------------------------------
// 3-component vector helpers (no external crate)
// ---------------------------------------------------------------------------

type V3 = [f64; 3];

#[inline]
fn v3_add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn v3_sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn v3_dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn v3_len(a: V3) -> f64 {
    v3_dot(a, a).sqrt()
}

#[inline]
fn v3_normalize(a: V3) -> V3 {
    let l = v3_len(a);
    if l < 1e-12 {
        [0.0, 1.0, 0.0] // fallback up vector
    } else {
        v3_scale(a, 1.0 / l)
    }
}

// ---------------------------------------------------------------------------
// Color conversion
// ---------------------------------------------------------------------------

fn color_to_rgb(color: ratatui::style::Color) -> [u8; 3] {
    match color {
        ratatui::style::Color::Rgb(r, g, b) => [r, g, b],
        _ => [180, 180, 180],
    }
}

// ---------------------------------------------------------------------------
// Catmull-Rom spline
// ---------------------------------------------------------------------------

/// Evaluate the Catmull-Rom spline between `p1` and `p2` at parameter `t` in
/// [0, 1], using `p0` and `p3` as the surrounding control points.
fn catmull_rom(p0: V3, p1: V3, p2: V3, p3: V3, t: f64) -> V3 {
    let t2 = t * t;
    let t3 = t2 * t;

    // q(t) = 0.5 * ((2*P1) + (-P0+P2)*t + (2*P0-5*P1+4*P2-P3)*t^2 + (-P0+3*P1-3*P2+P3)*t^3)
    let mut out = [0.0; 3];
    for i in 0..3 {
        out[i] = 0.5
            * ((2.0 * p1[i])
                + (-p0[i] + p2[i]) * t
                + (2.0 * p0[i] - 5.0 * p1[i] + 4.0 * p2[i] - p3[i]) * t2
                + (-p0[i] + 3.0 * p1[i] - 3.0 * p2[i] + p3[i]) * t3);
    }
    out
}

// ---------------------------------------------------------------------------
// Spline point with metadata
// ---------------------------------------------------------------------------

/// A single point along the backbone spline, carrying the local coordinate
/// frame and secondary-structure annotation.
struct SplinePoint {
    pos: V3,
    tangent: V3,
    normal: V3,
    binormal: V3,
    ss: SecondaryStructure,
    color: [u8; 3],
    /// True when this point lies within the arrowhead region at the end of a
    /// sheet run.  The arrowhead linearly widens over the last two original
    /// residue spans (2 * SPLINE_SUBDIVISIONS points).
    arrow_t: Option<f64>, // 0.0 = start of arrow, 1.0 = tip
}

// ---------------------------------------------------------------------------
// Cross-section generation
// ---------------------------------------------------------------------------

/// Build the cross-section ring for a given spline point.  Returns the
/// world-space positions of the cross-section vertices.
///
/// For ribbons (helix/sheet) we return 2 points (left, right) so that the
/// surface is a flat band.  For coils we return `COIL_SEGMENTS` points in a
/// circle.
fn cross_section(sp: &SplinePoint) -> Vec<V3> {
    match sp.ss {
        SecondaryStructure::Helix => ribbon_cross_section(sp, HELIX_HALF_WIDTH, HELIX_HALF_HEIGHT),
        SecondaryStructure::Sheet => {
            let hw = if let Some(t) = sp.arrow_t {
                // Linearly widen from normal sheet width to arrow width.
                let base = SHEET_HALF_WIDTH;
                let tip = SHEET_ARROW_HALF_WIDTH;
                base + (tip - base) * t
            } else {
                SHEET_HALF_WIDTH
            };
            ribbon_cross_section(sp, hw, SHEET_HALF_HEIGHT)
        }
        SecondaryStructure::Turn | SecondaryStructure::Coil => coil_cross_section(sp),
    }
}

/// Flat ribbon cross-section with 4 vertices (top-left, top-right,
/// bottom-right, bottom-left) forming a thin rectangular profile.
fn ribbon_cross_section(sp: &SplinePoint, half_w: f64, half_h: f64) -> Vec<V3> {
    let n = sp.normal;
    let b = sp.binormal;

    // Four corners of the rectangular cross-section.
    let tl = v3_add(sp.pos, v3_add(v3_scale(b, -half_w), v3_scale(n, half_h)));
    let tr = v3_add(sp.pos, v3_add(v3_scale(b, half_w), v3_scale(n, half_h)));
    let br = v3_add(sp.pos, v3_add(v3_scale(b, half_w), v3_scale(n, -half_h)));
    let bl = v3_add(sp.pos, v3_add(v3_scale(b, -half_w), v3_scale(n, -half_h)));

    vec![tl, tr, br, bl]
}

/// Circular tube cross-section for coil/turn regions.
fn coil_cross_section(sp: &SplinePoint) -> Vec<V3> {
    let n = sp.normal;
    let b = sp.binormal;
    let mut pts = Vec::with_capacity(COIL_SEGMENTS);
    for i in 0..COIL_SEGMENTS {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (COIL_SEGMENTS as f64);
        let (sin_a, cos_a) = angle.sin_cos();
        let offset = v3_add(v3_scale(n, cos_a * COIL_RADIUS), v3_scale(b, sin_a * COIL_RADIUS));
        pts.push(v3_add(sp.pos, offset));
    }
    pts
}

// ---------------------------------------------------------------------------
// Triangle normal
// ---------------------------------------------------------------------------

fn triangle_normal(v0: V3, v1: V3, v2: V3) -> V3 {
    let e1 = v3_sub(v1, v0);
    let e2 = v3_sub(v2, v0);
    v3_normalize(v3_cross(e1, e2))
}

// ---------------------------------------------------------------------------
// Emit triangle strip between two cross-sections
// ---------------------------------------------------------------------------

/// Connect two consecutive cross-section rings with a triangle strip.
/// Both rings must have the same number of vertices.
fn emit_strip(
    ring_a: &[V3],
    ring_b: &[V3],
    color: [u8; 3],
    out: &mut Vec<RibbonTriangle>,
) {
    let n = ring_a.len();
    debug_assert_eq!(n, ring_b.len());
    if n == 0 {
        return;
    }

    for i in 0..n {
        let j = (i + 1) % n;

        let a0 = ring_a[i];
        let a1 = ring_a[j];
        let b0 = ring_b[i];
        let b1 = ring_b[j];

        // Quad (a0, a1, b1, b0) -> two triangles.
        // Triangle 1: a0, a1, b0
        let n1 = triangle_normal(a0, a1, b0);
        out.push(RibbonTriangle {
            verts: [a0, a1, b0],
            color,
            normal: n1,
        });

        // Triangle 2: a1, b1, b0
        let n2 = triangle_normal(a1, b1, b0);
        out.push(RibbonTriangle {
            verts: [a1, b1, b0],
            color,
            normal: n2,
        });
    }
}

/// Emit a cap (disc) to close the end of a tube or ribbon.
/// `ring` is the cross-section ring, `center` is the spline point center,
/// and `facing_forward` controls the winding direction.
fn emit_cap(
    ring: &[V3],
    center: V3,
    color: [u8; 3],
    facing_forward: bool,
    out: &mut Vec<RibbonTriangle>,
) {
    let n = ring.len();
    if n < 3 {
        return;
    }
    for i in 0..n {
        let j = (i + 1) % n;
        let (v0, v1) = if facing_forward {
            (ring[i], ring[j])
        } else {
            (ring[j], ring[i])
        };
        let norm = triangle_normal(center, v0, v1);
        out.push(RibbonTriangle {
            verts: [center, v0, v1],
            color,
            normal: norm,
        });
    }
}

// ---------------------------------------------------------------------------
// Transition cross-sections between different secondary structure types
// ---------------------------------------------------------------------------

/// When the secondary structure type changes between two consecutive spline
/// points the cross-section vertex counts may differ.  We handle this by
/// building an intermediate ring that matches the *other* side's count so
/// that `emit_strip` always gets equal-length rings.
///
/// This simply re-samples the given ring to have `target_count` vertices by
/// linear interpolation around the perimeter.  For the common 4->6 and 6->4
/// transitions this produces a reasonable visual blend.
fn resample_ring(ring: &[V3], target_count: usize) -> Vec<V3> {
    let n = ring.len();
    if n == target_count {
        return ring.to_vec();
    }
    if n == 0 || target_count == 0 {
        return Vec::new();
    }

    let mut out = Vec::with_capacity(target_count);
    for i in 0..target_count {
        let frac = (i as f64) / (target_count as f64) * (n as f64);
        let idx = frac as usize;
        let t = frac - idx as f64;
        let a = ring[idx % n];
        let b = ring[(idx + 1) % n];
        out.push(v3_add(v3_scale(a, 1.0 - t), v3_scale(b, t)));
    }
    out
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Generate the complete ribbon/cartoon triangle mesh for a protein.
///
/// The returned triangles are in world space.  The caller should project each
/// vertex through the camera and then rasterize.
pub fn generate_ribbon_mesh(
    protein: &Protein,
    color_scheme: &ColorScheme,
) -> Vec<RibbonTriangle> {
    let mut triangles: Vec<RibbonTriangle> = Vec::new();

    for chain in &protein.chains {
        generate_chain_ribbon(chain, color_scheme, &mut triangles);
    }

    triangles
}

// ---------------------------------------------------------------------------
// Per-chain generation
// ---------------------------------------------------------------------------

/// C-alpha record extracted from a residue.
struct CaRecord {
    pos: V3,
    ss: SecondaryStructure,
    color: [u8; 3],
}

fn generate_chain_ribbon(
    chain: &crate::model::protein::Chain,
    color_scheme: &ColorScheme,
    out: &mut Vec<RibbonTriangle>,
) {
    // 1. Collect C-alpha positions, SS types, and colors.
    let cas: Vec<CaRecord> = chain
        .residues
        .iter()
        .filter_map(|res| {
            let ca = res.atoms.iter().find(|a| a.is_ca)?;
            let color = color_to_rgb(color_scheme.residue_color(res, chain));
            Some(CaRecord {
                pos: [ca.x, ca.y, ca.z],
                ss: res.secondary_structure,
                color,
            })
        })
        .collect();

    let n = cas.len();
    if n < 2 {
        return;
    }

    // 2. Identify arrowhead regions.  For each residue that is the last in a
    //    contiguous sheet run we want to widen the last two residue spans into
    //    an arrow.  We mark the *residue indices* where an arrow starts.
    //    (An arrow occupies residue indices [arrow_start..=last_sheet].)
    let mut arrow_start: Vec<usize> = Vec::new(); // residue index where arrow begins
    {
        let mut i = 0;
        while i < n {
            if cas[i].ss == SecondaryStructure::Sheet {
                // Find the end of this sheet run.
                let run_start = i;
                while i < n && cas[i].ss == SecondaryStructure::Sheet {
                    i += 1;
                }
                let run_end = i; // exclusive
                let run_len = run_end - run_start;
                // Arrow occupies the last 2 residues of the run (or fewer if
                // the run itself is shorter).
                let arrow_residues = run_len.min(2);
                let start = run_end - arrow_residues;
                arrow_start.push(start);
            } else {
                i += 1;
            }
        }
    }

    // 3. Generate spline points with Catmull-Rom interpolation.
    let mut spline_points: Vec<SplinePoint> = Vec::new();

    for seg in 0..n - 1 {
        // Indices for the four control points, clamping at endpoints.
        let i0 = if seg == 0 { 0 } else { seg - 1 };
        let i1 = seg;
        let i2 = seg + 1;
        let i3 = if seg + 2 >= n { n - 1 } else { seg + 2 };

        let p0 = cas[i0].pos;
        let p1 = cas[i1].pos;
        let p2 = cas[i2].pos;
        let p3 = cas[i3].pos;

        let subdivs = if seg == n - 2 {
            SPLINE_SUBDIVISIONS + 1 // include the last point
        } else {
            SPLINE_SUBDIVISIONS
        };

        for sub in 0..subdivs {
            let t = sub as f64 / SPLINE_SUBDIVISIONS as f64;
            let pos = catmull_rom(p0, p1, p2, p3, t);

            // Interpolated secondary structure: use the nearer residue.
            let ss = if t < 0.5 { cas[i1].ss } else { cas[i2].ss };
            // Interpolated color: use the nearer residue.
            let color = if t < 0.5 { cas[i1].color } else { cas[i2].color };

            // Determine if this point is in an arrowhead region.
            // Arrow region covers the last 2 residue spans of a sheet run.
            let arrow_t = arrow_start.iter().find_map(|&astart| {
                // Arrow spans from residue `astart` to the end of its sheet run.
                // Find end of sheet run from astart.
                let mut aend = astart;
                while aend < n && cas[aend].ss == SecondaryStructure::Sheet {
                    aend += 1;
                }
                // The arrow goes from the first spline sample of residue `astart`
                // to the last sample of residue `aend - 1`.
                let arrow_span = aend - astart; // number of residue segments
                if arrow_span == 0 {
                    return None;
                }

                // Current global spline position as a floating-point residue index.
                let global_pos = seg as f64 + t;
                let arrow_begin = astart as f64;
                let arrow_end = aend as f64; // exclusive in residue space but we approach it
                // We actually want the arrow to span [astart, aend-1] in segment indices,
                // so the last spline sample is at aend.
                if global_pos >= arrow_begin && global_pos <= arrow_end as f64 {
                    let frac = (global_pos - arrow_begin) / (arrow_end as f64 - arrow_begin);
                    Some(frac.clamp(0.0, 1.0))
                } else {
                    None
                }
            });

            spline_points.push(SplinePoint {
                pos,
                tangent: [0.0, 0.0, 0.0], // computed below
                normal: [0.0, 1.0, 0.0],
                binormal: [0.0, 0.0, 1.0],
                ss,
                color,
                arrow_t,
            });
        }
    }

    if spline_points.len() < 2 {
        return;
    }

    // 4. Compute tangents via finite differences.
    let sp_len = spline_points.len();
    for i in 0..sp_len {
        let prev = if i == 0 { 0 } else { i - 1 };
        let next = if i == sp_len - 1 { sp_len - 1 } else { i + 1 };
        let t = v3_normalize(v3_sub(spline_points[next].pos, spline_points[prev].pos));
        spline_points[i].tangent = t;
    }

    // 5. Compute Frenet-Serret frames with a propagated reference normal to
    //    avoid flipping.
    //
    //    We use the "parallel transport" variant: choose an initial normal
    //    perpendicular to the first tangent, then propagate it along the curve
    //    by projecting out the tangent component at each step.
    {
        // Choose initial normal perpendicular to first tangent.
        let t0 = spline_points[0].tangent;
        let arbitrary = if t0[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let mut prev_normal = v3_normalize(v3_cross(t0, arbitrary));

        for sp in spline_points.iter_mut() {
            let t = sp.tangent;
            // Project previous normal onto plane perpendicular to current tangent.
            let proj = v3_scale(t, v3_dot(prev_normal, t));
            let mut n = v3_sub(prev_normal, proj);
            let nl = v3_len(n);
            if nl < 1e-12 {
                // Degenerate: pick a new arbitrary normal.
                let arb = if t[0].abs() < 0.9 {
                    [1.0, 0.0, 0.0]
                } else {
                    [0.0, 1.0, 0.0]
                };
                n = v3_normalize(v3_cross(t, arb));
            } else {
                n = v3_scale(n, 1.0 / nl);
            }
            let b = v3_normalize(v3_cross(t, n));

            sp.normal = n;
            sp.binormal = b;
            prev_normal = n;
        }
    }

    // 6. Build cross-sections and emit triangle strips.
    let mut prev_ring = cross_section(&spline_points[0]);

    for i in 1..spline_points.len() {
        let mut curr_ring = cross_section(&spline_points[i]);
        let color = spline_points[i].color;

        // Handle cross-section vertex count mismatch at SS transitions.
        if prev_ring.len() != curr_ring.len() {
            let target = prev_ring.len().max(curr_ring.len());
            if prev_ring.len() != target {
                prev_ring = resample_ring(&prev_ring, target);
            }
            if curr_ring.len() != target {
                curr_ring = resample_ring(&curr_ring, target);
            }
        }

        emit_strip(&prev_ring, &curr_ring, color, out);
        prev_ring = curr_ring;
    }

    // 7. Cap the ends of the ribbon.
    let first_ring = cross_section(&spline_points[0]);
    emit_cap(
        &first_ring,
        spline_points[0].pos,
        spline_points[0].color,
        false,
        out,
    );

    let last_ring = cross_section(spline_points.last().unwrap());
    emit_cap(
        &last_ring,
        spline_points.last().unwrap().pos,
        spline_points.last().unwrap().color,
        true,
        out,
    );
}

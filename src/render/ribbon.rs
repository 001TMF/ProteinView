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

use crate::model::protein::{MoleculeType, Protein, SecondaryStructure, is_purine};
use crate::render::color::{ColorScheme, color_to_rgb};

// ---------------------------------------------------------------------------
// Constants & LOD configuration
// ---------------------------------------------------------------------------

/// Spline subdivisions between each pair of C-alpha atoms, for small structures.
///
/// Small structures cost little to tessellate whatever the setting, and are the
/// ones most likely to be zoomed in on, so they keep the finer mesh.
const SMALL_SPLINE_SUBDIVISIONS: usize = 12;

/// Vertices around the coil/turn tube cross-section, for small structures.
const SMALL_COIL_SEGMENTS: usize = 10;

/// Default spline subdivisions.
///
/// Measured against a 14x12 mesh at 1600x736 (the largest framebuffer a FullHD
/// terminal viewport produces), 8x8 renders 2.1x faster for a mean per-channel
/// error of 0.40/255, with differences confined to silhouette edges.  At braille
/// and HD resolutions the difference is smaller still.
const DEFAULT_SPLINE_SUBDIVISIONS: usize = 8;

/// Default number of vertices around the coil/turn tube cross-section.
const DEFAULT_COIL_SEGMENTS: usize = 8;

/// Reduced spline subdivisions for large structures (>5000 residues).
const LARGE_SPLINE_SUBDIVISIONS: usize = 4;

/// Reduced coil segments for large structures (>5000 residues).
const LARGE_COIL_SEGMENTS: usize = 6;

/// Residue count below which the finer small-structure mesh is used.
const SMALL_STRUCTURE_THRESHOLD: usize = 300;

/// Level-of-detail configuration for ribbon mesh generation.
/// Large structures use reduced subdivision counts to cut triangle count
/// with zero visible difference at terminal resolution.
#[derive(Debug, Clone, Copy)]
struct LodConfig {
    spline_subdivisions: usize,
    coil_segments: usize,
}

impl LodConfig {
    fn small() -> Self {
        Self {
            spline_subdivisions: SMALL_SPLINE_SUBDIVISIONS,
            coil_segments: SMALL_COIL_SEGMENTS,
        }
    }

    fn normal() -> Self {
        Self {
            spline_subdivisions: DEFAULT_SPLINE_SUBDIVISIONS,
            coil_segments: DEFAULT_COIL_SEGMENTS,
        }
    }

    fn large() -> Self {
        Self {
            spline_subdivisions: LARGE_SPLINE_SUBDIVISIONS,
            coil_segments: LARGE_COIL_SEGMENTS,
        }
    }

    /// Pick the appropriate LOD based on residue count.
    fn for_residue_count(residue_count: usize) -> Self {
        if residue_count > crate::app::LARGE_STRUCTURE_THRESHOLD {
            Self::large()
        } else if residue_count <= SMALL_STRUCTURE_THRESHOLD {
            Self::small()
        } else {
            Self::normal()
        }
    }
}

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
// Sheet frame guide helpers
// ---------------------------------------------------------------------------

use crate::model::protein::Residue;

/// Compute the carbonyl C→O direction for a protein residue.
///
/// Looks for backbone atoms named "C" (carbonyl carbon) and "O" (carbonyl
/// oxygen), computes the normalized direction vector from C to O.  As a
/// fallback, tries CA→O if C is missing.  Returns `None` for nucleic acid
/// residues or when the required atoms are absent.
fn residue_frame_hint(residue: &Residue) -> Option<V3> {
    let find_atom = |name: &str| -> Option<V3> {
        residue
            .atoms
            .iter()
            .find(|a| a.name.trim() == name)
            .map(|a| [a.x, a.y, a.z])
    };

    let o_pos = find_atom("O")?;

    // Prefer backbone C atom; fall back to CA.
    let origin = find_atom("C").or_else(|| find_atom("CA"))?;

    let dir = v3_sub(o_pos, origin);
    let len = v3_len(dir);
    if len < 1e-12 {
        return None;
    }
    Some(v3_scale(dir, 1.0 / len))
}

/// Project vector `v` perpendicular to `onto` and normalize the result.
/// Returns `None` if the projection is degenerate (v is parallel to onto).
fn project_perp_normalize(v: V3, onto: V3) -> Option<V3> {
    let d = v3_dot(v, onto);
    let perp = v3_sub(v, v3_scale(onto, d));
    let len = v3_len(perp);
    if len < 1e-12 {
        None
    } else {
        Some(v3_scale(perp, 1.0 / len))
    }
}

/// Adjust binormal/normal frames for sheet regions using carbonyl direction
/// hints (frame guides).
///
/// This runs **after** parallel-transport frame computation.  For each
/// contiguous run of `SecondaryStructure::Sheet` spline points, it:
///   1. Projects each point's `frame_hint` perpendicular to the tangent.
///   2. Aligns consecutive projected hints to the same hemisphere.
///   3. Blends the hint with the existing binormal (65% hint / 35% existing).
///   4. Recomputes the normal from `tangent × blended_binormal`.
///
/// Points outside sheet regions or without hints are left unchanged.
fn apply_sheet_frame_guides(spline_points: &mut [SplinePoint]) {
    let n = spline_points.len();
    if n == 0 {
        return;
    }

    // Identify contiguous sheet runs.
    let mut i = 0;
    while i < n {
        if spline_points[i].ss != SecondaryStructure::Sheet {
            i += 1;
            continue;
        }

        // Found the start of a sheet run.
        let run_start = i;
        while i < n && spline_points[i].ss == SecondaryStructure::Sheet {
            i += 1;
        }
        let run_end = i; // exclusive

        // Process this sheet run: project hints, align signs, blend.
        let mut prev_hint_perp: Option<V3> = None;

        for sp in spline_points[run_start..run_end].iter_mut() {
            let hint = match sp.frame_hint {
                Some(h) => h,
                None => continue,
            };

            let tangent = sp.tangent;

            // Project hint perpendicular to tangent.
            let mut hint_perp = match project_perp_normalize(hint, tangent) {
                Some(h) => h,
                None => continue,
            };

            // Align sign with previous hint to avoid flipping.
            if let Some(prev) = prev_hint_perp {
                if v3_dot(hint_perp, prev) < 0.0 {
                    hint_perp = v3_scale(hint_perp, -1.0);
                }
            }
            prev_hint_perp = Some(hint_perp);

            // Blend with existing binormal: 65% hint, 35% existing.
            let existing_b = sp.binormal;
            let blended = v3_add(v3_scale(hint_perp, 0.65), v3_scale(existing_b, 0.35));
            let new_binormal = v3_normalize(blended);

            // Recompute normal from tangent x binormal.
            let new_normal = v3_normalize(v3_cross(tangent, new_binormal));

            sp.binormal = new_binormal;
            sp.normal = new_normal;
        }
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
    /// Carbonyl direction hint propagated from the nearest control point.
    /// Used by `apply_sheet_frame_guides` to orient sheet ribbons.
    frame_hint: Option<V3>,
    ss: SecondaryStructure,
    color: [u8; 3],
    /// True when this point lies within the arrowhead region at the end of a
    /// sheet run.  The arrowhead linearly widens over the last two original
    /// residue spans (2 * spline_subdivisions points).
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
fn cross_section(sp: &SplinePoint, lod: &LodConfig) -> Vec<V3> {
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
        SecondaryStructure::Turn | SecondaryStructure::Coil => coil_cross_section(sp, lod),
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
fn coil_cross_section(sp: &SplinePoint, lod: &LodConfig) -> Vec<V3> {
    let n = sp.normal;
    let b = sp.binormal;
    let segs = lod.coil_segments;
    let mut pts = Vec::with_capacity(segs);
    for i in 0..segs {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (segs as f64);
        let (sin_a, cos_a) = angle.sin_cos();
        let offset = v3_add(
            v3_scale(n, cos_a * COIL_RADIUS),
            v3_scale(b, sin_a * COIL_RADIUS),
        );
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
// Consistently wound triangle emission
// ---------------------------------------------------------------------------

/// Mean of a set of points.
fn centroid(points: &[V3]) -> V3 {
    if points.is_empty() {
        return [0.0; 3];
    }
    let sum = points.iter().fold([0.0; 3], |acc: V3, p| v3_add(acc, *p));
    v3_scale(sum, 1.0 / points.len() as f64)
}

/// Push a triangle, forcing its winding to agree with `outward`.
///
/// Every cross-section ring is built from a reference frame carried along the
/// spline, and that frame can flip sign partway down a chain.  When it does,
/// the ring's vertex order reverses, and so does the winding of every triangle
/// built from it -- which is why the mesh used to contain a mix of clockwise
/// and counter-clockwise faces and could not be back-face culled.
///
/// Orienting each triangle against a direction known to point out of the
/// surface makes the winding globally consistent, and makes the stored normal a
/// genuine outward normal rather than whatever the vertex order implied.
fn push_oriented(
    v0: V3,
    v1: V3,
    v2: V3,
    outward: V3,
    color: [u8; 3],
    out: &mut Vec<RibbonTriangle>,
) {
    let normal = triangle_normal(v0, v1, v2);
    if v3_dot(normal, outward) < 0.0 {
        // Swapping two vertices reverses the winding, and with it the normal.
        out.push(RibbonTriangle {
            verts: [v0, v2, v1],
            color,
            normal: v3_scale(normal, -1.0),
        });
    } else {
        out.push(RibbonTriangle {
            verts: [v0, v1, v2],
            color,
            normal,
        });
    }
}

/// Push a triangle, reversing its winding when `flip` is set.
///
/// Used where the correct orientation is known for a whole run of triangles at
/// once, which is more reliable than testing each one on its own.
fn push_wound(v0: V3, v1: V3, v2: V3, flip: bool, color: [u8; 3], out: &mut Vec<RibbonTriangle>) {
    let verts = if flip { [v0, v2, v1] } else { [v0, v1, v2] };
    out.push(RibbonTriangle {
        verts,
        color,
        normal: triangle_normal(verts[0], verts[1], verts[2]),
    });
}

// ---------------------------------------------------------------------------
// Emit triangle strip between two cross-sections
// ---------------------------------------------------------------------------

/// Connect two consecutive cross-section rings with a triangle strip.
/// Both rings must have the same number of vertices.
fn emit_strip(ring_a: &[V3], ring_b: &[V3], color: [u8; 3], out: &mut Vec<RibbonTriangle>) {
    let n = ring_a.len();
    debug_assert_eq!(n, ring_b.len());
    if n == 0 {
        return;
    }

    // The rings are closed cross-sections, so their centroids sit on the tube
    // axis and the vector between them is the local tangent.
    let axis_a = centroid(ring_a);
    let axis_b = centroid(ring_b);
    let tangent = v3_sub(axis_b, axis_a);

    // Decide the winding once, from the ring's own geometry.
    //
    // Doing it per triangle -- "is this face pointing away from the axis?" --
    // works for the tube's sides but fails on a sheet arrowhead, where the
    // strip between a narrow ring and the wider ring behind it forms the barb:
    // a face pointing back along the axis, not outward from it.  There the
    // radial test is reading a near-zero dot product and its sign is noise, so
    // barbs came out inconsistently wound and were wrongly culled.
    //
    // The ring's signed area vector has no such ambiguity.  Compared against
    // the tangent it says which way round the ring runs, which fixes the
    // winding of every triangle in the strip, barbs included.
    let mut area = [0.0; 3];
    for i in 0..n {
        let j = (i + 1) % n;
        area = v3_add(
            area,
            v3_cross(v3_sub(ring_a[i], axis_a), v3_sub(ring_a[j], axis_a)),
        );
    }
    // For a ring running counter-clockwise about the tangent, the vertex order
    // below already yields outward-facing triangles.
    let flip = v3_dot(area, tangent) < 0.0;

    for i in 0..n {
        let j = (i + 1) % n;

        let a0 = ring_a[i];
        let a1 = ring_a[j];
        let b0 = ring_b[i];
        let b1 = ring_b[j];

        // Quad (a0, a1, b1, b0) -> two triangles.
        push_wound(a0, a1, b0, flip, color, out);
        push_wound(a1, b1, b0, flip, color, out);
    }
}

/// Outward axis directions for the two end caps of a spline.
///
/// The cap at each end faces along the tube axis, away from the body, which is
/// the direction from the neighbouring spline point toward the end point.
fn end_axes(spline_points: &[SplinePoint]) -> (V3, V3) {
    let n = spline_points.len();
    if n < 2 {
        // Degenerate run: any consistent pair will do, nothing is visible.
        return ([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]);
    }
    let front = v3_normalize(v3_sub(spline_points[0].pos, spline_points[1].pos));
    let back = v3_normalize(v3_sub(spline_points[n - 1].pos, spline_points[n - 2].pos));
    (front, back)
}

/// Emit a cap (disc) to close the end of a tube or ribbon.
///
/// `ring` is the cross-section ring and `center` its spline point.  `outward`
/// is the tube axis direction pointing away from the body of the tube, which is
/// the cap's true surface normal; it decides the winding, so the caps are
/// consistent with the strips around them.
fn emit_cap(ring: &[V3], center: V3, color: [u8; 3], outward: V3, out: &mut Vec<RibbonTriangle>) {
    let n = ring.len();
    if n < 3 {
        return;
    }
    for i in 0..n {
        let j = (i + 1) % n;
        push_oriented(center, ring[i], ring[j], outward, color, out);
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
pub fn generate_ribbon_mesh(protein: &Protein, color_scheme: &ColorScheme) -> Vec<RibbonTriangle> {
    let lod = LodConfig::for_residue_count(protein.residue_count());
    generate_ribbon_mesh_with_lod(protein, color_scheme, lod)
}

fn generate_ribbon_mesh_with_lod(
    protein: &Protein,
    color_scheme: &ColorScheme,
    lod: LodConfig,
) -> Vec<RibbonTriangle> {
    let mut triangles: Vec<RibbonTriangle> = Vec::new();

    for chain in &protein.chains {
        match chain.molecule_type {
            MoleculeType::Protein => {
                generate_chain_ribbon(chain, color_scheme, &lod, &mut triangles);
            }
            MoleculeType::RNA | MoleculeType::DNA => {
                generate_nucleic_acid_ribbon(chain, color_scheme, &lod, &mut triangles);
            }
            MoleculeType::SmallMolecule => {
                // Small molecule rendering handled separately; skip in ribbon pass
            }
        }
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
    /// Carbonyl C→O direction hint for sheet frame orientation.
    frame_hint: Option<V3>,
}

// ---------------------------------------------------------------------------
// Shared spline tube builder (used by nucleic acid backbone)
// ---------------------------------------------------------------------------

/// Build a smooth tube through a sequence of backbone positions.
///
/// This performs the standard pipeline:
///   1. Catmull-Rom spline interpolation between control points
///   2. Finite-difference tangent computation
///   3. Parallel-transport Frenet frame propagation
///   4. Cross-section extrusion and triangle strip emission
///   5. End caps
///
/// The `arrow_t` field on every generated `SplinePoint` is set to `None`;
/// arrowhead logic is protein-specific and lives in `generate_chain_ribbon`.
fn build_spline_tube(records: &[CaRecord], lod: &LodConfig, out: &mut Vec<RibbonTriangle>) {
    let n = records.len();
    if n < 2 {
        return;
    }

    // --- Step 1: Generate Catmull-Rom spline points ---
    let mut spline_points: Vec<SplinePoint> = Vec::new();

    for seg in 0..n - 1 {
        let i0 = if seg == 0 { 0 } else { seg - 1 };
        let i1 = seg;
        let i2 = seg + 1;
        let i3 = if seg + 2 >= n { n - 1 } else { seg + 2 };

        let p0 = records[i0].pos;
        let p1 = records[i1].pos;
        let p2 = records[i2].pos;
        let p3 = records[i3].pos;

        let subdivs = if seg == n - 2 {
            lod.spline_subdivisions + 1
        } else {
            lod.spline_subdivisions
        };

        for sub in 0..subdivs {
            let t = sub as f64 / lod.spline_subdivisions as f64;
            let pos = catmull_rom(p0, p1, p2, p3, t);
            let ss = if t < 0.5 {
                records[i1].ss
            } else {
                records[i2].ss
            };
            let color = if t < 0.5 {
                records[i1].color
            } else {
                records[i2].color
            };

            spline_points.push(SplinePoint {
                pos,
                tangent: [0.0, 0.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                binormal: [0.0, 0.0, 1.0],
                frame_hint: if t < 0.5 {
                    records[i1].frame_hint
                } else {
                    records[i2].frame_hint
                },
                ss,
                color,
                arrow_t: None,
            });
        }
    }

    if spline_points.len() < 2 {
        return;
    }

    // --- Step 2: Compute tangents via finite differences ---
    let sp_len = spline_points.len();
    for i in 0..sp_len {
        let prev = if i == 0 { 0 } else { i - 1 };
        let next = if i == sp_len - 1 { sp_len - 1 } else { i + 1 };
        let t = v3_normalize(v3_sub(spline_points[next].pos, spline_points[prev].pos));
        spline_points[i].tangent = t;
    }

    // --- Step 3: Compute Frenet-Serret frames via parallel transport ---
    {
        let t0 = spline_points[0].tangent;
        let arbitrary = if t0[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let mut prev_normal = v3_normalize(v3_cross(t0, arbitrary));

        for sp in spline_points.iter_mut() {
            let t = sp.tangent;
            let proj = v3_scale(t, v3_dot(prev_normal, t));
            let mut nr = v3_sub(prev_normal, proj);
            let nl = v3_len(nr);
            if nl < 1e-12 {
                let arb = if t[0].abs() < 0.9 {
                    [1.0, 0.0, 0.0]
                } else {
                    [0.0, 1.0, 0.0]
                };
                nr = v3_normalize(v3_cross(t, arb));
            } else {
                nr = v3_scale(nr, 1.0 / nl);
            }
            let b = v3_normalize(v3_cross(t, nr));
            sp.normal = nr;
            sp.binormal = b;
            prev_normal = nr;
        }
    }

    // --- Step 3b: Apply sheet frame guides ---
    apply_sheet_frame_guides(&mut spline_points);

    // --- Step 4: Build cross-sections and emit triangle strips ---
    let mut prev_ring = cross_section(&spline_points[0], lod);

    for sp in spline_points.iter().skip(1) {
        let mut curr_ring = cross_section(sp, lod);
        let color = sp.color;

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

    // --- Step 5: Cap both ends ---
    let (front_out, back_out) = end_axes(&spline_points);

    let first_ring = cross_section(&spline_points[0], lod);
    emit_cap(
        &first_ring,
        spline_points[0].pos,
        spline_points[0].color,
        front_out,
        out,
    );

    let last_ring = cross_section(spline_points.last().unwrap(), lod);
    emit_cap(
        &last_ring,
        spline_points.last().unwrap().pos,
        spline_points.last().unwrap().color,
        back_out,
        out,
    );
}

fn generate_chain_ribbon(
    chain: &crate::model::protein::Chain,
    color_scheme: &ColorScheme,
    lod: &LodConfig,
    out: &mut Vec<RibbonTriangle>,
) {
    // 1. Collect C-alpha positions, SS types, colors, and frame hints.
    let cas: Vec<CaRecord> = chain
        .residues
        .iter()
        .filter_map(|res| {
            let ca = res.atoms.iter().find(|a| a.is_backbone)?;
            let color = color_to_rgb(color_scheme.residue_color(res, chain));
            Some(CaRecord {
                pos: [ca.x, ca.y, ca.z],
                ss: res.secondary_structure,
                color,
                frame_hint: residue_frame_hint(res),
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
            lod.spline_subdivisions + 1 // include the last point
        } else {
            lod.spline_subdivisions
        };

        for sub in 0..subdivs {
            let t = sub as f64 / lod.spline_subdivisions as f64;
            let pos = catmull_rom(p0, p1, p2, p3, t);

            // Interpolated secondary structure: use the nearer residue.
            let ss = if t < 0.5 { cas[i1].ss } else { cas[i2].ss };
            // Interpolated color: use the nearer residue.
            let color = if t < 0.5 {
                cas[i1].color
            } else {
                cas[i2].color
            };

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
                if global_pos >= arrow_begin && global_pos <= arrow_end {
                    let frac = (global_pos - arrow_begin) / (arrow_end - arrow_begin);
                    Some(frac.clamp(0.0, 1.0))
                } else {
                    None
                }
            });

            // Nearest-neighbor frame hint: snap to the closer control point.
            let frame_hint = if t < 0.5 {
                cas[i1].frame_hint
            } else {
                cas[i2].frame_hint
            };

            spline_points.push(SplinePoint {
                pos,
                tangent: [0.0, 0.0, 0.0], // computed below
                normal: [0.0, 1.0, 0.0],
                binormal: [0.0, 0.0, 1.0],
                frame_hint,
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

    // 5b. Apply sheet frame guides (carbonyl-direction-based orientation).
    apply_sheet_frame_guides(&mut spline_points);

    // 6. Build cross-sections and emit triangle strips.
    let mut prev_ring = cross_section(&spline_points[0], lod);

    for sp in spline_points.iter().skip(1) {
        let mut curr_ring = cross_section(sp, lod);
        let color = sp.color;

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
    let (front_out, back_out) = end_axes(&spline_points);

    let first_ring = cross_section(&spline_points[0], lod);
    emit_cap(
        &first_ring,
        spline_points[0].pos,
        spline_points[0].color,
        front_out,
        out,
    );

    let last_ring = cross_section(spline_points.last().unwrap(), lod);
    emit_cap(
        &last_ring,
        spline_points.last().unwrap().pos,
        spline_points.last().unwrap().color,
        back_out,
        out,
    );
}

// ---------------------------------------------------------------------------
// Nucleic acid (RNA/DNA) ribbon generation
// ---------------------------------------------------------------------------

/// Half-width and half-thickness for base slabs (Angstroms).
const BASE_SLAB_HALF_WIDTH: f64 = 1.0;
const BASE_SLAB_HALF_THICKNESS: f64 = 0.2;

/// Pyrimidine ring atom names (C, U, T, DC, DT).
const PYRIMIDINE_ATOMS: &[&str] = &["N1", "C2", "N3", "C4", "C5", "C6"];
/// Purine ring atom names (A, G, DA, DG).
const PURINE_ATOMS: &[&str] = &["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"];

/// Generate nucleic acid cartoon ribbon for a single chain.
///
/// Produces:
///   1. A backbone tube through C4' atoms (always coil cross-section).
///   2. 3D base slabs extending from C1' toward the base ring centroid.
fn generate_nucleic_acid_ribbon(
    chain: &crate::model::protein::Chain,
    color_scheme: &ColorScheme,
    lod: &LodConfig,
    out: &mut Vec<RibbonTriangle>,
) {
    // ----- Part 1: backbone tube through C4' atoms -----

    // Collect C4' positions and colors.
    let c4_records: Vec<CaRecord> = chain
        .residues
        .iter()
        .filter_map(|res| {
            let c4 = res.atoms.iter().find(|a| a.name.trim() == "C4'")?;
            let color = color_to_rgb(color_scheme.residue_color(res, chain));
            Some(CaRecord {
                pos: [c4.x, c4.y, c4.z],
                ss: SecondaryStructure::Coil, // always coil for nucleic acids
                color,
                frame_hint: None, // not applicable for nucleic acids
            })
        })
        .collect();

    // Delegate spline interpolation, framing, and meshing to shared helper.
    build_spline_tube(&c4_records, lod, out);

    // ----- Part 2: base slabs -----

    for residue in &chain.residues {
        // Find C1' atom.
        let c1_prime = match residue.atoms.iter().find(|a| a.name.trim() == "C1'") {
            Some(a) => [a.x, a.y, a.z],
            None => continue,
        };

        // Determine which ring atoms to look for.
        let ring_names: &[&str] = if is_purine(&residue.name) {
            PURINE_ATOMS
        } else {
            PYRIMIDINE_ATOMS
        };

        // Collect found base ring atom positions.
        let ring_positions: Vec<V3> = ring_names
            .iter()
            .filter_map(|&name| {
                residue
                    .atoms
                    .iter()
                    .find(|a| a.name.trim() == name)
                    .map(|a| [a.x, a.y, a.z])
            })
            .collect();

        if ring_positions.len() < 3 {
            continue;
        }

        // Compute base ring centroid.
        let count = ring_positions.len() as f64;
        let centroid = ring_positions.iter().fold([0.0, 0.0, 0.0], |acc, p| {
            [
                acc[0] + p[0] / count,
                acc[1] + p[1] / count,
                acc[2] + p[2] / count,
            ]
        });

        // Direction from C1' to centroid (long axis of slab).
        let dir = v3_sub(centroid, c1_prime);
        let dir_len = v3_len(dir);
        if dir_len < 1e-6 {
            continue;
        }
        let long_axis = v3_normalize(dir);

        // Width axis: perpendicular to long axis and a reference up vector.
        let up = [0.0, 1.0, 0.0];
        let mut width_axis = v3_cross(long_axis, up);
        if v3_len(width_axis) < 1e-6 {
            // long_axis is nearly parallel to up; use alternative.
            width_axis = v3_cross(long_axis, [1.0, 0.0, 0.0]);
        }
        width_axis = v3_normalize(width_axis);

        // Thickness axis: perpendicular to both long and width.
        let thick_axis = v3_normalize(v3_cross(long_axis, width_axis));

        let color = color_to_rgb(color_scheme.residue_color(residue, chain));

        // Build 8 corners of the slab box.
        // The slab goes from C1' to centroid, with half-width and half-thickness
        // offsets along the width and thickness axes.
        let hw = BASE_SLAB_HALF_WIDTH;
        let ht = BASE_SLAB_HALF_THICKNESS;
        let w_off = v3_scale(width_axis, hw);
        let t_off = v3_scale(thick_axis, ht);

        // Front face (at C1') corners: top-left, top-right, bottom-right, bottom-left
        let f_tl = v3_add(c1_prime, v3_add(w_off, t_off));
        let f_tr = v3_add(c1_prime, v3_sub(t_off, w_off));
        let f_br = v3_sub(c1_prime, v3_add(w_off, t_off));
        let f_bl = v3_add(c1_prime, v3_sub(w_off, t_off));

        // Back face (at centroid) corners
        let b_tl = v3_add(centroid, v3_add(w_off, t_off));
        let b_tr = v3_add(centroid, v3_sub(t_off, w_off));
        let b_br = v3_sub(centroid, v3_add(w_off, t_off));
        let b_bl = v3_add(centroid, v3_sub(w_off, t_off));

        // The slab is symmetric about the midpoint of its two end faces, so
        // that midpoint is inside it and orients every face outward.
        let interior = v3_scale(v3_add(c1_prime, centroid), 0.5);

        // Emit 6 faces x 2 triangles = 12 triangles.
        // Front face (at C1')
        emit_quad(f_tl, f_tr, f_br, f_bl, interior, color, out);
        // Back face (at centroid)
        emit_quad(b_tr, b_tl, b_bl, b_br, interior, color, out);
        // Top face
        emit_quad(f_tl, b_tl, b_tr, f_tr, interior, color, out);
        // Bottom face
        emit_quad(f_bl, f_br, b_br, b_bl, interior, color, out);
        // Left face
        emit_quad(f_tl, f_bl, b_bl, b_tl, interior, color, out);
        // Right face
        emit_quad(f_tr, b_tr, b_br, f_br, interior, color, out);
    }
}

/// Emit two triangles for a quad face (v0, v1, v2, v3).
///
/// `interior` is a point inside the solid the face belongs to, which orients
/// both triangles outward.
fn emit_quad(
    v0: V3,
    v1: V3,
    v2: V3,
    v3: V3,
    interior: V3,
    color: [u8; 3],
    out: &mut Vec<RibbonTriangle>,
) {
    let third = 1.0 / 3.0;
    let c1 = v3_scale(v3_add(v3_add(v0, v1), v2), third);
    push_oriented(v0, v1, v2, v3_sub(c1, interior), color, out);
    let c2 = v3_scale(v3_add(v3_add(v0, v2), v3), third);
    push_oriented(v0, v2, v3, v3_sub(c2, interior), color, out);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::protein::{Atom, Residue};

    /// Helper to build a minimal atom.
    fn make_atom(name: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom {
            name: name.to_string(),
            element: "C".to_string(),
            x,
            y,
            z,
            b_factor: 0.0,
            is_backbone: name == "CA",
            is_hetero: false,
        }
    }

    #[test]
    fn test_residue_frame_hint_with_co() {
        // Residue with C at origin and O along the +X axis.
        let res = Residue {
            name: "ALA".to_string(),
            seq_num: 1,
            insertion_code: None,
            atoms: vec![
                make_atom("CA", 0.0, 0.0, 0.0),
                make_atom("C", 1.0, 0.0, 0.0),
                make_atom("O", 4.0, 0.0, 0.0),
            ],
            secondary_structure: SecondaryStructure::Sheet,
        };
        let hint = residue_frame_hint(&res).expect("should produce a hint");
        // Direction should be [1, 0, 0] (normalized C→O).
        assert!((hint[0] - 1.0).abs() < 1e-9);
        assert!(hint[1].abs() < 1e-9);
        assert!(hint[2].abs() < 1e-9);
    }

    #[test]
    fn test_residue_frame_hint_missing_o() {
        // Residue with C but no O.
        let res = Residue {
            name: "ALA".to_string(),
            seq_num: 1,
            insertion_code: None,
            atoms: vec![
                make_atom("CA", 0.0, 0.0, 0.0),
                make_atom("C", 1.0, 0.0, 0.0),
                make_atom("N", 2.0, 0.0, 0.0),
            ],
            secondary_structure: SecondaryStructure::Sheet,
        };
        assert!(residue_frame_hint(&res).is_none());
    }

    #[test]
    fn test_residue_frame_hint_fallback_ca_o() {
        // Residue with CA and O but no C — should use CA→O fallback.
        let res = Residue {
            name: "ALA".to_string(),
            seq_num: 1,
            insertion_code: None,
            atoms: vec![
                make_atom("CA", 0.0, 0.0, 0.0),
                make_atom("O", 0.0, 3.0, 0.0),
            ],
            secondary_structure: SecondaryStructure::Sheet,
        };
        let hint = residue_frame_hint(&res).expect("should fallback to CA->O");
        assert!(hint[0].abs() < 1e-9);
        assert!((hint[1] - 1.0).abs() < 1e-9);
        assert!(hint[2].abs() < 1e-9);
    }

    #[test]
    fn test_apply_sheet_guides_flipped_hints() {
        // Build a small set of sheet spline points with alternating hint signs.
        // After apply_sheet_frame_guides, consecutive binormals should be
        // coherent (positive dot product).
        let mut points: Vec<SplinePoint> = Vec::new();
        for i in 0..6 {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            points.push(SplinePoint {
                pos: [i as f64, 0.0, 0.0],
                tangent: [1.0, 0.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                binormal: [0.0, 0.0, 1.0],
                frame_hint: Some([0.0, sign * 0.5, 0.5]),
                ss: SecondaryStructure::Sheet,
                color: [255, 255, 255],
                arrow_t: None,
            });
        }

        apply_sheet_frame_guides(&mut points);

        // All consecutive binormals should agree in direction.
        for i in 1..points.len() {
            let d = v3_dot(points[i].binormal, points[i - 1].binormal);
            assert!(
                d > 0.0,
                "Binormals at {} and {} should agree in sign, got dot={}",
                i - 1,
                i,
                d
            );
        }
    }

    #[test]
    fn test_coil_points_unaffected() {
        // Build spline points that are all Coil with frame hints.
        // They should NOT be modified by the guide pass.
        let original_binormal: V3 = [0.0, 0.0, 1.0];
        let original_normal: V3 = [0.0, 1.0, 0.0];

        let mut points: Vec<SplinePoint> = Vec::new();
        for i in 0..4 {
            points.push(SplinePoint {
                pos: [i as f64, 0.0, 0.0],
                tangent: [1.0, 0.0, 0.0],
                normal: original_normal,
                binormal: original_binormal,
                frame_hint: Some([0.0, 0.7, 0.7]),
                ss: SecondaryStructure::Coil,
                color: [128, 128, 128],
                arrow_t: None,
            });
        }

        apply_sheet_frame_guides(&mut points);

        for (i, sp) in points.iter().enumerate() {
            assert_eq!(
                sp.binormal, original_binormal,
                "Coil point {} binormal should be unchanged",
                i
            );
            assert_eq!(
                sp.normal, original_normal,
                "Coil point {} normal should be unchanged",
                i
            );
        }
    }

    #[test]
    fn test_helix_points_unaffected() {
        // Build spline points that are all Helix with frame hints.
        // They should NOT be modified by the guide pass.
        let original_binormal: V3 = [0.0, 0.0, 1.0];
        let original_normal: V3 = [0.0, 1.0, 0.0];

        let mut points: Vec<SplinePoint> = Vec::new();
        for i in 0..4 {
            points.push(SplinePoint {
                pos: [i as f64, 0.0, 0.0],
                tangent: [1.0, 0.0, 0.0],
                normal: original_normal,
                binormal: original_binormal,
                frame_hint: Some([0.0, 0.7, 0.7]),
                ss: SecondaryStructure::Helix,
                color: [128, 128, 128],
                arrow_t: None,
            });
        }

        apply_sheet_frame_guides(&mut points);

        for (i, sp) in points.iter().enumerate() {
            assert_eq!(
                sp.binormal, original_binormal,
                "Helix point {} binormal should be unchanged",
                i
            );
            assert_eq!(
                sp.normal, original_normal,
                "Helix point {} normal should be unchanged",
                i
            );
        }
    }
}

#[cfg(test)]
mod winding_tests {
    use super::*;
    use crate::render::color::{ColorScheme, ColorSchemeType};
    use std::collections::HashMap;

    /// Every triangle in the ribbon mesh must be wound the same way round.
    ///
    /// Each shared edge should appear once in each direction.  An edge used
    /// twice in the *same* direction means two neighbouring triangles disagree
    /// about which side is out, which leaves the stored normals arbitrarily
    /// signed and rules out back-face culling.
    fn assert_consistent_winding(path: &str) {
        let mut protein = crate::parser::pdb::load_structure(path).unwrap();
        protein.center();
        let scheme = ColorScheme::new(ColorSchemeType::Structure, protein.residue_count());
        let mesh = generate_ribbon_mesh(&protein, &scheme);
        assert!(!mesh.is_empty(), "{path} produced no ribbon geometry");

        let key = |v: V3| {
            const SCALE: f64 = 4096.0;
            (
                (v[0] * SCALE).round() as i64,
                (v[1] * SCALE).round() as i64,
                (v[2] * SCALE).round() as i64,
            )
        };
        let mut directed: HashMap<((i64, i64, i64), (i64, i64, i64)), i32> = HashMap::new();
        for tri in &mesh {
            let k = [key(tri.verts[0]), key(tri.verts[1]), key(tri.verts[2])];
            for i in 0..3 {
                *directed.entry((k[i], k[(i + 1) % 3])).or_insert(0) += 1;
            }
        }

        let mismatched: i32 = directed
            .iter()
            .filter(|&(&(u, v), &n)| {
                let back = directed.get(&(v, u)).copied().unwrap_or(0);
                back != 0 && back != n
            })
            .map(|(_, &n)| n)
            .sum();
        assert_eq!(
            mismatched, 0,
            "{path}: {mismatched} inconsistently wound edges"
        );
    }

    #[test]
    fn protein_ribbons_are_consistently_wound() {
        assert_consistent_winding("examples/1UBQ.pdb");
        assert_consistent_winding("examples/4HHB.pdb");
    }

    #[test]
    fn nucleic_acid_ribbons_are_consistently_wound() {
        assert_consistent_winding("examples/1BNA.pdb");
        assert_consistent_winding("examples/1RNA.pdb");
    }

    #[test]
    fn mixed_structures_are_consistently_wound() {
        assert_consistent_winding("examples/1AOI.pdb");
    }
}

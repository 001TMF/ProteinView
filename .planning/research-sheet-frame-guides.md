# Research: Beta Sheet Frame Guides (Carbonyl-Direction-Based Ribbon Orientation)

## 1. Algorithm: How Carbonyl Direction Guides Work

### 1.1 The Problem

Our current ribbon renderer computes local coordinate frames using **parallel transport**: an initial normal is chosen perpendicular to the first tangent, then propagated along the spline by projecting out the tangent component at each step. This produces a smooth, twist-minimizing frame that works well for helices and coils, but fails for beta sheets.

In a real beta sheet, the peptide planes alternate orientation in a pleated pattern. The flat ribbon should align with the peptide plane so the viewer sees the characteristic flat, wide appearance of a sheet. Parallel transport is ignorant of the physical peptide geometry and can produce arbitrary rotations of the ribbon cross-section within sheet regions -- the ribbon may appear edge-on when it should be face-on, or twist gradually along a strand.

Molecular graphics programs (PyMOL, ChimeraX, VMD) solve this by using the **carbonyl (C=O) bond direction** as a frame guide for sheet residues.

### 1.2 `residue_frame_hint()` -- Extracting the C->O Direction

The fork's implementation (lines 666-680 of `desperadus-master:src/render/ribbon.rs`):

```rust
fn residue_frame_hint(residue: &Residue) -> Option<V3> {
    let carbonyl = match (atom_pos(residue, "C"), atom_pos(residue, "O")) {
        (Some(c), Some(o)) => Some(v3_sub(o, c)),
        _ => None,
    };
    let ca_to_o = match (atom_pos(residue, "CA"), atom_pos(residue, "O")) {
        (Some(ca), Some(o)) => Some(v3_sub(o, ca)),
        _ => None,
    };

    carbonyl.or(ca_to_o).and_then(|hint| {
        let len = v3_len(hint);
        (len >= 1e-8).then(|| v3_scale(hint, 1.0 / len))
    })
}
```

**How it works:**
1. **Primary**: Compute `O - C` (the carbonyl bond direction vector). This is the most physically meaningful direction -- it points perpendicular to the backbone within the peptide plane.
2. **Fallback**: If the backbone C atom is missing but CA and O exist, use `O - CA` instead. This is less accurate but still roughly perpendicular to the strand direction.
3. **Normalization**: The result is normalized to a unit vector. If the vector length is degenerate (< 1e-8), returns `None` to indicate no hint is available.
4. **Return type**: `Option<V3>` -- `None` when atoms are missing or degenerate.

The helper `atom_pos()` simply searches the residue's atom list by name:
```rust
fn atom_pos(residue: &Residue, name: &str) -> Option<V3> {
    residue.atoms.iter()
        .find(|atom| atom.name.trim() == name)
        .map(|atom| [atom.x, atom.y, atom.z])
}
```

### 1.3 `apply_sheet_frame_guides()` -- The Core Algorithm

This function runs as a post-processing pass after parallel transport frames are computed. It modifies only the frames within sheet segments.

**Step 1: Identify contiguous sheet runs.**
```rust
let mut i = 0;
while i < spline_points.len() {
    if spline_points[i].ss != SecondaryStructure::Sheet {
        i += 1;
        continue;
    }
    let run_start = i;
    while i < spline_points.len() && spline_points[i].ss == SecondaryStructure::Sheet {
        i += 1;
    }
    let run_end = i;
    // ... process this run
}
```

**Step 2: For each spline point in the run, compute a guided binormal.**

For each spline point at index `idx`, the algorithm looks at a neighborhood of +/-2 spline points (clamped to the run boundaries). It collects the `frame_hint` from each neighbor that has one, projects each hint perpendicular to the local tangent, and aligns the sign of each projected hint for consistency:

```rust
for neighbor_idx in idx.saturating_sub(2)..=(idx + 2).min(run_end - 1) {
    let hint = spline_points[neighbor_idx].frame_hint;
    let projected = project_perpendicular(hint, tangent);
    let aligned = align_vector_sign(projected, reference);
    acc = v3_add(acc, aligned);
    count += 1;
}
guided_binormals[local_idx] = Some(v3_normalize(acc));
```

The first hint in a run is seeded against the parallel-transported binormal (`base`) to establish a consistent sign convention. Subsequent hints align against the previous accepted hint (`guide_ref`).

**Step 3: Blend guided binormal with parallel-transported binormal.**

```rust
let blended = v3_normalize(
    v3_add(v3_scale(transported, 0.35), v3_scale(guided, 0.65))
);
```

The **65/35 blend ratio** means:
- 65% weight on the carbonyl-guided direction (physical reality)
- 35% weight on the parallel-transported direction (geometric smoothness)

This blend prevents the ribbon from snapping abruptly to the guide direction and maintains some of the smooth propagation behavior. A 100% guide weight would produce a more physically accurate but potentially noisier ribbon.

**Step 4: Ensure sign consistency along the run.**

After blending, the final binormal is sign-aligned against the previous point's binormal:
```rust
let final_binormal = align_vector_sign(blended, prev_binormal);
let final_normal = align_vector_sign(
    v3_normalize(v3_cross(final_binormal, tangent)),
    prev_normal
);
```

The normal is recomputed from `binormal x tangent` and also sign-aligned. Both `prev_binormal` and `prev_normal` are propagated forward to maintain a consistent orientation along the entire sheet run.

### 1.4 `project_perpendicular()` and `align_vector_sign()`

**`project_perpendicular(vec, tangent)`**: Removes the tangent component from `vec`, leaving only the component perpendicular to the backbone direction. Returns `None` if the result is degenerate (the hint is nearly parallel to the tangent, which can happen at chain termini).

```rust
fn project_perpendicular(vec: V3, tangent: V3) -> Option<V3> {
    let projected = v3_sub(vec, v3_scale(tangent, v3_dot(vec, tangent)));
    let len = v3_len(projected);
    (len >= 1e-8).then(|| v3_scale(projected, 1.0 / len))
}
```

**`align_vector_sign(vec, reference)`**: Flips `vec` if it points away from `reference` (negative dot product). This resolves the fundamental ambiguity that the C=O direction in a beta sheet alternates between adjacent residues (the pleating), and we want a consistent binormal direction, not one that flips every residue.

```rust
fn align_vector_sign(vec: V3, reference: V3) -> V3 {
    if v3_dot(vec, reference) < 0.0 {
        v3_scale(vec, -1.0)
    } else {
        vec
    }
}
```

### 1.5 Why Sign Alignment Matters

In a beta sheet, the C=O bonds of successive residues point in **alternating directions** due to the pleated geometry. If we naively used the raw C=O direction as the binormal, the ribbon would flip 180 degrees at every residue boundary, producing a zig-zag surface instead of a flat sheet. The sign alignment step resolves this by ensuring all guide vectors point in the same half-space, choosing the half-space consistent with the initial parallel-transported frame.

---

## 2. Integration Points in Our Codebase

### 2.1 Where `CaRecord` is Constructed

**File**: `src/render/ribbon.rs`, lines 513-527 (`generate_chain_ribbon`)

Currently:
```rust
let cas: Vec<CaRecord> = chain.residues.iter().filter_map(|res| {
    let ca = res.atoms.iter().find(|a| a.is_backbone)?;
    let color = color_to_rgb(color_scheme.residue_color(res, chain));
    Some(CaRecord {
        pos: [ca.x, ca.y, ca.z],
        ss: res.secondary_structure,
        color,
    })
}).collect();
```

**Change needed**: Add `frame_hint: Option<V3>` to the `CaRecord` struct (line 353-357) and populate it via `residue_frame_hint(res)` in the construction closure. The `CaRecord` struct definition currently has three fields (`pos`, `ss`, `color`); a fourth `frame_hint` field adds negligible memory.

### 2.2 Where `compute_frames()` Would Be Called

Our codebase does **not** currently have a standalone `compute_frames()` function. The tangent computation and parallel transport are done inline in `generate_chain_ribbon()` (lines 637-684) and `build_spline_tube()` (lines 427-465).

**Change needed**: Extract our inline tangent + parallel transport code into a `compute_frames()` function (as the fork does). Then add the `apply_sheet_frame_guides()` call at the end of `compute_frames()`. This refactor also benefits `build_spline_tube()` which has identical logic.

The fork's `compute_frames()` function (lines 427-470) does exactly three things:
1. Finite-difference tangent computation (same as our step 4)
2. Parallel transport normal propagation (same as our step 5)
3. `apply_sheet_frame_guides()` as a post-processing pass

### 2.3 Where `SplinePoint` is Constructed

**File**: `src/render/ribbon.rs`, lines 620-628

Currently:
```rust
spline_points.push(SplinePoint {
    pos,
    tangent: [0.0, 0.0, 0.0],
    normal: [0.0, 1.0, 0.0],
    binormal: [0.0, 0.0, 1.0],
    ss,
    color,
    arrow_t,
});
```

**Change needed**: Add `frame_hint: Option<V3>` to `SplinePoint` struct (line 129) and populate it during spline construction. The fork uses nearest-residue assignment:
```rust
frame_hint: if t < 0.5 { cas[i1].frame_hint } else { cas[i2].frame_hint },
```

This is a **snap** interpolation (nearest-neighbor), not a smooth interpolation. This is intentional: the hint is a **direction** not a position, and linearly interpolating two unit vectors that might have opposite signs can produce a zero-length vector at the midpoint. The snap approach avoids this pitfall entirely.

### 2.4 Data Flow Summary

```
Residue (has C, O atoms)
  |
  v  residue_frame_hint(residue) -> Option<V3>
  |
CaRecord { pos, frame_hint, ss, color }
  |
  v  Catmull-Rom spline interpolation (position is interpolated, hint is snapped)
  |
SplinePoint { pos, tangent, normal, binormal, frame_hint, ss, color, arrow_t }
  |
  v  compute_frames():
  |    1. Finite-difference tangents
  |    2. Parallel transport normals/binormals
  |    3. apply_sheet_frame_guides() post-pass
  |
SplinePoint (frames now corrected for sheets)
  |
  v  cross_section() -> triangle mesh
```

---

## 3. Pitfall Analysis

### 3.1 Missing C or O Atoms

**Risk**: Some residues lack the backbone C or O atom. This can happen for:
- Terminal residues (especially C-terminal where OXT may replace O)
- Incomplete/truncated structures
- Non-standard residues (MSE for selenomethionine, etc.)
- Low-resolution structures with missing atoms

**Fork's mitigation**: `residue_frame_hint()` returns `Option<V3>`, and `apply_sheet_frame_guides()` simply skips spline points with `frame_hint: None`. The averaging over a +/-2 neighborhood means a single missing hint among 5 neighbors still produces a guided frame. Only if ALL neighbors in the window lack hints does the algorithm fall back to pure parallel transport.

**Our additional consideration**: We should verify that pdbtbx preserves C and O atoms in the atom list (confirmed: our parser at `src/parser/pdb.rs` line 37-46 maps ALL atoms from each residue, not just CA). The `is_backbone` flag is only true for CA and C4', but C and O atoms are present in `residue.atoms` with `is_backbone: false`.

### 3.2 Coil Regions Between Sheets

**Risk**: The frame guide should only apply to sheet (and potentially helix) regions. Applying it to coils would be wrong -- coil regions have no consistent peptide plane orientation relative to the backbone direction.

**Fork's mitigation**: `apply_sheet_frame_guides()` explicitly checks `spline_points[i].ss != SecondaryStructure::Sheet` and skips non-sheet points. The function identifies contiguous sheet runs and only modifies frames within those runs. Coil and helix frames remain purely parallel-transported.

**Boundary behavior**: At the transition between sheet and coil, there may be a small discontinuity in the frame. The 35% parallel-transport weight in the blend helps smooth this. The first and last points of a sheet run will have fewer neighbors in the averaging window (the +/-2 window is clamped to the run boundaries), which naturally tapers the guide influence at boundaries.

### 3.3 Hint Interpolation Through the Spline

**Risk**: The fork uses snap (nearest-neighbor) interpolation for hints rather than smooth interpolation. This means the hint value changes discontinuously at `t = 0.5` between two residues.

**Why this is acceptable**:
- The hints are only used as a guide direction for orientation, not for geometry.
- Linear interpolation of two unit vectors with potentially opposite signs can produce a zero-length vector at the midpoint (a classic problem with quaternion interpolation).
- The +/-2 neighbor averaging in `apply_sheet_frame_guides()` already provides smoothing.
- The 65/35 blend with parallel transport provides additional smoothness.

**Alternative considered**: One could use slerp (spherical linear interpolation) for the hints through the spline, but this adds complexity and requires sign-coherent hints at control points, which is exactly the problem the post-pass solves.

### 3.4 Performance Impact

**Per-residue cost**: `residue_frame_hint()` performs two linear searches through the residue's atom list (for "C" and "O"). Typical residues have ~5-20 atoms, so each search is O(n) with small n.

**Per-spline-point cost in `apply_sheet_frame_guides()`**: For each spline point in a sheet run, the algorithm looks at up to 5 neighbors (the point itself +/-2). Each neighbor requires a perpendicular projection (a few vector ops) and a sign alignment (one dot product). This is O(1) per spline point.

**Overall impact**: For a typical protein with ~300 residues, ~14 spline subdivisions each, and maybe 30% in sheet regions, this adds roughly:
- 300 calls to `residue_frame_hint()` (negligible: ~300 * 2 linear searches of ~10 atoms)
- ~1260 spline points in sheet regions, each doing ~5 neighbor lookups with vector math
- Total: a few thousand floating-point operations. **Negligible** compared to the triangle generation and rasterization cost.

### 3.5 Interaction with Mesh Caching

Our mesh cache (`app.rs` line 49: `pub mesh_cache: Vec<RibbonTriangle>`) is regenerated when the color scheme changes (line 156). The frame hints are derived from atom positions, which are static after loading. Adding frame guides does not change the cache invalidation strategy -- the hints will produce the same frames every time `generate_ribbon_mesh()` is called for a given protein structure.

---

## 4. Potential Improvements Over the Fork

### 4.1 Extend Hints to Helices

Alpha helices also have a canonical C=O direction (the hydrogen bond runs roughly along the helix axis, with C=O pointing ~parallel to the axis). However, in helices the parallel transport frame already works reasonably well because the backbone follows a smooth helical path. Applying carbonyl guides to helices could:

- **Pro**: Make the helix ribbon orientation more physically accurate (the flat face of the ribbon would align with the C=O direction, showing hydrogen bond directionality)
- **Con**: Helix ribbons already look correct with parallel transport; adding guides might introduce unnecessary complexity
- **Recommendation**: Start with sheets only (matching the fork). Consider helices as a future enhancement if visual artifacts are observed.

### 4.2 Peptide Plane Normal via Cross Product

Instead of using the raw C=O direction, we could compute the **peptide plane normal** as:

```
hint = normalize( (CA->C) x (CA->N) )
```

This cross product gives the normal to the plane defined by CA, C, and N atoms -- the peptide plane. This is arguably more robust because:
- It uses three atoms instead of two, reducing sensitivity to individual atom position errors
- It directly gives the plane normal rather than a vector within the plane
- It works even if the O atom is missing (as long as CA, C, N are present)

**However**: The fork's C->O approach is the standard used by PyMOL and other viewers. The C=O direction, after perpendicular projection, gives the same information as the peptide plane normal (they differ only by a 90-degree rotation within the plane perpendicular to the tangent, which the algorithm compensates for). Using C->O is simpler and matches established practice.

**Recommendation**: Use the C->O approach (matching the fork and standard practice). Consider the cross-product alternative as a fallback when O is missing but C and N are present.

### 4.3 RNA/DNA Chain Safety

Our `generate_nucleic_acid_ribbon()` already sets `frame_hint: None` for all nucleic acid CaRecords in the fork:
```rust
Some(CaRecord {
    pos: [c4.x, c4.y, c4.z],
    frame_hint: None,
    ss: SecondaryStructure::Coil, // always coil for nucleic acids
    color,
})
```

Nucleic acids are safe because:
1. Their SS is always `Coil`, so `apply_sheet_frame_guides()` will skip them entirely
2. Their `frame_hint` is `None`, so even if the SS check were removed, no guide would be applied
3. They don't have C=O carbonyl bonds in the backbone sense (nucleotide bases do have C=O but `residue_frame_hint()` looks for backbone "C" and "O" atom names, which nucleotides don't have)

**Our code already handles this correctly** through the `MoleculeType` dispatch in `generate_ribbon_mesh()` -- nucleic acid chains go through `generate_nucleic_acid_ribbon()` which uses `build_spline_tube()`, not `generate_chain_ribbon()`.

### 4.4 Fallback Hint: CA->C x CA->N Cross Product

As an enhancement beyond the fork, we could add a second fallback in `residue_frame_hint()`:

```rust
fn residue_frame_hint(residue: &Residue) -> Option<V3> {
    // Primary: C->O carbonyl direction
    let carbonyl = match (atom_pos(residue, "C"), atom_pos(residue, "O")) {
        (Some(c), Some(o)) => Some(v3_sub(o, c)),
        _ => None,
    };
    // Fallback 1: CA->O direction
    let ca_to_o = match (atom_pos(residue, "CA"), atom_pos(residue, "O")) {
        (Some(ca), Some(o)) => Some(v3_sub(o, ca)),
        _ => None,
    };
    // Fallback 2: peptide plane normal from CA, C, N
    let plane_normal = match (atom_pos(residue, "CA"), atom_pos(residue, "C"), atom_pos(residue, "N")) {
        (Some(ca), Some(c), Some(n)) => {
            let v1 = v3_sub(c, ca);
            let v2 = v3_sub(n, ca);
            Some(v3_cross(v1, v2))
        },
        _ => None,
    };

    carbonyl.or(ca_to_o).or(plane_normal).and_then(|hint| {
        let len = v3_len(hint);
        (len >= 1e-8).then(|| v3_scale(hint, 1.0 / len))
    })
}
```

**Note**: The plane normal is perpendicular to the C=O direction (it's the normal TO the plane, while C=O lies IN the plane). This means the sign alignment would still work, but the guide would orient the ribbon differently. We would need to verify that the final visual result is correct, or rotate the plane normal 90 degrees to match the C=O convention.

**Recommendation**: Start with the fork's approach (C->O primary, CA->O fallback). Add the plane normal fallback only if testing reveals residues that lack both C->O and CA->O but have CA, C, N.

---

## 5. Test Strategy

### 5.1 Unit Tests

**Test `residue_frame_hint()`:**
- Construct a `Residue` with known C, O, CA, N atom positions. Verify the returned hint direction matches the expected C->O unit vector.
- Construct a `Residue` missing the C atom but having CA and O. Verify fallback to CA->O direction.
- Construct a `Residue` missing both C and O. Verify `None` is returned.
- Construct a `Residue` where C and O are at the same position (degenerate). Verify `None` is returned.

**Test `project_perpendicular()`:**
- Verify that projecting a vector perpendicular to a tangent returns the same vector (up to normalization).
- Verify that projecting a vector parallel to the tangent returns `None`.
- Verify that an oblique vector is correctly decomposed.

**Test `align_vector_sign()`:**
- Verify that a vector in the same direction as reference is unchanged.
- Verify that a vector opposite to reference is negated.
- Verify behavior with orthogonal vectors (zero dot product -- currently returns the original vector, which is fine).

**Test `apply_sheet_frame_guides()` integration:**
- Construct a small set of `SplinePoint` values with known positions, tangents, and frame_hints, all marked as `Sheet`. Run `apply_sheet_frame_guides()` and verify:
  - The binormals are roughly aligned with the projected hints
  - The binormals are sign-consistent along the run
  - The normals, binormals, and tangents form a right-handed orthonormal frame

**Test that non-sheet regions are unaffected:**
- Construct a mixed sequence (Coil, Sheet, Coil). Run `compute_frames()`. Verify that the Coil regions have the same frames as they would without the guide pass.

### 5.2 Visual/Regression Tests

**Side-by-side comparison:**
- Render a protein with prominent beta sheets (e.g., immunoglobulin fold like 1ZVH, or a beta-barrel) with and without frame guides.
- Verify that sheet ribbons appear flat and face-on rather than edge-on.
- Verify that the ribbon does not flip or twist within a sheet strand.

**Specific test structures:**
- **1ZVH.cif** (already in our examples): Has both helices and sheets. Good for checking that sheets are improved while helices are unchanged.
- **1CRN** (crambin): Small protein with a beta sheet. Good for quick visual verification.
- **1TIM** (TIM barrel): 8 parallel beta strands surrounding 8 alpha helices. Excellent stress test for sheet guides because the parallel strands should all appear flat.

**Edge cases to test visually:**
- Very short sheet strands (1-2 residues): Verify no visual artifacts.
- Sheet strands at chain termini: Verify missing atom fallbacks work.
- Proteins with no sheets (e.g., all-alpha): Verify no regression.
- NMR structures (may have unusual geometry): Verify graceful handling.

### 5.3 Quantitative Frame Validation

For automated testing beyond visual inspection:
- For each spline point in a sheet region, verify that `|dot(normal, tangent)| < epsilon` and `|dot(binormal, tangent)| < epsilon` (orthogonality).
- Verify that `|dot(normal, binormal)| < epsilon` (orthogonality).
- Verify that `|len(normal) - 1.0| < epsilon` and `|len(binormal) - 1.0| < epsilon` (unit length).
- Verify that `cross(normal, binormal)` is approximately equal to `tangent` or `-tangent` (right-handed or left-handed frame, but consistent).

### 5.4 Performance Benchmarking

- Time `generate_ribbon_mesh()` before and after the change on a large protein (e.g., 1AOI with ~500 residues, or a larger structure).
- The expected overhead is < 1% of mesh generation time. If it exceeds 5%, investigate.

---

## 6. Implementation Plan (Summary)

The changes required are localized to `src/render/ribbon.rs`:

1. **Add `frame_hint: Option<V3>` to `CaRecord`** (struct definition, line ~353)
2. **Add `frame_hint: Option<V3>` to `SplinePoint`** (struct definition, line ~129)
3. **Add `atom_pos()` helper function** (~5 lines)
4. **Add `residue_frame_hint()` function** (~15 lines)
5. **Add `project_perpendicular()` function** (~5 lines)
6. **Add `align_vector_sign()` function** (~5 lines)
7. **Add `apply_sheet_frame_guides()` function** (~65 lines)
8. **Extract `compute_frames()` function** from inline code in `generate_chain_ribbon()` and `build_spline_tube()` (~40 lines, refactored from existing code)
9. **Wire up**: Populate `frame_hint` in `CaRecord` construction, propagate through spline interpolation, call `compute_frames()` which calls `apply_sheet_frame_guides()`
10. **Set `frame_hint: None`** in nucleic acid `CaRecord` construction (already the case in `build_spline_tube` since hints come from residues)

Total new code: approximately 100-120 lines. Total refactored code: approximately 40 lines extracted into `compute_frames()`. No changes needed outside `src/render/ribbon.rs`.

# Desperadus Fork Analysis

Analysis of the `desperadus-master` branch (Desperadus/ProteinView fork) compared
to our `master` branch. The fork contains 25 commits on top of an older version of
our codebase (pre-RNA/DNA, pre-small-molecule support).

---

## 1. Overview

The fork introduces several meaningful features: secondary structure inference from
backbone geometry (hydrogen bonding + torsion angles), pLDDT confidence coloring for
AlphaFold models, a proper orientation-matrix camera replacing Euler angles (fixing a
mirroring bug), beta-sheet rendering improvements using carbonyl direction guides,
adaptive spline tessellation, tiled parallel rasterization, cartoon edge outlines,
a chain-focus zoom feature, CI/CD for PyPI publishing via maturin, and ligand
rendering. It also includes a large volume of `ruff`-style autoformatting that
inflates the apparent diff size.

The fork **does not** have our RNA/DNA support, our `Ligand`/`LigandType` data model,
our separate `ligands` vec on `Protein`, our binding-pocket analysis, our NMR
multi-model fix, or our `is_hetero` atom field. The fork's ligand approach is
fundamentally different -- it stuffs ligands into chain residues and filters them via
`is_ligand_residue()` at render time rather than separating them at parse time.

---

## 2. Feature Inventory

### 2.1 Secondary Structure Inference from Coordinates

- **Description**: Falls back to inferring helices/sheets from backbone geometry when
  PDB HELIX/SHEET records and CIF SS annotations are absent (common with AlphaFold
  models). Uses DSSP-style hydrogen bond energy calculation (`E = 27.888 * (1/r_ON +
  1/r_CH - 1/r_OH - 1/r_CN)`) with amide H estimation, combined with phi/psi torsion
  angle windows. Assigns helices from i->i+3/4/5 H-bond patterns, sheets from
  parallel/antiparallel H-bond ladders. Includes gap-filling and minimum-run-length
  filtering.
- **Files**: `src/model/secondary.rs` (+350 lines), `src/parser/pdb.rs` (+3 lines)
- **Quality**: Clean. Well-structured with separate functions for torsion computation,
  H-bond mapping, helix/sheet assignment, gap filling, and run filtering. The H-bond
  energy formula matches DSSP literature. Torsion angle windows are reasonable
  (strong + weak variants for both helix and sheet). The `infer_protein_secondary_structure`
  function correctly skips chains that already have SS assignments and non-protein
  chains.
- **Tests**: Includes test against `AF3_TNFa.pdb` (AlphaFold model) and verification
  that explicit PDB SS records still take precedence (`1UBQ.pdb`).
- **Do we have equivalent?**: No. This is the highest-value feature in the fork.
  AlphaFold PDBs are a major use case and they look terrible without SS assignment.

### 2.2 pLDDT Confidence Coloring

- **Description**: Adds a `Plddt` color scheme variant that colors residues by
  AlphaFold confidence score stored in the B-factor column. Uses the standard AF
  four-band scheme: dark blue (>=90), light blue (>=70), yellow (>=50), orange (<50).
  Includes a `has_plddt()` heuristic on `Protein` that detects whether B-factors look
  like pLDDT scores (all in 0-100, mean >= 50, >= 25% above 70).
- **Files**: `src/render/color.rs` (+30 lines), `src/model/protein.rs` (+35 lines),
  `src/app.rs` (color cycle change), `src/main.rs` (CLI `--color plddt`)
- **Quality**: Clean. The heuristic is reasonable and conservative. The color bands
  match the standard AlphaFold visualization. The `next()` method on `ColorSchemeType`
  conditionally includes pLDDT only when detected, which is a nice UX touch.
- **Tests**: Includes unit test verifying all four bands return correct colors, and
  tests for the `has_plddt()` heuristic.
- **Do we have equivalent?**: No. Worth adopting.

### 2.3 Camera Orientation Matrix (Mirroring Bug Fix)

- **Description**: Replaces the three Euler angle fields (`rot_x`, `rot_y`, `rot_z`)
  with a 3x3 orientation matrix and a pivot point. Rotations are applied as local
  axis rotations via Rodrigues' rotation formula composed with the existing matrix.
  The `project()` method flips screen-space X (`-x_view`) to fix a mirroring issue.
  Adds `rotate_vector()` for normal transformation without zoom/pan.
- **Files**: `src/render/camera.rs` (complete rewrite, 88 -> 188 lines)
- **Quality**: Clean and mathematically correct. The `mat_mul`, `mat_vec_mul`, and
  `rotation_matrix` functions are standard. The pivot support enables the interface
  center-of-mass and chain-focus features. Tests verify composition order and pivot
  behavior.
- **Impact**: This is a **breaking change** that affects every caller of
  `camera.rot_x`/`rot_y`/`rot_z` (they no longer exist) and `rotate_normal()` in
  `hd.rs`. The mirroring fix (flipping X in projection) is genuinely important --
  our current code renders structures mirrored left-to-right. The Euler angle approach
  also suffers from gimbal lock, though in practice that is rarely encountered.
- **Do we have equivalent?**: No. The mirroring bug exists in our codebase. The
  orientation matrix is a strictly better representation.

### 2.4 Beta Sheet Rendering Improvements

- **Description**: Adds `frame_hint` guidance to the spline frame computation so
  beta sheets are oriented along the real peptide plane direction (C=O carbonyl
  vector) rather than relying solely on parallel-transported normals. The
  `apply_sheet_frame_guides()` function collects nearby carbonyl direction hints,
  projects them perpendicular to the tangent, aligns signs for consistency, and
  blends 65/35 with the transported frame. Also adds `residue_frame_hint()` which
  extracts the C->O direction from backbone atoms.
- **Files**: `src/render/ribbon.rs` (+170 lines in `compute_frames`,
  `apply_sheet_frame_guides`, `project_perpendicular`, `align_vector_sign`,
  `residue_frame_hint`)
- **Quality**: Good. The approach of using carbonyl directions as orientation guides
  is physically motivated -- it is how real molecular graphics programs (PyMOL,
  ChimeraX) orient sheet ribbons. The sign-coherence logic prevents flipping between
  adjacent residues. The 65/35 blend between guide and transport provides smoothness.
- **Tests**: Test that carbonyl direction is preferred, test that flipped hints
  produce coherent binormals across a sheet run.
- **Do we have equivalent?**: No. Our sheets use pure parallel transport which can
  produce arbitrary rotation of the ribbon plane. This is a visual improvement.

### 2.5 Adaptive Spline Tessellation

- **Description**: Adds `generate_ribbon_mesh_adaptive()` which subdivides spline
  segments recursively based on projected screen-space error (chord length and
  midpoint deviation thresholds). This replaces the fixed `SPLINE_SUBDIVISIONS = 14`
  per segment with camera-aware subdivision (max depth 8). Also adapts coil tube
  cross-section segment count based on projected radius.
- **Files**: `src/render/ribbon.rs` (+100 lines: `subdivide_segment`,
  `adaptive_segment_needs_split`, `adaptive_coil_segments`,
  `generate_chain_ribbon_adaptive`)
- **Quality**: Reasonable. The screen-space error metric is sound. However, this
  introduces a dependency on `Camera` into the ribbon mesh generation, which means
  the mesh must be regenerated every frame (no caching). The fork compensates by
  removing the mesh cache entirely and regenerating per-frame. This is a performance
  tradeoff -- fewer triangles when zoomed out, but no reuse when only rotating.
- **Do we have equivalent?**: No. Our approach is fixed tessellation with mesh
  caching. The fork's approach produces better quality at the cost of per-frame
  mesh regeneration.

### 2.6 Tiled Parallel Rasterization

- **Description**: Replaces per-triangle `rasterize_triangle_depth()` with a tiled
  approach (`rasterize_triangles_tiled`). Triangles are binned into 64x64 pixel tiles,
  then tiles are rasterized in parallel using `std::thread::scope`. Each tile has its
  own local color/depth buffer, blitted back to the main framebuffer. Also parallelizes
  `apply_depth_tint()` and `to_rgb_image()` conversion.
- **Files**: `src/render/framebuffer.rs` (+250 lines)
- **Quality**: Solid engineering. The tile binning is correct. The `PreparedTriangle`
  struct pre-computes shading and bounding boxes. Thread count is auto-detected via
  `available_parallelism()`. Fallback to single-threaded for small workloads
  (<16384 pixels or 1 tile). The `blit_tile` method copies scanlines efficiently.
- **Concern**: Introduces `std::thread` dependency into the framebuffer module. The
  original `rasterize_triangle_depth` is preserved only under `#[cfg(test)]`.
- **Do we have equivalent?**: No. Our rasterizer is single-threaded. This would help
  with HD rendering performance on large structures.

### 2.7 Cartoon Edge Outlines

- **Description**: Adds a post-processing pass in HD cartoon mode that detects edge
  pixels (background neighbors or depth discontinuities) and darkens them by 65% to
  create a cel-shading/outline effect.
- **Files**: `src/render/hd.rs` (`apply_cartoon_outline`, +40 lines)
- **Quality**: Simple and effective. The depth threshold (0.14) and darken factor
  (0.35) are hardcoded but reasonable. Only applied in HD cartoon mode.
- **Do we have equivalent?**: No. Visually appealing addition.

### 2.8 Chain Focus Mode (`/` keybinding)

- **Description**: Pressing `/` advances to the next chain, colors the focused chain
  normally while graying out others, sets the camera pivot to the chain's center of
  mass, and fits the zoom to the chain's bounding radius. Adds `Focus` color scheme
  variant and `new_focus_chain()` constructor. Adds `chain_center_of_mass()`,
  `chain_bounding_radius_from_pivot()`, and `focus_next_chain()` to `App`.
- **Files**: `src/app.rs` (+80 lines), `src/render/color.rs` (+20 lines),
  `src/main.rs` (+10 lines)
- **Quality**: Clean. The center-of-mass calculation uses atomic masses (C=12.011,
  N=14.007, etc.) which is more accurate than a simple geometric centroid. The zoom
  calculation reuses the auto-zoom logic cleanly.
- **Do we have equivalent?**: No. Useful for multi-chain structures.

### 2.9 HD Mode Toggle Preserves View

- **Description**: When toggling between braille and HD mode (`m`), the camera zoom
  and pan are scaled by the ratio of the old/new auto-zoom values so the structure
  stays at approximately the same visual size and position.
- **Files**: `src/app.rs` (`toggle_hd_mode_preserve_view`, +15 lines)
- **Quality**: Clean. Simple ratio-based rescaling.
- **Do we have equivalent?**: No. Currently toggling HD mode resets zoom.

### 2.10 Interface Pivot Point

- **Description**: When interface mode is active, the camera pivot is set to the
  center of mass of interface residues on the current focus chain, so rotation orbits
  around the binding site.
- **Files**: `src/app.rs` (`interface_pivot_for_chain`, `refresh_interface_pivot`,
  +30 lines)
- **Quality**: Clean. Integrates well with the camera pivot system.
- **Do we have equivalent?**: No. Requires the orientation-matrix camera.

### 2.11 Ligand Rendering (Fork's Approach)

- **Description**: Renders non-amino-acid, non-water, non-nucleotide residues as cyan
  sticks. Adds `is_ligand_residue()`, `is_amino_acid()`, `is_water()` predicates.
  Ligands stay in the chain's residue list and are filtered at render time. In
  wireframe mode, ligand residues are skipped from normal rendering and rendered
  separately.
- **Files**: `src/model/protein.rs` (+30 lines), `src/render/hd.rs` (+55 lines),
  `src/render/braille.rs` (+40 lines)
- **Quality**: Functional but inferior to our approach. All ligands are colored cyan
  regardless of element type. No ion detection or differentiated rendering. No
  element-based coloring of ligand atoms. Constant pixel radii rather than
  element-scaled sizes. The `is_ligand_residue` check runs on every residue during
  rendering, which is less clean than separating at parse time.
- **Do we have equivalent?**: Yes, and ours is better. We have separate `Ligand`
  structs with `LigandType::Ligand` vs `LigandType::Ion`, element-based CPK coloring
  for ligand atoms, element-scaled atom radii, binding pocket analysis, and a toggle
  (`show_ligands`). **Reject the fork's approach.**

### 2.12 Braille Renderer Removal

- **Description**: The fork effectively abandons the ratatui Canvas-based braille
  renderer. In non-HD mode, the viewport now calls `render_hd_framebuffer()` with
  `is_hd_output: false` and converts the result to braille via
  `framebuffer_to_braille_widget()`. The Canvas-based `render_protein()` in
  `braille.rs` is retained but no longer called from the viewport.
- **Files**: `src/ui/viewport.rs` (-10, +12 lines), `src/render/mod.rs` (removes
  `pub mod braille`)
- **Quality**: Debatable. This unifies the rendering path (cartoon mode works in
  braille now, which it didn't before in the fork), but it means braille mode does
  full software rasterization and then downsamples to 2x4 braille cells. The original
  Canvas approach was lighter-weight for backbone/wireframe modes.
- **Do we have equivalent?**: We still use the Canvas braille renderer. Keeping both
  paths is reasonable.

### 2.13 CI/CD Workflow + PyPI Packaging

- **Description**: Adds a GitHub Actions workflow that checks if the current
  Cargo.toml version is already on PyPI, builds Linux x86_64 wheels using maturin,
  and publishes to PyPI. Adds `pyproject.toml` for maturin integration.
- **Files**: `.github/workflows/workflow.yml` (90 lines), `pyproject.toml` (23 lines)
- **Quality**: The workflow is functional but x86-only. Uses trusted publisher (OIDC
  token) for PyPI, which is the modern best practice. The `pyproject.toml` uses
  `binaries = ["proteinview"]` to package the Rust binary. However, this is a
  distribution/packaging concern that we should set up ourselves rather than
  cherry-picking.
- **Do we have equivalent?**: No CI/CD. Worth creating our own, but not by
  cherry-picking theirs (our version is already 0.2.6, theirs is 0.1.4).

### 2.14 CLI Color/Mode Arguments

- **Description**: Adds `from_cli()` methods to `ColorSchemeType` and `VizMode` so
  the `--color` and `--mode` CLI arguments are properly parsed. Changes default
  `--mode` from `backbone` to `cartoon`.
- **Files**: `src/render/color.rs` (+15 lines), `src/app.rs` (+8 lines),
  `src/main.rs` (+5 lines)
- **Quality**: Clean. A small but useful UX improvement. The default mode change to
  cartoon makes sense since cartoon is the most visually informative.
- **Do we have equivalent?**: Partially. We have the CLI args but the default mode
  change and `from_cli` parsing are missing in our codebase.

### 2.15 Rotation Direction Fix

- **Description**: Swaps `j`/`k` rotation direction (`j` now rotates X by -1.0
  instead of +1.0, `k` by +1.0 instead of -1.0) to match "natural scroll" behavior
  where `j` (down) tilts the top of the molecule toward you.
- **Files**: `src/main.rs` (2 lines changed)
- **Quality**: Subjective. The fork also removed the rotation switching toggle that
  was briefly added, keeping only the new direction.
- **Do we have equivalent?**: No. Whether this is an improvement depends on user
  preference. The conventional direction (matching vim) might be more intuitive.

### 2.16 Heteroatom Coloring in Interface Wireframe

- **Description**: When in interface color mode and wireframe visualization, interface
  residues show N (blue) and O (red) atom colors instead of the uniform interface
  green/orange.
- **Files**: `src/render/color.rs` (`atom_color` match arm, +8 lines)
- **Quality**: Nice UX detail.
- **Do we have equivalent?**: No.

---

## 3. Bug Fixes

### 3.1 Mirroring Bug (CRITICAL)

- **Commit**: `a829aa5`
- **Description**: The projection was rendering structures mirrored left-to-right.
  Fixed by negating the X component in the projection output (`x: -x_view * self.zoom`).
- **Impact**: This is a real bug in our codebase too. Structures appear as mirror
  images of their true chirality. For a protein viewer this matters -- L-amino acids
  will look like D-amino acids.
- **Recommendation**: **Must fix.** Can be done independently of the camera rewrite
  by adding `x: -x3 * self.zoom` in our `project()` method.

### 3.2 Chain Cycling Fix

- **Commit**: `c7ca26a`
- **Description**: Fixes chain cycling to also update interface pivot when interface
  mode is active.
- **Impact**: Minor. Only affects interface mode UX.
- **Recommendation**: Adopt when integrating pivot support.

---

## 4. Architecture Changes

### 4.1 Camera: Euler Angles -> Orientation Matrix

The most impactful architectural change. Our camera stores three f64 angles; the
fork stores a 3x3 rotation matrix. This requires changing:
- `camera.rs`: Complete rewrite
- `hd.rs`: `rotate_normal()` becomes `camera.rotate_vector()`
- `main.rs`: Remove direct access to `rot_x`/`rot_y`/`rot_z`
- All tests referencing camera rotation angles

The orientation matrix is strictly superior (no gimbal lock, correct composition
order, enables pivot support), but it is a large change touching many files.

### 4.2 Mesh Cache Removal

The fork removes `mesh_cache` and `mesh_dirty` from `App` and regenerates the ribbon
mesh every frame. This is a consequence of the adaptive tessellation (which needs
camera info). Performance-wise, this trades mesh generation cost for fewer triangles
when zoomed out. For large structures in HD mode, the tiled parallel rasterizer
compensates.

Our current approach (cache mesh, regenerate only on color change) is more efficient
for rotation-only interactions. A hybrid approach (cache at a given LOD, invalidate
on significant zoom change) would be ideal but neither codebase implements that.

### 4.3 Ligand Data Model

The fork keeps ligands as chain residues and filters them in rendering code.
Our approach separates them into `Protein.ligands: Vec<Ligand>` at parse time with
a `LigandType` enum. Our approach is cleaner and supports richer functionality
(binding pocket analysis, element-based coloring, ion vs ligand distinction).

### 4.4 Braille Rendering Path

The fork routes non-HD rendering through the same framebuffer pipeline as HD,
then converts to braille characters. This eliminates the Canvas-based braille
renderer and unifies the code path. The downside is doing full rasterization even
for the lower-resolution braille output.

### 4.5 Protein Model Extensions

The fork adds `has_plddt()`, `is_ligand_residue()`, `is_amino_acid()`, `is_water()`
to the protein model. It also removes `backbone_atoms()` (replacing call sites with
inline iteration). We should adopt `has_plddt()` but keep our existing ligand model.

---

## 5. Conflict Analysis

### 5.1 High Conflict Files

| File | Conflict Reason |
|------|----------------|
| `src/model/protein.rs` | Fork has no `Ligand` struct, no `ligands` vec, no `is_hetero` field. We have `LigandType`, `WATER_NAMES`, `COMMON_IONS`, `Ligand` struct, binding-pocket fields. Completely different data models. |
| `src/app.rs` | Fork has no `show_ligands`, no `analyze_binding_pockets`. We have no `focus_next_chain`, no `toggle_hd_mode_preserve_view`, no `chain_center_of_mass`. Signatures differ (`App::new` has different args). |
| `src/render/hd.rs` | Fork has no `render_ligands_fb` (uses `render_ligand_sticks_fb` inline). We have separate `render_ligands_fb` with `LigandType` dispatch. Different `render_hd_framebuffer` signatures. |
| `src/render/braille.rs` | Fork strips out ligand rendering entirely (uses framebuffer path). We have `render_ligands` with `LigandType` support. |
| `src/render/camera.rs` | Complete rewrite. No shared API surface. |
| `src/render/color.rs` | Fork changes `next()` signature to take `has_plddt: bool`. We have the same `next()` but no pLDDT. Fork has no `ligand_atom_color()`. |
| `src/parser/pdb.rs` | Fork calls `infer_protein_secondary_structure()`. We don't. Otherwise similar. |

### 5.2 Low Conflict Files

| File | Notes |
|------|-------|
| `src/model/secondary.rs` | Fork adds inference functions. Our version is identical for existing code (just formatting diffs). Clean additive change. |
| `src/render/ribbon.rs` | Fork adds adaptive tessellation and frame guides. Existing code is structurally the same. Additive with minor refactoring. |
| `src/render/framebuffer.rs` | Fork adds tiled rasterizer. Existing code preserved under `#[cfg(test)]`. Mostly additive. |
| `src/ui/*` | Mostly formatting changes. Help overlay/helpbar add new keybinding entries. Low conflict. |

### 5.3 Features We Have That the Fork Lacks

- RNA/DNA support (MoleculeType::RNA/DNA, nucleic acid ribbon, backbone C4')
- Proper small molecule data model (Ligand struct, LigandType enum)
- Binding pocket analysis
- Ion-specific rendering (spheres vs ball-and-stick)
- Element-based ligand coloring (CPK colors)
- `show_ligands` toggle
- `is_hetero` field on Atom
- NMR multi-model handling (loading only first MODEL)
- SmallMolecule variant in MoleculeType

---

## 6. Cherry-Pick Recommendations

### 6.1 Adopt (High Value)

| Feature | Approach |
|---------|----------|
| **SS inference from coordinates** | Cherry-pick the inference functions from `secondary.rs` and the `infer_protein_secondary_structure()` call in `pdb.rs`. Clean additive change with minimal conflicts. |
| **pLDDT coloring** | Add `Plddt` variant to `ColorSchemeType`, add `plddt_color()` to `ColorScheme`, add `has_plddt()` to `Protein`. Adapt `next()` to conditionally include it. Adapt CLI parsing. |
| **Mirroring bug fix** | Fix in our `camera.rs` by negating X in projection output. One-line fix that doesn't require the full camera rewrite. |
| **Beta sheet frame guides** | Port `residue_frame_hint()`, `frame_hint` field on `SplinePoint` and `CaRecord`, `apply_sheet_frame_guides()`, and `compute_frames()` refactor. Moderate effort but significant visual improvement. |

### 6.2 Adopt (Medium Value, Higher Effort)

| Feature | Approach |
|---------|----------|
| **Camera orientation matrix** | Full rewrite of camera.rs. Better architecture but touches many files. Do this as a dedicated PR. |
| **Tiled parallel rasterization** | Port `rasterize_triangles_tiled()` and the parallel `apply_depth_tint`. Clean additive to framebuffer.rs. Keep single-triangle method for non-batched use. |
| **Cartoon edge outlines** | Port `apply_cartoon_outline()` to hd.rs. Small, self-contained. |
| **HD mode toggle preserves view** | Port `toggle_hd_mode_preserve_view()`. Requires refactoring `recalculate_zoom` into `auto_zoom_for_radius`. |

### 6.3 Consider (Low Priority)

| Feature | Notes |
|---------|-------|
| **Chain focus mode** | Nice UX but requires camera pivot support. Adopt after camera rewrite. |
| **Adaptive tessellation** | Good quality but kills mesh caching. Consider a hybrid approach instead. |
| **Interface pivot** | Requires camera pivot. Adopt after camera rewrite. |
| **CLI from_cli methods** | Minor convenience. Easy to add. |
| **Rotation direction swap** | Subjective. Test both and decide. |

### 6.4 Reject

| Feature | Reason |
|---------|--------|
| **Fork's ligand rendering** | Our `Ligand`/`LigandType` approach is strictly better. Element-based coloring, ion detection, size variation, binding pocket analysis. |
| **Braille renderer removal** | Keep our Canvas braille renderer. It is lighter-weight for non-HD mode and supports our ligand rendering. |
| **Fork's version numbers** | They use 0.1.x. We are at 0.2.6. Irrelevant. |
| **CI/CD workflow** | Create our own. The fork's is x86-only and references their PyPI package. |
| **pyproject.toml** | Create our own. Different version, different metadata. |
| **Example CIF files** | 1849+1824 lines of test data. Only useful if we want to test pLDDT on Boltz2 models. Consider adding an AlphaFold PDB test file instead. |
| **Removal of `backbone_atoms()`** | We still use it in braille.rs. Keep it. |

---

## 7. Bloat Assessment

### 7.1 Autoformatting Noise

Approximately 40-50% of the diff line count is `ruff`/`rustfmt` reformatting:
- Breaking long argument lists into multi-line format
- Reordering `use` statements alphabetically
- Expanding single-line `if` blocks to multi-line
- Moving `ratatui::Frame` import position

This inflates the apparent 6500-line diff. The actual feature code is roughly
3000-3500 lines of substantive changes.

### 7.2 Example Files

The two CIF files (`OBP5_model_0.cif` at 1824 lines, `1-heptanol_model_0.cif` at
1849 lines) add 3673 lines to the repository. These are Boltz2 predictions used for
testing pLDDT coloring. Only one is needed for tests; both is overkill.

### 7.3 Dead Code

- `generate_ribbon_mesh()` (non-adaptive) is marked `#[allow(dead_code)]` but
  retained. The fork only calls `generate_ribbon_mesh_adaptive()`.
- `generate_chain_ribbon()` similarly marked `#[allow(dead_code)]`.
- The braille module's `render_protein()` is no longer called but remains exported.

### 7.4 Unnecessary Complexity

- The `atomic_mass()` function in `App` for center-of-mass calculation is overkill.
  A geometric centroid (equal weights) would give essentially the same result for
  protein chains. Mass-weighting matters only for small molecules with heavy atoms.
- The `estimate_amide_h()` function uses a bisector of C->N and CA->N directions.
  This is a simplification -- the actual amide H position depends on hybridization
  and nearby atoms. However, it is good enough for DSSP-style classification and
  matches the approach used in other lightweight tools.

### 7.5 Scope Creep

- The rotation direction toggle was added and then removed in subsequent commits
  (commits `061e103` then `5d49b91`). Development churn that should have been
  squashed.
- Version bumps appear four times in the commit history (`1532386`, `5d96d11`,
  `128945e`, `c4ebcf1`). Unnecessary noise.

---

## Summary: Priority Implementation Order

1. **Mirroring bug fix** -- one-line fix, critical correctness issue
2. **SS inference from coordinates** -- high value, clean additive change
3. **pLDDT coloring** -- high value, moderate effort, clean integration
4. **Beta sheet frame guides** -- moderate effort, significant visual improvement
5. **Camera orientation matrix** -- large effort, prerequisite for pivot/focus features
6. **Tiled parallel rasterization** -- moderate effort, performance improvement
7. **Cartoon edge outlines** -- small effort, visual polish
8. **HD mode toggle preserve view** -- small effort, UX improvement

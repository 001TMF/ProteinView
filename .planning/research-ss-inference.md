# Research: Secondary Structure Inference from Coordinates

This document analyzes the secondary structure inference implementation in the
`desperadus-master` fork and plans how to integrate equivalent functionality
into our codebase, with improvements.

---

## 1. Algorithm Deep-Dive: Fork's Implementation

The fork adds ~350 lines to `src/model/secondary.rs` implementing a simplified
DSSP-style secondary structure inference. The algorithm has five stages:

### 1.1 Torsion Angle Computation (`compute_torsions`)

For each residue `i` (except the first and last), computes the phi/psi backbone
dihedral angles:

- **phi**: dihedral(C[i-1], N[i], CA[i], C[i])
- **psi**: dihedral(N[i], CA[i], C[i], N[i+1])

Requires atoms `C`, `N`, `CA` on the current residue, `C` on the previous, and
`N` on the next. If any atom is missing, the residue gets `None` torsions.

The `dihedral()` function uses the standard four-atom dihedral angle calculation
via cross products of bond vectors and `atan2`, returning degrees in [-180, 180].

**Assessment**: Correct and standard. The implementation handles degenerate cases
(zero-length bond vectors) by returning 0.0, which is safe enough for
classification purposes.

### 1.2 Hydrogen Bond Map (`compute_hbond_map`)

Builds an `n x n` boolean matrix `hbonds[acceptor][donor]` indicating backbone
N-H...O=C hydrogen bonds. This is the core of the DSSP method.

#### 1.2.1 Amide Hydrogen Estimation (`estimate_amide_h`)

PDB files typically lack hydrogen coordinates. The fork estimates the amide H
position using the bisector method:

```
dir_prev = normalize(N[i] - C[i-1])
dir_ca   = normalize(N[i] - CA[i])
bisector = normalize(dir_prev + dir_ca)
H_pos    = N[i] + bisector * 1.0 Angstrom
```

This places H along the bisector of the two bond directions away from N, at a
fixed 1.0 A N-H bond length.

**Assessment**: This is a reasonable simplification. The actual DSSP program
(Kabsch & Sander, 1983) uses the same basic idea but with the C=O direction of
the *previous* residue rather than a bisector. Specifically, DSSP estimates H as:

```
H = N + normalize(N - O[i-1])   (using O of previous residue, not C)
```

The fork's bisector approach using C[i-1] and CA[i] is a variant that produces
similar results in most cases. The 1.0 A bond length is standard. The main
limitation is that it cannot estimate H for the first residue in a chain (no
previous C atom), which is fine since the first residue cannot be an H-bond
donor anyway.

**Potential issue**: The function takes `_n_atoms` and `_ca_atoms` as unused
parameters (prefixed with underscore). This suggests the author considered
alternative estimation approaches but settled on the bisector. The unused
parameters are dead code but harmless.

#### 1.2.2 Hydrogen Bond Energy (`hbond_energy`)

Uses the DSSP electrostatic energy formula:

```
E = 27.888 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) kcal/mol
```

Where:
- `r_ON` = distance(O, N)
- `r_CH` = distance(C, H)
- `r_OH` = distance(O, H)
- `r_CN` = distance(C, N)

All distances are clamped to a minimum of 0.5 A to avoid division by zero.

A bond is accepted when `E < -0.5 kcal/mol`.

**Assessment**: This exactly matches the DSSP paper (Kabsch & Sander, 1983).
The constant 27.888 = q1 * q2 * 332 where q1 = 0.42e (N-H partial charge),
q2 = 0.20e (C=O partial charge), and 332 is the conversion factor from e^2/A
to kcal/mol. The -0.5 threshold is the standard DSSP cutoff. The 0.5 A minimum
distance prevents numerical explosion but is extremely conservative (no real
interatomic distance should be that small in a reasonable structure).

**Good**: Formula is correct and matches published DSSP.

**Potential improvement**: The minimum distance clamping at 0.5 A is fine but we
could add a maximum distance pre-filter (e.g., skip pairs where N...O > 5.2 A)
to avoid computing the full energy for obviously non-bonded pairs. This is where
spatial indexing would help (see Section 4).

#### 1.2.3 H-Bond Search: O(n^2) All-Pairs

The fork iterates all `(acceptor, donor)` pairs with `acceptor != donor` and
`|acceptor - donor| > 1` (adjacent residues excluded). For each pair, it
extracts atom positions, estimates the amide H, and computes the energy.

**Assessment**: This is the major performance concern. For a protein with `n`
residues, this loop is `O(n^2)`. For a typical AlphaFold model (~300-600
residues), this means 90K-360K energy evaluations. For very large structures
(e.g., ribosomes, ~5000+ residues), this becomes 25M+ evaluations.

The atom position lookups are pre-extracted into Vecs at the top of the function
(`c_atoms`, `o_atoms`, `n_atoms`, `ca_atoms`), which avoids repeated linear
scans through residue atom lists. This is a good optimization.

### 1.3 Helix Assignment (`assign_helices_from_hbonds`)

Identifies helices using the DSSP hydrogen bond pattern:

- For turn sizes 4, 3, and 5 (alpha-helix, 3_10-helix, pi-helix respectively):
  check if `hbonds[i][i+turn]` is true (C=O of residue i bonds to N-H of
  residue i+turn).
- If the H-bond exists, check that at least half the residues in the spanned
  region (i+1 through i+turn) have torsion angles compatible with helix geometry.
- Accumulate a "support" count per residue. Any residue with support > 0 is
  assigned Helix.

**Assessment**: This is a reasonable simplification of DSSP. True DSSP identifies
"turns" (single i->i+n H-bonds) and then requires two consecutive overlapping
turns to form a helix. The fork's approach is more permissive -- a single H-bond
with torsion-angle support is enough. This will occasionally over-assign helices,
but the minimum-run-length filter (Section 1.5) compensates.

The torsion angle filtering is a smart addition that pure DSSP lacks. It prevents
spurious H-bond matches in irregular regions from being classified as helices.

### 1.4 Sheet Assignment (`assign_sheets_from_hbonds`)

Identifies beta sheets using parallel and antiparallel H-bond ladder patterns:

- **Antiparallel**: `(hbonds[i][j] && hbonds[j][i])` or
  `(hbonds[i-1][j+1] && hbonds[j-1][i+1])`
- **Parallel**: `(hbonds[i-1][j] && hbonds[j][i+1])` or
  `(hbonds[j-1][i] && hbonds[i][j+1])`

Only pairs where both `i` and `j` have sheet-compatible torsion angles are
considered. Both residues in a detected pair get support, and only residues still
assigned as Coil (not already Helix) can become Sheet.

**Assessment**: The H-bond ladder patterns match DSSP's bridge detection
(parallel and antiparallel). True DSSP identifies "bridges" (isolated paired
H-bonds) vs "ladders" (consecutive bridges that form sheets), then assigns
extended strand (E) to residues in ladders and bridge (B) to isolated bridges.
The fork simplifies this by treating any paired H-bond match as sheet-worthy,
relying on torsion angles and minimum-run-length to filter noise. This is
adequate for visualization purposes.

The inner loop iterates `j` from `i+2` to `n-1`, which is O(n^2) per chain.
Combined with the H-bond map construction, total complexity is O(n^2).

### 1.5 Post-Processing

#### Torsion Angle Windows

Four categories of torsion angle compatibility:

| Type | phi range | psi range |
|------|-----------|-----------|
| Strong helix | [-140, -80] | [-170, -100] |
| Weak helix | [-170, -40] | [-180, -60] |
| Strong sheet | [-100, -40] | [20, 90] **or** [80, 180] | [120, 180] |
| Weak sheet | [-140, -20] | [0, 180] **or** [60, 180] | [90, 180] |

The `torsions_match_target()` function accepts either strong or weak for the
target type. The distinction between strong/weak is only used in naming; both
pass the match test.

**Assessment**: These windows are reasonable but somewhat wider than the
Ramachandran plot core regions:
- Alpha helix core: phi ~ -57, psi ~ -47; the strong window [-140,-80] x
  [-170,-100] is a generous box that will also capture 3_10-helix and some turns.
- Beta sheet core: phi ~ -120, psi ~ 120; the strong window captures this well.
- The weak windows are very permissive. The weak helix window
  [-170,-40] x [-180,-60] covers most of the left half of the Ramachandran plot.
  This permissiveness is acceptable because the torsion check is combined with
  H-bond evidence -- it is not used alone.

**Concern**: The strong sheet window has a second case `[80,180] x [120,180]`
which captures left-handed conformations. This is unusual for standard beta
sheets but might catch a few edge cases. Could lead to false positives for
polyproline II helices in this region.

#### Gap Filling (`fill_single_residue_gaps`)

If residues i-1 and i+1 are both assigned the same SS type (Helix or Sheet) but
residue i is Coil, and residue i has compatible torsion angles (or missing
torsions), then residue i is promoted to match its neighbors.

**Assessment**: Good heuristic. Single-residue gaps in SS runs are usually
artifacts of marginal H-bond energies or missing backbone atoms, not real
structural breaks. Gap-filling is important for clean visualization because
the ribbon renderer looks bad with single-residue SS interruptions.

The `torsions[i].is_none()` case (missing backbone atoms) is also filled, which
is reasonable -- if the flanking residues are both helical, the middle one is
almost certainly helical too even if we cannot compute its torsions.

#### Minimum Run-Length Filtering (`retain_runs`)

Removes SS assignments from runs shorter than:
- Helix: minimum 3 residues
- Sheet: minimum 2 residues

**Assessment**: Appropriate minimums. Real alpha helices need at least ~4
residues for one full turn; 3 is generous but acceptable (captures 3_10 helices
and helix caps). Real beta strands can be as short as 2 residues if they are
part of a sheet with another strand providing the hydrogen bonds.

### 1.6 Summary: Fork Strengths and Weaknesses

**Strengths:**
- Correct DSSP energy formula
- Combined H-bond + torsion angle approach reduces false positives
- Skip logic for non-protein chains and chains with existing SS
- Gap-filling and run-length filtering produce clean assignments
- Pre-extraction of atom positions avoids repeated lookups
- Handles 3_10, alpha, and pi helices via the turn=[3,4,5] loop

**Weaknesses:**
- O(n^2) H-bond search with no distance pre-filtering or spatial indexing
- O(n^2) memory for the hbond boolean matrix (n * n bytes)
- Amide H estimation differs slightly from canonical DSSP
- Torsion windows are permissive (could over-assign in edge cases)
- No proline handling (proline has no amide H; should skip as donor)
- No beta-bridge vs beta-ladder distinction
- Sheet detection ignores the "register" of H-bond ladders
- Unused function parameters in `estimate_amide_h`

---

## 2. Integration Points in Our Codebase

### 2.1 Where to Insert the Inference Call

The inference should be called in `src/parser/pdb.rs` in the `load_structure()`
function, after the existing PDB HELIX/SHEET and CIF SS assignment attempts, and
before returning the `Protein`. Specifically, after line 97 in our current code:

```rust
// Current code in src/parser/pdb.rs, lines 82-99:
let mut protein = Protein { name, chains, ligands };

assign_from_pdb_file(&mut protein, path);  // Try PDB HELIX/SHEET

let all_coil = protein.chains.iter()
    .flat_map(|c| &c.residues)
    .all(|r| r.secondary_structure == SecondaryStructure::Coil);
if all_coil {
    let lower = path.to_lowercase();
    if lower.ends_with(".cif") || lower.ends_with(".mmcif") {
        assign_from_cif_file(&mut protein, path);
    }
}

// INSERT HERE: infer_protein_secondary_structure(&mut protein);

Ok(protein)
```

The fork does exactly this -- calls `infer_protein_secondary_structure(&mut protein)`
as the last step before `Ok(protein)`.

### 2.2 What the Inference Needs (Input)

Per residue, the inference needs these backbone atom positions:
- `N` (nitrogen) -- for H-bond donor, torsion angles
- `CA` (C-alpha) -- for torsion angles, amide H estimation
- `C` (carbonyl carbon) -- for H-bond acceptor, torsion angles, amide H estimation
- `O` (carbonyl oxygen) -- for H-bond acceptor

These are accessed via `atom_pos(residue, "N")` etc., which does a linear scan
of `residue.atoms` looking for `atom.name == name`.

The inference also needs:
- `chain.molecule_type` -- to skip non-Protein chains
- `residue.secondary_structure` -- to check if SS is already assigned

### 2.3 What the Inference Produces (Output)

The inference sets `residue.secondary_structure` to one of:
- `SecondaryStructure::Helix`
- `SecondaryStructure::Sheet`
- `SecondaryStructure::Coil` (default, unchanged)

It does NOT use `SecondaryStructure::Turn`.

### 2.4 How SS Types Drive Ribbon Rendering

In `src/render/ribbon.rs`, the `cross_section()` function dispatches on
`SecondaryStructure`:

- `Helix` -> flat ribbon cross-section (elliptical, 1.30 x 0.40 A)
- `Sheet` -> flat ribbon with optional arrowhead widening (1.50 x 0.20 A)
- `Turn | Coil` -> circular tube cross-section (0.40 A radius)

The spline interpolation assigns the nearer residue's SS type to each
subdivision point (`if t < 0.5 { cas[i1].ss } else { cas[i2].ss }`).

Sheet arrowheads are detected by scanning for contiguous sheet runs and marking
the last 2 residues as arrow regions.

**Key insight**: The ribbon renderer only needs the `SecondaryStructure` enum
value per residue. It does not care how the assignment was made (PDB records,
CIF records, or inference). This means the inference output plugs in
transparently.

### 2.5 Required Import Changes

In `src/parser/pdb.rs`, add to the import:
```rust
use crate::model::secondary::{assign_from_pdb_file, assign_from_cif_file, infer_protein_secondary_structure};
```

In `src/model/secondary.rs`, the new function must be `pub`:
```rust
pub fn infer_protein_secondary_structure(protein: &mut Protein) { ... }
```

---

## 3. Pitfall Analysis

### 3.1 Performance: O(n^2) H-Bond Search

**Problem**: The H-bond map construction iterates all `(acceptor, donor)` pairs.
For n residues this is O(n^2) time and O(n^2) memory (the `Vec<Vec<bool>>` matrix).

| n (residues) | Pairs evaluated | Memory (bool matrix) |
|---|---|---|
| 100 | ~10K | ~10 KB |
| 500 | ~250K | ~250 KB |
| 1000 | ~1M | ~1 MB |
| 5000 | ~25M | ~25 MB |
| 10000 | ~100M | ~100 MB |

For typical AlphaFold models (200-1000 residues), this is fine. For very large
structures (full ribosomes, viral capsids), it becomes problematic. A 10K-residue
structure would require ~100M energy evaluations and 100 MB for the boolean matrix.

**Mitigation**: Most H-bonds have N...O distances < 5.2 A. We can apply a distance
pre-filter to skip pairs that are obviously too far apart. Better yet, use spatial
indexing to find only nearby N-O pairs.

### 3.2 Performance: Memory for H-Bond Matrix

**Problem**: `Vec<Vec<bool>>` is an n x n matrix. For 5000 residues, this is 25
million booleans = 25 MB. Rust `bool` is 1 byte, so this is not compact.

**Mitigation options**:
1. Use `bitvec` or manual bit-packing to reduce to n^2/8 bytes
2. Do not store the matrix at all; instead process H-bonds on-the-fly, accumulating
   only the per-residue support counts needed for helix/sheet assignment
3. Use a sparse representation (e.g., `HashSet<(usize, usize)>`)

Option 2 is the most memory-efficient but requires restructuring the algorithm.
Option 3 is a good balance: typical proteins have O(n) H-bonds, not O(n^2).

### 3.3 Edge Case: Incomplete Backbone Atoms

**Problem**: Some residues may lack N, CA, C, or O atoms due to:
- Disordered regions with missing density (X-ray)
- Partial occupancy
- Truncated models
- Non-standard residues (modified amino acids)

**Current handling**: The fork returns `None` from `atom_pos()` and propagates
with early returns (`let else` / `continue`). Residues with missing backbone
atoms get `None` torsions and are skipped as H-bond donors/acceptors.

**Assessment**: This is correct. Missing backbone atoms mean we cannot compute
meaningful torsion angles or H-bond energies for those residues. They will
default to Coil, which is the safest visual representation. Gap-filling may
promote them if flanked by structured residues.

### 3.4 Edge Case: Non-Protein Chains (RNA/DNA)

**Problem**: Our codebase has RNA and DNA chains (via `MoleculeType`) that must
not be subjected to protein SS inference.

**Current handling**: The fork's `infer_protein_secondary_structure()` checks
`chain.molecule_type != MoleculeType::Protein` and skips non-protein chains.
This is exactly correct.

**Our advantage**: We already have `MoleculeType` with `RNA`, `DNA`, and
`SmallMolecule` variants. The fork only has `Protein`, `RNA`, `DNA` (no
`SmallMolecule`). Our skip logic is already covered.

### 3.5 Edge Case: Structures with Partial SS Records

**Problem**: Some PDB files have HELIX records but no SHEET records, or vice
versa. The fork's per-chain check (`has_existing_ss`) looks for ANY non-Coil
residue in the chain. If even one residue has a Helix assignment from PDB
records, the entire chain is skipped for inference.

**Assessment**: This is the correct behavior. PDB HELIX/SHEET records, when
present, are authoritative. If a file has HELIX records but no SHEET records,
the helix assignments are correct and the remaining Coil residues are genuinely
coil -- we should not second-guess them with inference.

However, there is a subtle case: AlphaFold3 CIF files sometimes have
`_struct_conf` sections that contain only a few HELX_P entries and miss many
others. In this case, per-chain checking is correct (if the chain has some
assignments, respect them), but per-protein checking would be wrong (other
chains may need inference).

The fork checks per-chain, which is correct.

### 3.6 Deviation from DSSP: Helix Detection

**DSSP algorithm**: A helix is identified when there are two consecutive
n-turns. An n-turn exists when residue i is H-bonded to residue i+n. Two
consecutive turns at positions i and i+1 define a helix from i+1 to i+n.

**Fork algorithm**: A single H-bond from i to i+n, combined with torsion angle
support for at least half the spanned residues, is enough. No requirement for
consecutive overlapping turns.

**Impact**: The fork will identify helices more aggressively than DSSP. A single
i->i+4 H-bond with compatible torsion angles will assign 3-4 residues as helix,
whereas DSSP would require a second overlapping H-bond (e.g., i+1 -> i+5) to
confirm. The minimum-run-length filter (3 residues for helix) partially
compensates but does not fully replicate DSSP's rigor.

For visualization purposes, this over-assignment is generally acceptable. It is
better to show a slightly longer helix than to miss one entirely.

### 3.7 Deviation from DSSP: Sheet Detection

**DSSP algorithm**: Identifies "bridges" (paired H-bonds between residues i and
j) and then groups consecutive bridges into "ladders." Residues participating
in ladders get extended strand (E) assignment.

**Fork algorithm**: Detects the same bridge H-bond patterns but does not
require bridges to be consecutive (no ladder grouping). Any single bridge pair
with compatible torsion angles counts.

**Impact**: This could produce isolated "sheet" residues that are bridges but
not parts of a real beta sheet. The minimum-run-length filter (2 residues for
sheet) helps but does not fully prevent this. In practice, isolated bridges are
rare in well-folded proteins.

### 3.8 Missing: Proline as H-Bond Donor

**Problem**: Proline lacks an amide hydrogen (the nitrogen is part of the
pyrrolidine ring). It cannot donate backbone H-bonds. The fork does not
explicitly skip proline as an H-bond donor.

**Impact**: The `estimate_amide_h` function will still estimate an H position
for proline using the bisector method, producing a fictitious H atom. This may
lead to false-positive H-bonds involving proline donors. However, the energy
threshold (-0.5 kcal/mol) and the geometric placement of the fictitious H
usually prevent this from being a real problem in practice.

**Recommendation**: Add a check for proline (residue name "PRO") and skip it
as a donor. This is a minor improvement.

---

## 4. Our Improvements Over the Fork

### 4.1 RNA/DNA Chain Skipping

We already have `MoleculeType::RNA`, `MoleculeType::DNA`, and
`MoleculeType::SmallMolecule`. The inference function's check
`chain.molecule_type != MoleculeType::Protein` will correctly skip all of these.
The fork lacks `SmallMolecule`, so adopting their code directly is safe (our
extra variant is handled by the `!=` check).

### 4.2 Ligands Stored Separately

Our ligands are in `protein.ligands: Vec<Ligand>`, completely separate from
chain residues. The inference iterates `chain.residues` only, so ligands
will never interfere. The fork stores ligands inside chain residues and needs
`is_ligand_residue()` checks in various places; we do not.

### 4.3 Spatial Indexing with rstar

We have `pdbtbx` with the `rstar` feature enabled, which brings in `rstar 0.12.2`
as a transitive dependency. We can use `rstar::RTree` to build a spatial index
of backbone nitrogen positions, then for each carbonyl (C, O) acceptor, query
for nearby N atoms within a distance cutoff (e.g., 5.2 A for N...O distance).

This would reduce the H-bond search from O(n^2) to O(n log n) average case:
- Build R-tree of N atom positions: O(n log n)
- For each of n acceptors, query nearby donors: O(log n) per query
- Total: O(n log n)

Memory for the R-tree is O(n) vs O(n^2) for the boolean matrix.

**Implementation sketch**:
```rust
use rstar::{RTree, PointDistance, AABB};

struct NAtomEntry {
    pos: [f64; 3],
    residue_idx: usize,
}

impl rstar::RTreeObject for NAtomEntry {
    type Envelope = AABB<[f64; 3]>;
    fn envelope(&self) -> Self::Envelope {
        AABB::from_point(self.pos)
    }
}

impl rstar::PointDistance for NAtomEntry {
    fn distance_2(&self, point: &[f64; 3]) -> f64 {
        let dx = self.pos[0] - point[0];
        let dy = self.pos[1] - point[1];
        let dz = self.pos[2] - point[2];
        dx*dx + dy*dy + dz*dz
    }
}
```

Then query with `rtree.locate_within_distance(o_pos, 5.2 * 5.2)` for each
acceptor O atom, only evaluating H-bond energy for nearby pairs.

**Trade-off**: Adding the spatial index increases code complexity. For typical
AlphaFold models (200-600 residues), the O(n^2) approach takes <1ms and the
optimization is unnecessary. We should consider whether to:

(a) Ship the simple O(n^2) version first and optimize later if profiling shows
    it matters, or
(b) Implement spatial indexing from the start for correctness at scale.

**Recommendation**: Ship O(n^2) first with a distance pre-filter
(`if distance(o, n_atom) > 5.2 { continue; }`). This simple check avoids the
expensive H estimation and energy calculation for distant pairs, reducing the
constant factor dramatically. Add R-tree only if profiling shows it is needed
for real user workloads. The distance pre-filter alone should make structures
up to ~2000 residues fast enough.

### 4.4 Eliminate the Boolean Matrix

Instead of building the full `n x n` hbond matrix and then scanning it for
patterns, we can process H-bonds on-the-fly:

1. Compute torsions (same as fork)
2. For each H-bond found, immediately check if it matches a helix pattern
   (i -> i+3/4/5) or could be a bridge pattern, and accumulate per-residue
   support counts directly.
3. For sheet bridges, use a `Vec<(usize, usize)>` sparse list of detected
   bridges instead of the full matrix.

This reduces memory from O(n^2) to O(n) + O(bridges) which is typically O(n).

### 4.5 Proline Handling

Add a check: if the donor residue name is "PRO", skip it as a donor. This is
a one-line addition in the H-bond search loop.

### 4.6 Unit Tests for H-Bond Energy Formula

The fork has no unit test for `hbond_energy()` itself. We should add tests with
known values:

- A canonical alpha-helix H-bond (known geometry) should produce E < -0.5
- A pair of residues with H...O distance > 4.0 A should produce E ~ 0
- The function should handle degenerate inputs (overlapping atoms) gracefully
  due to the 0.5 A minimum distance clamping

### 4.7 Explicit Distance Pre-Filter

Add an N...O distance check before computing the full H-bond energy. In
real protein structures, backbone H-bonds have N...O distances of 2.7-3.5 A.
Using a generous cutoff of 5.2 A (same as DSSP) eliminates ~95% of pairs
without false negatives.

---

## 5. Test Strategy

### 5.1 Test Files

| File | Purpose | Why |
|------|---------|-----|
| `examples/AF3_TNFa.pdb` | **Primary test**: AlphaFold model without HELIX/SHEET records | This is the exact use case motivating the feature. TNF-alpha is a beta-sandwich protein with some helices. Inference should find both. |
| `examples/1UBQ.pdb` | **Regression test**: PDB with explicit HELIX/SHEET records | Verify that inference does NOT overwrite existing SS assignments. Ubiquitin has well-characterized SS from X-ray. |
| `examples/4HHB.pdb` | **Multi-chain test**: Hemoglobin with 4 chains and ligands | Verify inference works per-chain, skips chains with existing SS, and ligands (HEM) do not interfere. |
| `examples/1RNA.pdb` | **Non-protein skip test**: RNA structure | Verify inference skips RNA chains entirely. |
| `examples/1BNA.pdb` | **Non-protein skip test**: DNA structure | Verify inference skips DNA chains entirely. |
| `examples/1AOI.pdb` | **Mixed content test**: Protein-DNA complex | Verify inference runs on protein chains but skips DNA chains in the same structure. |

### 5.2 Unit Tests

#### 5.2.1 `hbond_energy` Unit Test
```
Test: Known alpha-helix H-bond geometry
Input: Typical i -> i+4 H-bond with r_ON=2.9A, r_OH=1.9A, r_CH=3.5A, r_CN=3.8A
Expected: E < -0.5 kcal/mol

Test: Distant pair
Input: r_ON=8.0A, r_OH=7.5A, r_CH=9.0A, r_CN=8.5A
Expected: E close to 0

Test: Degenerate input (atoms very close)
Input: r_ON=0.1A (will be clamped to 0.5)
Expected: No panic, finite result
```

#### 5.2.2 `dihedral` Unit Test
```
Test: Known planar geometry
Input: Four atoms forming a 180-degree dihedral (trans conformation)
Expected: ~180 or ~-180 degrees

Test: Known geometry
Input: Four atoms forming a 0-degree dihedral (cis conformation)
Expected: ~0 degrees
```

#### 5.2.3 `estimate_amide_h` Unit Test
```
Test: Known geometry
Input: N, C[i-1], CA with known positions
Expected: H placed at ~1.0A from N along bisector direction
```

#### 5.2.4 Torsion Angle Windows
```
Test: Alpha helix core (phi=-57, psi=-47)
Expected: is_strong_helix_torsion returns true

Test: Beta sheet core (phi=-120, psi=120)
Expected: is_strong_sheet_torsion returns true

Test: Left-handed helix (phi=60, psi=40)
Expected: Both helix and sheet torsion checks return false
```

#### 5.2.5 `retain_runs` Unit Test
```
Test: Short helix run (2 residues) filtered out
Input: [Coil, Helix, Helix, Coil]
Expected: [Coil, Coil, Coil, Coil] (min_len=3 for helix)

Test: Long enough helix run kept
Input: [Coil, Helix, Helix, Helix, Coil]
Expected: [Coil, Helix, Helix, Helix, Coil]
```

#### 5.2.6 `fill_single_residue_gaps` Unit Test
```
Test: Single-residue gap between helices filled
Input: [Helix, Coil, Helix], torsions[1] = Some((-60, -45)) (helix-compatible)
Expected: [Helix, Helix, Helix]

Test: Single-residue gap with incompatible torsions NOT filled
Input: [Helix, Coil, Helix], torsions[1] = Some((60, 60)) (not helix)
Expected: [Helix, Coil, Helix]
```

### 5.3 Integration Tests

#### 5.3.1 AlphaFold Model Gets SS Assigned
```rust
#[test]
fn test_infer_ss_alphafold_model() {
    let protein = load_structure("examples/AF3_TNFa.pdb").unwrap();
    let helix_count = protein.chains.iter()
        .flat_map(|c| &c.residues)
        .filter(|r| r.secondary_structure == SecondaryStructure::Helix)
        .count();
    let sheet_count = protein.chains.iter()
        .flat_map(|c| &c.residues)
        .filter(|r| r.secondary_structure == SecondaryStructure::Sheet)
        .count();
    // TNF-alpha should have both helices and substantial beta-sheets
    assert!(helix_count > 0, "Expected inferred helices");
    assert!(sheet_count > 0, "Expected inferred sheets");
    // Total structured residues should be a significant fraction
    let total = protein.residue_count();
    let structured = helix_count + sheet_count;
    assert!(structured as f64 / total as f64 > 0.2,
        "Expected >20% structured residues, got {}/{}", structured, total);
}
```

#### 5.3.2 Existing PDB SS Not Overwritten
```rust
#[test]
fn test_existing_ss_preserved() {
    let protein = load_structure("examples/1UBQ.pdb").unwrap();
    let chain = &protein.chains[0];
    // 1UBQ residue 10 is in a beta sheet (SHEET record)
    let res10 = chain.residues.iter().find(|r| r.seq_num == 10).unwrap();
    assert_eq!(res10.secondary_structure, SecondaryStructure::Sheet);
    // 1UBQ residue 23 is in a helix (HELIX record)
    let res23 = chain.residues.iter().find(|r| r.seq_num == 23).unwrap();
    assert_eq!(res23.secondary_structure, SecondaryStructure::Helix);
}
```

#### 5.3.3 RNA/DNA Chains Skipped
```rust
#[test]
fn test_rna_chain_not_inferred() {
    let protein = load_structure("examples/1RNA.pdb").unwrap();
    // All residues should remain Coil (no inference applied to RNA)
    let all_coil = protein.chains.iter()
        .filter(|c| c.molecule_type != MoleculeType::Protein)
        .flat_map(|c| &c.residues)
        .all(|r| r.secondary_structure == SecondaryStructure::Coil);
    assert!(all_coil, "RNA residues should remain Coil");
}
```

#### 5.3.4 Ligands Do Not Interfere
```rust
#[test]
fn test_ligands_dont_interfere_with_inference() {
    // 4HHB has HEM ligands; verify they are in protein.ligands
    // and chain residues only contain amino acids
    let protein = load_structure("examples/4HHB.pdb").unwrap();
    assert!(protein.ligand_count() > 0);
    // Inference should still produce structured residues
    let structured = protein.chains.iter()
        .flat_map(|c| &c.residues)
        .filter(|r| r.secondary_structure != SecondaryStructure::Coil)
        .count();
    assert!(structured > 0);
}
```

### 5.4 Performance Benchmark (Manual)

For development validation, manually time the inference on a large structure.
Consider downloading a ribosome structure (~5000+ residues) to verify that
inference completes in a reasonable time (<5 seconds). This is not an automated
test but a manual checkpoint.

---

## 6. Implementation Plan Summary

### Phase 1: Direct Port (MVP)

1. Add the inference functions to `src/model/secondary.rs`:
   - `infer_protein_secondary_structure()` (public entry point)
   - `infer_chain_secondary_structure()` (per-chain logic)
   - `compute_torsions()`
   - `compute_hbond_map()` with distance pre-filter
   - `estimate_amide_h()`
   - `hbond_energy()`
   - `assign_helices_from_hbonds()`
   - `assign_sheets_from_hbonds()`
   - `fill_single_residue_gaps()`
   - `retain_runs()`
   - Torsion angle window functions
   - Vector math helpers (`sub`, `add`, `scale`, `dot`, `cross`, `norm`,
     `distance`, `normalize`, `dihedral`)

2. Add the call in `src/parser/pdb.rs` after the CIF fallback.

3. Add unit and integration tests.

### Phase 2: Improvements (Follow-up)

1. Add N...O distance pre-filter (5.2 A cutoff) in H-bond search
2. Add proline skip for H-bond donors
3. Replace boolean matrix with sparse H-bond list or on-the-fly processing
4. Consider R-tree spatial indexing if profiling shows need
5. Clean up unused parameters in `estimate_amide_h`
6. Add the `AF3_TNFa.pdb` example file (already present in our repo)

### Estimated Scope

- Lines of new code: ~350 (functions) + ~150 (tests) = ~500 lines
- Files modified: 2 (`src/model/secondary.rs`, `src/parser/pdb.rs`)
- Files created: 0
- Risk: Low (purely additive; no changes to existing functionality)
- Dependencies: None new (rstar available but not needed for Phase 1)

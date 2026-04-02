# Research: pLDDT Confidence Coloring

## 1. Algorithm

### Standard AlphaFold pLDDT 4-Band Color Scheme

pLDDT (predicted Local Distance Difference Test) is AlphaFold's per-residue confidence
metric, scored 0-100. The standard visualization uses four discrete color bands:

| pLDDT Range | Color       | RGB              | Meaning                    |
|-------------|-------------|------------------|----------------------------|
| >= 90       | Dark blue   | (0, 83, 214)     | Very high confidence       |
| 70 - 89     | Light blue  | (101, 203, 243)  | Confident                  |
| 50 - 69     | Yellow      | (255, 219, 19)   | Low confidence             |
| < 50        | Orange      | (255, 125, 69)   | Very low confidence        |

These exact RGB values are used in the `desperadus-master` fork's `plddt_color()` function
and match the official AlphaFold color scheme used by the EBI/AlphaFold DB website and
Mol* viewer.

The implementation averages the B-factor values across all atoms in a residue:
```
avg_plddt = sum(atom.b_factor for atom in residue.atoms) / len(residue.atoms)
```
This averaging is appropriate because pLDDT is a per-residue score -- all atoms within
the same residue should have identical B-factor values in a well-formed AlphaFold PDB.
The averaging provides robustness if values differ slightly.

### `has_plddt()` Detection Heuristic

The fork implements a heuristic on `Protein` to auto-detect whether B-factors are
actually pLDDT scores. It iterates all atoms across all chains and computes three
statistics:

```
in_range_fraction  = count(0 <= b_factor <= 100) / total_atoms
mean               = sum(b_factor) / total_atoms
high_conf_fraction = count(b_factor >= 70) / total_atoms
```

Detection triggers when ALL THREE conditions are met:

1. `in_range_fraction >= 0.95` -- nearly all values in [0, 100]
2. `mean >= 50.0` -- average is above midpoint
3. `high_conf_fraction >= 0.25` -- at least 25% of atoms have high confidence

**Rationale**: Crystallographic B-factors (Debye-Waller factors) measure atomic
displacement and are typically 5-80 A^2 for well-refined structures, but can exceed 100.
The key discriminator is the mean: experimental B-factors for well-resolved structures
average around 15-30 A^2, while pLDDT scores from AlphaFold typically average 70-90
since the model assigns high confidence to most of the structure.


## 2. Integration Points

### 2a. ColorSchemeType Enum (`src/render/color.rs`)

Add a `Plddt` variant to `ColorSchemeType`:

```rust
pub enum ColorSchemeType {
    Structure,
    Plddt,     // <-- new
    Chain,
    Element,
    BFactor,
    Rainbow,
    Interface,
}
```

Add a `plddt_color()` method to `ColorScheme`:

```rust
fn plddt_color(&self, residue: &Residue) -> Color {
    let avg_plddt = if residue.atoms.is_empty() {
        0.0
    } else {
        residue.atoms.iter().map(|a| a.b_factor).sum::<f64>() / residue.atoms.len() as f64
    };
    if avg_plddt >= 90.0 {
        Color::Rgb(0, 83, 214)
    } else if avg_plddt >= 70.0 {
        Color::Rgb(101, 203, 243)
    } else if avg_plddt >= 50.0 {
        Color::Rgb(255, 219, 19)
    } else {
        Color::Rgb(255, 125, 69)
    }
}
```

Wire it into `residue_color()`:
```rust
ColorSchemeType::Plddt => self.plddt_color(residue),
```

Add `name()` entry:
```rust
Self::Plddt => "pLDDT",
```

For ligand coloring under pLDDT mode, we should fall back to a neutral color or the
structure-mode ligand colors, since ligands do not have meaningful pLDDT scores.

### 2b. `next()` Cycling Method (`src/render/color.rs`)

The fork changes `next()` to accept `has_plddt: bool` and conditionally inserts pLDDT
into the cycle:

```rust
pub fn next(&self, has_plddt: bool) -> Self {
    match self {
        Self::Structure => if has_plddt { Self::Plddt } else { Self::Chain },
        Self::Plddt => Self::Chain,
        Self::Chain => Self::Element,
        // ... rest unchanged
    }
}
```

This means pLDDT only appears in the rotation when auto-detected (or forced via CLI).
The signature change propagates to `App::cycle_color()`.

### 2c. App State (`src/app.rs`)

`cycle_color()` must pass the detection result:

```rust
pub fn cycle_color(&mut self) {
    let next = self.color_scheme.scheme_type.next(self.protein.has_plddt());
    self.color_scheme = ColorScheme::new(next, self.protein.residue_count());
    self.mesh_dirty = true;
}
```

The fork also adds a `from_cli()` constructor on `ColorSchemeType` that maps the
`--color plddt` argument. It gates on `has_plddt` -- if the user requests pLDDT but the
file does not look like an AlphaFold model, it falls back to Structure.

### 2d. CLI Argument (`src/main.rs`)

The `--color` argument help text needs updating to include `plddt`:

```rust
/// Color scheme: structure, plddt, chain, element, bfactor, rainbow
#[arg(long, default_value = "structure")]
color: String,
```

Our current codebase does not have a `from_cli()` method on `ColorSchemeType`. We need
to add one. The fork's `main.rs` calls `has_plddt()` on the protein before constructing
the App, then passes the result through:

```rust
let has_plddt = protein.has_plddt();
let initial_color = ColorSchemeType::from_cli(&cli.color, has_plddt);
```

### 2e. `has_plddt()` Method (`src/model/protein.rs`)

Add the detection heuristic as a method on `Protein`. In our codebase, this needs to
iterate both chain atoms AND ligand atoms (or, better, only chain atoms -- see Section 3
for rationale).

### 2f. Ribbon Mesh Cache

Our codebase has `mesh_dirty` and `mesh_cache` on App. Any color scheme change already
sets `mesh_dirty = true`, so pLDDT will get correct ribbon mesh regeneration with no
additional work.


## 3. Pitfall Analysis

### 3a. False Positives: Crystallographic B-factors Mimicking pLDDT

**Scenario**: A poorly refined or low-resolution X-ray structure might have B-factors
clustered in 0-100 with a mean above 50.

**Risk assessment**: Moderate. The three-condition heuristic is reasonably robust:
- Well-refined structures: mean B-factor is typically 15-30 A^2 -- fails `mean >= 50`
- Low-resolution structures (>3A): B-factors can be high (40-80), but they also tend
  to have more variance and values exceeding 100 -- fails `in_range >= 0.95`
- Cryo-EM structures: B-factors are often set to 0 or a uniform value -- fails either
  `mean >= 50` or `high_conf >= 0.25`

**Known failure case**: Structures refined with TLS (Translation-Libration-Screw) groups
can have B-factor columns that store only the *residual* B-factor while the TLS
contribution is separate. These residuals can sometimes cluster in plausible pLDDT ranges.
However, this is rare in practice.

**Mitigation**: The CLI override `--color plddt` (force) and `--color bfactor` (force
standard B-factor) let users override the heuristic. We should also consider checking
the file's REMARK records or mmCIF metadata for AlphaFold provenance (see 3d).

### 3b. Structures with Mixed B-factor Semantics

**Scenario**: A hybrid model where some chains are from AlphaFold prediction and others
are experimental, or a structure with deposited AlphaFold models alongside experimental
ligand coordinates.

**Risk**: The heuristic examines ALL atoms globally, so a minority of non-pLDDT chains
could be diluted by the majority, or vice versa.

**Mitigation**: A per-chain variant of `has_plddt()` could be considered, but adds
complexity. For an initial implementation, the global heuristic is sufficient. Users
can always force `--color bfactor` if the auto-detection is wrong.

### 3c. mmCIF Files and Alternative pLDDT Fields

In mmCIF (PDBx/mmCIF) format, pLDDT may be stored in different locations:
- `_atom_site.B_iso_or_equiv` -- the standard B-factor column (most common; AlphaFold
  DB deposits use this)
- `_ma_qa_metric_local.metric_value` -- ModelArchive quality assessment metric. This is
  the "proper" location in the MA (ModelArchive) extension dictionary, but few tools
  read it.

Our parser (`pdbtbx`) reads `B_iso_or_equiv` via `atom.b_factor()`, which covers the
vast majority of AlphaFold PDB and mmCIF files. The MA metric field would require
custom mmCIF parsing and is out of scope for the initial implementation.

### 3d. B-factors Per-Atom vs. Per-Residue

pLDDT is conceptually per-residue, but in PDB files it is written identically to every
atom in the residue (AlphaFold writes the same value for all atoms in a residue).
Our implementation averages across atoms, which correctly collapses per-atom values
that should already be identical. If a file has slight per-atom variation (e.g., from
post-processing), the average still produces the correct residue-level score.

No special handling is needed, but it is worth noting: if we ever display per-atom
pLDDT tooltips, we should show the residue average, not individual atom values.

### 3e. B-factor = 0 Everywhere (Unrefined Structures)

**Scenario**: Many theoretical models, homology models, and unrefined structures set
all B-factors to 0.0.

**Heuristic behavior**: `mean = 0.0` fails `mean >= 50.0`, so `has_plddt()` correctly
returns false. No problem here.

### 3f. AlphaFold3 Multi-Chain Predictions

AlphaFold3 can predict complexes with proteins, nucleic acids, and ligands. In these
files:
- Protein/nucleic acid chains have per-residue pLDDT in the B-factor column
- Ligand atoms may have pLDDT-like confidence values OR default values
- The current parser already uses loose mode fallback for AF3 files
  (`set_level(StrictnessLevel::Loose)` and `set_only_atomic_coords(true)`)

Our ligand atoms already get read with their B-factors. The coloring for ligands under
pLDDT mode should use their own pLDDT if present (AF3 does provide per-atom confidence
for ligands), or fall back to a neutral color.


## 4. Our Improvements Over the Fork

### 4a. Auto-Detect with CLI Override

The fork only allows pLDDT through `from_cli()` with a guard. We should support:
- **Auto-detect default**: When `--color` is not specified, if `has_plddt()` is true,
  default to pLDDT coloring instead of Structure. This gives the best out-of-box
  experience for AlphaFold models.
- **Explicit force**: `--color plddt` should enable pLDDT coloring even if the
  heuristic does not trigger (user knows their file has pLDDT). This means removing
  the `has_plddt` gate in `from_cli()` for the explicit case.
- **Explicit standard**: `--color bfactor` always uses the gradient, even on AF models.

This differs from the fork, where `--color plddt` falls back to Structure if the
heuristic fails. We should trust the user's explicit request.

### 4b. Ligand Coloring Integration

Our codebase has a `Ligand` struct and `ligand_color()`/`ligand_atom_color()` methods
that the fork does not (the fork has no small molecule support). We need to handle pLDDT
mode in these methods:

```rust
ColorSchemeType::Plddt => {
    // Ligands don't have meaningful pLDDT -- use structure-mode colors
    match ligand.ligand_type {
        LigandType::Ligand => Color::Rgb(255, 0, 255),  // magenta
        LigandType::Ion => Color::Rgb(0, 255, 255),     // cyan
    }
},
```

For `ligand_atom_color()`, pLDDT mode should fall through to the ligand base color
(not element coloring) to keep the visual distinction clear.

### 4c. SmallMolecule Chain Exclusion from Detection

Our `MoleculeType` enum includes `SmallMolecule` (used for ligand-only chains). These
chains should be excluded from the `has_plddt()` heuristic because:
- Ligand B-factors may be set to arbitrary values (often 0 or a uniform placeholder)
- Including them could skew the mean downward, causing false negatives
- pLDDT is only defined for polymer chains (protein, RNA, DNA)

Implementation: filter to only `MoleculeType::Protein`, `MoleculeType::RNA`, and
`MoleculeType::DNA` chains in `has_plddt()`.

However, note that our current architecture stores ligands in `Protein.ligands` (a
separate Vec), not as chains with `MoleculeType::SmallMolecule`. So the `SmallMolecule`
variant is not yet used. The `has_plddt()` method should iterate only
`self.chains.iter().flat_map(|c| &c.residues).flat_map(|r| &r.atoms)` and explicitly
exclude `self.ligands`. This is exactly what the fork does (it has no `ligands` field).

### 4d. Provenance Detection Enhancement (Future)

For a more robust detection in the future, we could check PDB REMARK records or mmCIF
metadata for AlphaFold provenance strings:
- PDB: `REMARK 1` or `TITLE` containing "ALPHAFOLD" or "PREDICTED"
- mmCIF: `_ma_protocol_step.software_group_id` or `_software.name` containing
  "AlphaFold"

This is out of scope for the initial implementation but would eliminate false
positives entirely.


## 5. Test Strategy

### 5a. Unit Tests for `plddt_color()`

Test color band boundaries:

| Input Score | Expected Color          | Band     |
|-------------|-------------------------|----------|
| 100.0       | (0, 83, 214)            | >= 90    |
| 95.0        | (0, 83, 214)            | >= 90    |
| 90.0        | (0, 83, 214)            | >= 90    |
| 89.9        | (101, 203, 243)         | 70-89    |
| 80.0        | (101, 203, 243)         | 70-89    |
| 70.0        | (101, 203, 243)         | 70-89    |
| 69.9        | (255, 219, 19)          | 50-69    |
| 60.0        | (255, 219, 19)          | 50-69    |
| 50.0        | (255, 219, 19)          | 50-69    |
| 49.9        | (255, 125, 69)          | < 50     |
| 25.0        | (255, 125, 69)          | < 50     |
| 0.0         | (255, 125, 69)          | < 50     |

Test with a residue that has multiple atoms (all same B-factor) to confirm averaging works.

Test with an empty residue (no atoms) -- should produce orange (avg = 0.0 < 50).

### 5b. Unit Tests for `has_plddt()`

**Positive cases** (should return true):
- Typical AlphaFold output: scores like [95, 92, 88, 76, 67, 54] (mean ~79)
- All atoms at 70.0 (boundary: mean=70, in_range=1.0, high_conf=1.0)
- Mixed high/low: 80% at 90, 20% at 30 (mean=78, high_conf=0.8)

**Negative cases** (should return false):
- Typical crystallographic: [12, 18, 22, 30, 16, 25] (mean ~20.5)
- All zeros: [0, 0, 0, 0] (mean=0)
- High B-factors exceeding 100: [120, 85, 95, 110] (in_range < 0.95)
- Uniform low values: [20, 20, 20, 20] (mean=20)
- Empty protein (no atoms)

**Edge cases**:
- Single atom at 90.0 (should pass all three conditions)
- Mean exactly at threshold: 50.0 (should pass)
- in_range exactly at threshold: 95% in range, 5% out (should pass at 0.95)
- high_conf exactly at threshold: 25% above 70 (should pass at 0.25)

### 5c. Integration with Ligand Exclusion

- Protein with pLDDT-like B-factors on chains + ligands with B-factor=0 should still
  detect as pLDDT (if ligands are excluded from the heuristic)
- Protein with ligands only having high B-factors should NOT trigger pLDDT

### 5d. Color Cycling Tests

- When `has_plddt() == true`: Structure -> Plddt -> Chain -> Element -> BFactor -> Rainbow -> Structure
- When `has_plddt() == false`: Structure -> Chain -> Element -> BFactor -> Rainbow -> Structure (pLDDT skipped)
- From Interface mode: always returns to Structure

### 5e. CLI Parsing Tests

- `--color plddt` with AF model -> Plddt
- `--color plddt` with non-AF model -> Plddt (user override, trust user)
- `--color structure` with AF model -> Structure (user chose explicitly)
- `--color bfactor` with AF model -> BFactor (user wants gradient, not bands)
- `--color invalid` -> Structure (fallback)

### 5f. Ligand Color Under pLDDT Mode

- Ligand type `Ligand` in pLDDT mode -> magenta (same as structure mode)
- Ligand type `Ion` in pLDDT mode -> cyan (same as structure mode)
- `ligand_atom_color()` in pLDDT mode -> uses base ligand color, not element

### 5g. Real File Tests (if example AF files are added)

Consider adding a small AlphaFold PDB to `examples/` for integration testing:
- Verify `has_plddt()` returns true
- Verify default color scheme selection
- Verify color output at specific residues matches expected bands


## 6. Summary of Required Changes

| File                     | Change                                                        |
|--------------------------|---------------------------------------------------------------|
| `src/model/protein.rs`   | Add `has_plddt()` method on `Protein`                        |
| `src/render/color.rs`    | Add `Plddt` variant to `ColorSchemeType`                     |
| `src/render/color.rs`    | Add `plddt_color()` method to `ColorScheme`                  |
| `src/render/color.rs`    | Wire pLDDT into `residue_color()`, `ligand_color()`, `name()`|
| `src/render/color.rs`    | Change `next()` signature to accept `has_plddt: bool`        |
| `src/render/color.rs`    | Add `from_cli()` constructor on `ColorSchemeType`            |
| `src/app.rs`             | Update `cycle_color()` to pass `has_plddt()`                 |
| `src/main.rs`            | Call `has_plddt()`, pass to `from_cli()` and `App::new()`    |
| `src/main.rs`            | Update `--color` help text to include `plddt`                |

Estimated scope: ~80-120 lines of new code (excluding tests), ~20 lines of modified code.

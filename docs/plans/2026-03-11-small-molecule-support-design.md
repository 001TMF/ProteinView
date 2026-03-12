# Design: Small Molecule Support (Ligands + Ions)

**Date:** 2026-03-11
**Branch:** feature/small-molecule-support

## Overview & Motivation

proteinview currently parses all residues through the same pipeline: protein amino acids and nucleotides are recognized, but HETATM records (ligands, cofactors, metal ions) are silently treated as protein residues. This means ligands either render as broken backbone traces or get lumped into the nearest chain with incorrect inter-residue bonding.

Proper small molecule support lets users see drug-protein interactions, cofactor binding sites, and metal coordination -- the most common reason to inspect a PDB structure in practice.

## Scope

### In scope
- Parse HETATM records and classify them as ligands or ions
- Filter out water molecules (HOH, WAT, DOD) as visual noise
- Store ligands separately from polymer chains in the model
- Render ligands as ball-and-stick (all viz modes, both backends)
- Render ions as spheres
- Extend interface analysis to detect ligand binding pockets (polymer residues within contact distance of each ligand)
- Binding pocket summary in the interface panel
- Ligand info in the status bar
- Toggle ligand visibility

### Out of scope
- Water molecule rendering (filtered out)
- Ligand-specific force-field bond tables (use distance-based bonding)
- SDF/MOL2 file format support
- Ligand 2D chemical diagrams
- Modified amino acids (MSE, HYP, etc.) -- these are already handled as regular residues

## Architecture

### Model Changes (`src/model/protein.rs`)

**MoleculeType enum** -- add `SmallMolecule` variant:

```rust
pub enum MoleculeType {
    Protein,
    RNA,
    DNA,
    SmallMolecule, // new
}
```

**LigandType enum** -- distinguish ligands from ions:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LigandType {
    Ligand, // multi-atom small molecule (ATP, HEM, etc.)
    Ion,    // single-atom metal ion (ZN, MG, FE, etc.)
}
```

**Ligand struct** -- a single hetero group:

```rust
#[derive(Debug, Clone)]
pub struct Ligand {
    pub name: String,       // residue name (e.g., "ATP", "HEM", "ZN")
    pub chain_id: String,   // parent chain ID
    pub seq_num: i32,       // residue sequence number
    pub atoms: Vec<Atom>,   // all atoms in this hetero residue
    pub ligand_type: LigandType,
}
```

**Protein struct** -- add ligands field:

```rust
pub struct Protein {
    pub name: String,
    pub chains: Vec<Chain>,
    pub ligands: Vec<Ligand>,  // new -- all non-water HETATM groups
}
```

**Protein methods** -- update to include ligand atoms:

- `atom_count()`: add `self.ligands.iter().flat_map(|l| &l.atoms).count()`
- `bounding_radius()`: include ligand atom positions in the max-distance calculation (use all atoms, not just backbone)
- `center()`: include ligand atoms in centroid computation and apply the same translation
- Add `ligand_count() -> usize`
- Add `ligand_atom_count() -> usize`

### Common Ions List

Single-atom HETATM residues with names matching common metal/halide ions:

```rust
const COMMON_IONS: &[&str] = &[
    "ZN", "MG", "CA", "FE", "MN", "CO", "CU", "NI", "CD",
    "NA", "K", "CL", "BR", "I", "F",
    "HG", "PT", "AU", "AG", "PB",
];
```

Classification rule: if the HETATM residue has exactly 1 non-hydrogen atom, or its name is in `COMMON_IONS`, it is an `Ion`. Otherwise it is a `Ligand`.

### Water Residues to Filter

```rust
const WATER_NAMES: &[&str] = &["HOH", "WAT", "DOD", "H2O", "OH2"];
```

### Parser Changes (`src/parser/pdb.rs`)

pdbtbx 0.12 exposes `Atom::hetero() -> bool` which returns `true` for HETATM records. The parser restructuring:

1. Iterate chains as before for polymer residues (ATOM records)
2. In a second pass (or within the same loop), collect HETATM residues:
   - For each `chain.residues()`, check if all atoms in the residue are `hetero()`
   - Skip residues with names matching `WATER_NAMES`
   - Classify as `Ion` or `Ligand` using the rules above
   - Collect into `Vec<Ligand>`
3. Exclude HETATM residues from the polymer `Chain.residues` so they do not break backbone tracing or secondary structure assignment

**Implementation approach**: During the existing residue iteration loop, check `atom.hetero()` on the first atom of each residue. If all atoms are hetero, divert the residue into the ligands list instead of the chain's residues list. This avoids a second pass.

```
for residue in chain.residues() {
    let atoms: Vec<pdbtbx::Atom> = residue.atoms().collect();
    let all_hetero = atoms.iter().all(|a| a.hetero());
    let res_name = residue.name().unwrap_or("UNK");

    if all_hetero && !WATER_NAMES.contains(&res_name) {
        // -> build Ligand, push to ligands vec
    } else if !all_hetero {
        // -> build Residue, push to chain residues vec (existing path)
    }
    // else: water -- skip entirely
}
```

### Atom Struct Extension

Add an `is_hetero` field to `Atom`:

```rust
pub struct Atom {
    pub name: String,
    pub element: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub b_factor: f64,
    pub is_backbone: bool,
    pub is_hetero: bool, // new -- true for HETATM records
}
```

This field is informational and helps downstream code (e.g., color schemes) without needing to cross-reference the ligands list.

## Rendering Strategy

### Ligand Rendering Primitives

**Ball-and-stick** (ligands):
- Atoms: spheres with radius proportional to element (scaled to ~0.4 A for display)
- Bonds: cylinders/lines between atoms within 1.9 A (same distance threshold used for intra-residue wireframe bonds)
- Bond radius: ~0.15 A (thinner than atom spheres)

**Sphere** (ions):
- Single sphere with radius ~0.8 A
- Colored by element (CPK)

### Per Viz Mode

| Viz Mode | Polymer Rendering | Ligand Rendering |
|---|---|---|
| Backbone | CA/C4' trace (thick lines) | Ball-and-stick overlay |
| Cartoon | Ribbon mesh | Ball-and-stick overlay |
| Wireframe | All-atom bonds | All-atom bonds (same as wireframe, already natural) |

In Backbone and Cartoon modes, ligands are always rendered as ball-and-stick regardless of the polymer representation. This is the standard approach in molecular viewers (PyMOL, ChimeraX).

In Wireframe mode, ligands render identically to polymer residues -- all atoms with distance-based bonds -- which already works correctly since ligands are self-contained residues.

### HD Mode (`src/render/hd.rs`)

Add `render_ligands_fb()` function called after polymer rendering in all three viz modes:

```rust
fn render_ligands_fb(
    fb: &mut Framebuffer,
    protein: &Protein,
    camera: &Camera,
    color_scheme: &ColorScheme,
    half_w: f64,
    half_h: f64,
) {
    for ligand in &protein.ligands {
        match ligand.ligand_type {
            LigandType::Ion => {
                // Single sphere, larger radius (~3.0 px)
                render_ion_sphere(fb, ligand, camera, color_scheme, half_w, half_h);
            }
            LigandType::Ligand => {
                // Ball-and-stick: atom spheres + bond sticks
                render_ball_and_stick(fb, ligand, camera, color_scheme, half_w, half_h);
            }
        }
    }
}
```

Ball-and-stick rendering uses existing `Framebuffer` primitives:
- `draw_circle_z()` for atom spheres (radius 2.5-3.5 px depending on element)
- `draw_thick_line_3d()` for bond sticks (thickness 1.5 px)

Ion spheres use `draw_circle_z()` with a larger radius (4.0-5.0 px).

The z-buffer ensures correct depth ordering with polymer geometry. The post-pass `apply_depth_tint()` applies uniformly.

### Braille Mode (`src/render/braille.rs`)

Add `render_ligands()` function called after polymer rendering:

- Ligand atoms: `ctx.draw()` points (braille dots) at projected positions
- Ligand bonds: `ctx.draw(&Line { ... })` for bonded atom pairs within 1.9 A
- Ions: single braille dot at the ion position

This mirrors the wireframe rendering approach but uses the ligand data.

### Cartoon Ribbon Mesh (`src/render/ribbon.rs`)

No changes needed. Ligands are not part of the ribbon mesh. They are rendered as a separate overlay pass in `hd.rs` and `braille.rs`.

## Color Scheme Extensions

### Element Scheme (`ColorSchemeType::Element`)

Already works correctly for ligands via `element_color()`. CPK coloring applies to any atom regardless of residue type. No changes needed.

### Structure Scheme (`ColorSchemeType::Structure`)

Add ligand/ion cases. For ligands, since they are not `Residue` objects, add a `ligand_color()` method to `ColorScheme`:

```rust
pub fn ligand_color(&self, ligand: &Ligand) -> Color {
    match self.scheme_type {
        ColorSchemeType::Structure => match ligand.ligand_type {
            LigandType::Ligand => Color::Rgb(255, 0, 255),  // magenta
            LigandType::Ion => Color::Rgb(0, 255, 255),     // cyan
        },
        ColorSchemeType::Chain => self.chain_color_by_id(&ligand.chain_id),
        ColorSchemeType::Element => Color::Rgb(144, 144, 144), // default, overridden per-atom
        ColorSchemeType::BFactor => self.bfactor_color_from_atoms(&ligand.atoms),
        ColorSchemeType::Rainbow => Color::Rgb(255, 0, 255),  // fixed, no sequence position
        ColorSchemeType::Interface => self.ligand_interface_color(ligand),
    }
}

pub fn ligand_atom_color(&self, atom: &Atom, ligand: &Ligand) -> Color {
    match self.scheme_type {
        ColorSchemeType::Element => Self::element_color(atom),
        _ => self.ligand_color(ligand),
    }
}
```

### Interface Scheme

When interface mode is active and binding pocket analysis is available:
- Ligand atoms: bright white `[255, 255, 255]` to stand out
- Binding pocket residues (polymer residues contacting the ligand): bright magenta `[255, 100, 255]`
- Non-pocket residues: dimmed as in current interface scheme

## Interface Analysis Extensions

### New: Ligand Binding Pocket Detection (`src/model/interface.rs`)

Add a `LigandContact` struct and binding pocket analysis:

```rust
/// A contact between a ligand and a polymer residue.
#[derive(Debug, Clone)]
pub struct LigandContact {
    pub ligand_idx: usize,     // index into protein.ligands
    pub chain_idx: usize,      // index into protein.chains
    pub residue_idx: usize,    // index into chain.residues
    pub min_distance: f64,     // minimum heavy-atom distance
}

/// Binding pocket analysis for all ligands.
#[derive(Debug, Clone)]
pub struct BindingPocketAnalysis {
    pub contacts: Vec<LigandContact>,
    /// Per-ligand: set of (chain_idx, residue_idx) forming the binding pocket.
    pub pockets: Vec<HashSet<(usize, usize)>>,
}
```

Add `analyze_binding_pockets()`:

```rust
pub fn analyze_binding_pockets(protein: &Protein, cutoff: f64) -> BindingPocketAnalysis {
    // For each ligand, find all polymer residues with any heavy atom
    // within cutoff distance of any ligand heavy atom.
    // Use the same dist_sq() helper and hydrogen-skipping logic
    // as analyze_interface().
}
```

Default cutoff: 4.5 A (same as chain-chain interface).

### InterfaceAnalysis Integration

Extend `InterfaceAnalysis` to optionally hold `BindingPocketAnalysis`:

```rust
pub struct InterfaceAnalysis {
    pub contacts: Vec<Contact>,
    pub interface_residues: HashSet<(usize, usize)>,
    pub chain_interface_counts: Vec<usize>,
    pub total_interface_residues: usize,
    pub binding_pockets: Option<BindingPocketAnalysis>, // new
}
```

Update `summary()` to append ligand binding information when pockets are present:

```
Ligand HEM (Chain A): 12 pocket residues, min dist 2.1A
Ligand ATP (Chain B): 8 pocket residues, min dist 2.4A
Ion ZN (Chain A): 4 coordinating residues, min dist 2.0A
```

## UI Changes

### Status Bar

Add ligand count to the info line when ligands are present:

```
│ Chain A │ 248 res │ 3 ligands │ Cartoon │ Structure │ HD
```

### Interface Panel

When binding pockets are detected, add a "Ligand Contacts" section below the existing chain-chain interface section:

```
 Ligand Contacts
 HEM (A:154): 12 res, min 2.1A
   A:HIS87 (2.1A), A:CYS102 (2.3A), ...
 ZN (A:301): 4 res, min 2.0A
   A:HIS69 (2.0A), A:CYS96 (2.1A), ...
```

### Keybindings

Add `g` key to toggle ligand visibility:

| Key | Action |
|---|---|
| `g` | Toggle ligand/ion visibility |

### App State (`src/app.rs`)

Add fields:

```rust
pub struct App {
    // ... existing fields ...
    pub show_ligands: bool,            // new -- default true
    pub binding_pockets: Option<BindingPocketAnalysis>, // new
}
```

`show_ligands` defaults to `true` so ligands are visible on load.

## Implementation Phases

### Phase 1: Model & Parser (no rendering changes)
**Files:** `src/model/protein.rs`, `src/parser/pdb.rs`
**Dependencies:** None

1. Add `SmallMolecule` to `MoleculeType`, add `LigandType` enum
2. Add `Ligand` struct
3. Add `is_hetero` to `Atom`
4. Add `ligands: Vec<Ligand>` to `Protein`
5. Update `Protein::center()`, `bounding_radius()`, `atom_count()` to include ligands
6. Update parser to detect HETATM via `atom.hetero()`, filter waters, classify ligands vs ions, populate `protein.ligands`
7. Exclude HETATM residues from polymer chain residue lists

### Phase 2: Color Scheme Extensions
**Files:** `src/render/color.rs`
**Dependencies:** Phase 1

1. Add `ligand_color()` and `ligand_atom_color()` methods to `ColorScheme`
2. Add Structure scheme colors for ligands (magenta) and ions (cyan)
3. Add `element_color()` entries for additional elements common in ligands

### Phase 3: HD Rendering
**Files:** `src/render/hd.rs`
**Dependencies:** Phase 1, Phase 2

1. Add `render_ligands_fb()` with ball-and-stick and ion sphere rendering
2. Call from `render_hd_framebuffer()` after polymer rendering, gated on `show_ligands`
3. Z-buffer handles depth ordering automatically

### Phase 4: Braille Rendering
**Files:** `src/render/braille.rs`
**Dependencies:** Phase 1, Phase 2

1. Add `render_ligands()` with point/line rendering for ligands and ions
2. Call from `render_protein()` after polymer rendering, gated on `show_ligands`

### Phase 5: Interface Analysis & Binding Pockets
**Files:** `src/model/interface.rs`
**Dependencies:** Phase 1

1. Add `LigandContact`, `BindingPocketAnalysis` structs
2. Implement `analyze_binding_pockets()`
3. Extend `InterfaceAnalysis` with optional `binding_pockets`
4. Update `summary()` to include ligand pocket info

### Phase 6: App State & UI
**Files:** `src/app.rs`, `src/main.rs`, `src/ui/statusbar.rs`, `src/ui/interface_panel.rs`, `src/ui/help_overlay.rs`, `src/ui/helpbar.rs`
**Dependencies:** Phase 1-5

1. Add `show_ligands` and `binding_pockets` to `App`
2. Compute binding pockets in `App::new()`
3. Add `g` keybinding to toggle ligand visibility
4. Update status bar to show ligand count
5. Update interface panel to show binding pocket info
6. Update help overlay with new keybinding

## Files Changed Summary

| File | Change |
|---|---|
| `src/model/protein.rs` | `SmallMolecule` variant, `LigandType` enum, `Ligand` struct, `is_hetero` on `Atom`, `ligands` on `Protein`, update `center()`/`bounding_radius()`/`atom_count()` |
| `src/parser/pdb.rs` | HETATM detection via `atom.hetero()`, water filtering, ligand/ion classification, exclude HETATM from chain residues |
| `src/render/color.rs` | `ligand_color()`, `ligand_atom_color()`, Structure scheme ligand/ion colors, extended CPK table |
| `src/render/hd.rs` | `render_ligands_fb()` with ball-and-stick + ion spheres, called gated on `show_ligands` |
| `src/render/braille.rs` | `render_ligands()` with point/line rendering, called gated on `show_ligands` |
| `src/model/interface.rs` | `LigandContact`, `BindingPocketAnalysis`, `analyze_binding_pockets()`, pocket info in `summary()` |
| `src/app.rs` | `show_ligands`, `binding_pockets` fields, `toggle_ligands()` method |
| `src/main.rs` | `g` keybinding |
| `src/ui/statusbar.rs` | Ligand count display |
| `src/ui/interface_panel.rs` | Binding pocket section |
| `src/ui/help_overlay.rs` | `g` keybinding entry |
| `src/ui/helpbar.rs` | Updated hint text |

## Testing Strategy

### Unit Tests

**Model (`src/model/protein.rs`):**
- `Ligand` construction and field access
- `LigandType` classification: multi-atom -> `Ligand`, single-atom metal -> `Ion`
- `Protein::atom_count()` includes ligand atoms
- `Protein::bounding_radius()` includes ligand atoms
- `Protein::center()` shifts ligand atom coordinates

**Parser (`src/parser/pdb.rs`):**
- Water residues (HOH, WAT) are filtered out
- HETATM residues are excluded from chain residue lists
- Single-atom metal HETATM -> `Ion`
- Multi-atom HETATM -> `Ligand`
- Mixed ATOM/HETATM residues handled correctly (edge case)
- `classify_chain_type()` still returns correct types (regression)

**Color (`src/render/color.rs`):**
- `ligand_color()` returns magenta for `Ligand` in Structure scheme
- `ligand_color()` returns cyan for `Ion` in Structure scheme
- `ligand_atom_color()` returns CPK color in Element scheme

**Interface (`src/model/interface.rs`):**
- `analyze_binding_pockets()` detects residues within cutoff of ligand atoms
- Hydrogen atoms are skipped in distance calculations
- Empty ligands list produces empty pockets
- Pocket summary formatting

### Integration Tests

- Load PDB files and verify ligand parsing round-trip
- Load a PDB with known ligands, verify correct count and classification
- Verify existing tests pass (no regressions in polymer handling)

### Manual Testing

- Fetch 1HHO (hemoglobin with HEM): verify HEM renders as ball-and-stick, iron as sphere
- Fetch 4HHB (deoxyhemoglobin): verify multiple HEM groups across chains
- Fetch 1ATP (ATP-bound kinase): verify ATP ball-and-stick + Mg ion sphere
- Toggle ligand visibility with `g` key
- Verify binding pocket highlighting in interface mode
- Test all three viz modes (Backbone, Cartoon, Wireframe) with ligands
- Test both rendering backends (Braille, HD)

<p align="center">
  <b>P R O T E I N V I E W</b>
</p>

<p align="center">
  <em>Explore molecular structures in your terminal</em>
</p>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT"></a>
  <a href="https://www.rust-lang.org/"><img src="https://img.shields.io/badge/Rust-1.85%2B-orange.svg" alt="Rust"></a>
  <img src="https://img.shields.io/badge/version-0.3.0-green.svg" alt="Version">
  <a href="https://github.com/001TMF/ProteinView/pulls"><img src="https://img.shields.io/badge/PRs-welcome-brightgreen.svg" alt="PRs Welcome"></a>
</p>

<p align="center">
  <img src="assets/hero-histone.png" alt="Nucleosome core particle with histone proteins and DNA rendered in FullHD mode" width="700">
</p>

<p align="center">
  <sub>Nucleosome core particle — histone octamer wrapped in DNA, rendered with Kitty graphics protocol</sub>
</p>

---

Terminal molecular structure viewer — load, rotate, and explore proteins, nucleic acids, and small molecules from PDB/CIF files right in your terminal. No browser, no GUI, no dependencies.

## Features

- **3-tier render modes** — Braille, HD, and FullHD (Sixel/Kitty) with automatic SSH detection
- **PNG-compressed Kitty protocol** — ~10-20x smaller than raw RGBA, making FullHD viable over SSH
- **Cartoon ribbon visualization** — Lambert-shaded ribbons with depth fog for helices, sheets, and coils
- **RNA/DNA support** — backbone, wireframe, and cartoon modes with base-type coloring
- **Small molecule rendering** — ligands as ball-and-stick, ions as spheres
- **Interface analysis** — inter-chain contacts, binding pockets, and interaction visualization (H-bonds, salt bridges, hydrophobic contacts)
- **7 color schemes** — structure, chain, element (CPK), B-factor, rainbow, pLDDT (AlphaFold)
- **Interactive controls** — vim-style rotation, zoom, pan with auto-rotation
- **PDB & mmCIF** — both formats supported, with RCSB PDB fetch (`--fetch`)
- **Single static binary** — zero runtime dependencies

## Render Modes

Three rendering tiers to match your terminal and connection:

<p align="center">
  <img src="assets/render-modes-grid.png" alt="Braille vs HD vs FullHD rendering comparison" width="700">
</p>

<p align="center">
  <sub>Left: Braille (works everywhere, including SSH/tmux) · Middle: HD (Lambert-shaded braille) · Right: FullHD (Kitty pixel graphics)</sub>
</p>

| Mode | Key | Quality | SSH Performance |
|------|-----|---------|-----------------|
| **Braille** | default | Text-based, monochrome per cell | Excellent |
| **HD** | `m` | Shaded braille with lighting + depth fog | Excellent |
| **FullHD** | `M` | Sixel/Kitty pixel graphics | Good (PNG compressed) |

`--hd` is SSH-aware: defaults to HD over SSH, FullHD locally. Use `--fullhd` to force pixel graphics.

## Visualization Modes

<p align="center">
  <img src="assets/viz-modes-grid.png" alt="Cartoon, Wireframe, and Backbone visualization modes" width="700">
</p>

<p align="center">
  <sub>Left: Cartoon (ribbon) · Middle: Wireframe (all-atom) · Right: Backbone (CA trace)</sub>
</p>

| Mode | Description |
|------|-------------|
| **Cartoon** | Smooth ribbon rendering — helices, beta-sheets, and coils with Lambert shading. Default. |
| **Wireframe** | All-atom bonds including inter-residue peptide and phosphodiester linkages. |
| **Backbone** | CA trace (proteins) / C4' trace (nucleic acids) with spheres and thick connecting lines. |

## Interface Analysis & Interactions

<p align="center">
  <img src="assets/interface-grid.png" alt="Interface analysis with interaction visualization" width="700">
</p>

<p align="center">
  <sub>Left: Interface residue coloring with sidebar panel · Right: Dashed interaction lines (H-bonds, salt bridges, hydrophobic contacts)</sub>
</p>

Press `f` to toggle interface mode — highlights contact residues across chain boundaries with a detailed sidebar. Press `I` to overlay interaction lines:

| Color | Interaction | Distance |
|-------|-------------|----------|
| Cyan | Hydrogen bond | &le; 3.5 &Aring; |
| Red | Salt bridge | &le; 4.0 &Aring; |
| Yellow | Hydrophobic contact | &le; 4.5 &Aring; |
| Gray | Other | &le; 4.5 &Aring; |

## Nucleic Acids

<p align="center">
  <img src="assets/dna-element.png" alt="B-DNA double helix with element (CPK) coloring" width="500">
</p>

<p align="center">
  <sub>B-DNA dodecamer in wireframe mode with CPK element coloring</sub>
</p>

Full support for DNA and RNA structures — backbone traces, wireframe bonds, and cartoon ribbons with nucleotide base-type coloring (A=red, U/T=blue, G=green, C=yellow).

## AlphaFold & pLDDT

<p align="center">
  <img src="assets/plddt-cartoon.png" alt="AlphaFold prediction with pLDDT confidence coloring" width="500">
</p>

<p align="center">
  <sub>AlphaFold prediction with pLDDT confidence coloring — blue (high confidence) to orange/yellow (low confidence)</sub>
</p>

Automatically detects AlphaFold structures and offers pLDDT confidence coloring. Cycle through color schemes with `c`.

## Installation

```bash
# Basic install
cargo install --path .

# With RCSB PDB fetch support
cargo install --path . --features fetch

# Update an existing installation
cargo install --path . --force
```

## Quick Start

```bash
# View a local PDB file
proteinview examples/1AOI.pdb

# HD mode (fast text-based shading)
proteinview examples/4HHB.pdb --hd

# FullHD pixel mode (Kitty/Sixel terminals)
proteinview examples/4HHB.pdb --fullhd

# Fetch from RCSB PDB
proteinview --fetch 1UBQ

# Choose color scheme and visualization
proteinview examples/1UBQ.pdb --color rainbow --mode wireframe
```

## Keybindings

| Key | Action |
|-----|--------|
| `h`/`l` | Rotate Y |
| `j`/`k` | Rotate X |
| `u`/`i` | Roll |
| `+`/`-` | Zoom |
| `w`/`a`/`s`/`d` | Pan |
| `r` | Reset view |
| `c` | Cycle color scheme |
| `v` | Cycle viz mode |
| `m` | Braille / HD |
| `M` | HD / FullHD |
| `f` | Interface analysis |
| `I` | Interface interactions |
| `g` | Toggle ligands |
| `[`/`]` | Prev/next chain |
| `Space` | Auto-rotate |
| `?` | Help |
| `q` | Quit |

## Color Schemes

| Scheme | Description |
|--------|-------------|
| **Structure** | Helix (red), sheet (yellow), coil (green). Default. |
| **Chain** | Distinct color per chain. |
| **Element** | CPK coloring (C, N, O, S, P, metals). |
| **B-factor** | Blue (rigid) to red (flexible). |
| **Rainbow** | N-terminus (blue) to C-terminus (red). |
| **pLDDT** | AlphaFold confidence (blue=high, orange=low). |

## Terminal Support

| Terminal | Braille | HD | FullHD |
|----------|---------|-----|--------|
| Any Unicode terminal | Yes | Yes | -- |
| Kitty | Yes | Yes | Yes (PNG) |
| WezTerm | Yes | Yes | Yes (Sixel) |
| iTerm2 | Yes | Yes | Yes |
| foot | Yes | Yes | Yes (Sixel) |
| tmux/screen | Yes | Yes | -- |

## Building

```bash
cargo build --release

# With RCSB fetch support
cargo build --release --features fetch
```

## Contributing

Contributions are welcome! Here's how to get started:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make your changes and add tests
4. Run `cargo test` to verify
5. Open a pull request against `develop`

Please open an issue first for major changes to discuss the approach.

## License

[MIT](LICENSE)

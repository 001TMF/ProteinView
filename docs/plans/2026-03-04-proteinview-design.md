# ProteinView — TUI Protein Structure Viewer

**Date:** 2026-03-04
**Status:** Approved

## Overview

A terminal-based protein structure viewer written in Rust. Load PDB/CIF files and interactively view 3D protein structures directly in the terminal over SSH. Two rendering modes: braille characters (default, universal) and HD pixel graphics (sixel/kitty, opt-in).

## CLI Interface

```bash
# Default (braille)
proteinview myprotein.pdb

# HD pixel mode
proteinview myprotein.pdb --hd

# mmCIF support
proteinview structure.cif

# Fetch from RCSB PDB
proteinview --fetch 1UBQ

# Options
--hd / --pixel        Sixel/Kitty pixel rendering (auto-detected)
--color <scheme>      structure|chain|element|bfactor|rainbow (default: structure)
--mode <mode>         backbone|ballandstick|wireframe (default: backbone)
--fetch <PDB_ID>     Download from RCSB and display
```

## TUI Layout

```
╭─── ProteinView ─── 1UBQ: Ubiquitin ──────────────────────────────╮
│                                                                    │
│                       3D VIEWPORT                                  │
│                    (braille or HD)                                  │
│                                                                    │
├────────────────────────────────────────────────────────────────────┤
│ Chain A │ 76 res │ Backbone │ Structure Color │ Braille           │
╰── h/l: rotate Y  j/k: rotate X  +/-: zoom  wasd: pan  ?: help ──╯
```

- **Title bar** — retro ASCII header with protein name/PDB ID
- **Viewport** — 3D rendering area (~85% of screen)
- **Status bar** — chain, residue count, mode, color scheme, render mode
- **Keybinding bar** — context-sensitive hints

## Keybindings

| Key | Action |
|-----|--------|
| h/l | Rotate Y-axis |
| j/k | Rotate X-axis |
| u/i | Rotate Z-axis (roll) |
| +/- | Zoom in/out |
| w/a/s/d | Pan |
| r | Reset view |
| c | Cycle color scheme |
| v | Cycle visualization mode |
| m | Toggle braille/HD |
| [/] | Previous/next chain |
| space | Toggle auto-rotation |
| ? | Help overlay |
| q | Quit |

## Architecture

```
proteinview/
├── src/
│   ├── main.rs              # CLI parsing (clap), app entry
│   ├── app.rs               # App state, event loop, mode management
│   ├── parser/
│   │   ├── mod.rs
│   │   ├── pdb.rs           # PDB/CIF loading via pdbtbx
│   │   └── fetch.rs         # RCSB PDB download
│   ├── model/
│   │   ├── mod.rs
│   │   ├── protein.rs       # Protein data: atoms, residues, chains
│   │   └── secondary.rs     # Secondary structure assignment
│   ├── render/
│   │   ├── mod.rs
│   │   ├── camera.rs        # Rotation, zoom, pan, projection
│   │   ├── braille.rs       # Braille renderer (ratatui Canvas)
│   │   ├── pixel.rs         # Sixel/Kitty pixel renderer
│   │   └── color.rs         # Color schemes (CPK, structure, chain, etc.)
│   ├── ui/
│   │   ├── mod.rs
│   │   ├── header.rs        # Retro ASCII header/title bar
│   │   ├── viewport.rs      # Main 3D viewport widget
│   │   ├── statusbar.rs     # Status bar
│   │   ├── helpbar.rs       # Keybinding hints
│   │   └── help_overlay.rs  # Full help popup
│   └── event.rs             # Keyboard/mouse event handling
├── Cargo.toml
└── assets/
    └── header.txt           # Retro ASCII art header
```

## Tech Stack

- **ratatui** + **crossterm** — TUI framework + terminal backend
- **pdbtbx** — PDB/mmCIF parsing with R*-tree spatial indexing
- **clap** — CLI argument parsing
- **reqwest** (optional) — HTTP for --fetch

## Rendering

### Braille Mode (default)
- Unicode braille characters (U+2800..U+28FF), 2x4 dots per character cell
- ~160x96 effective resolution in 80x24 terminal
- One foreground color per character cell
- Uses ratatui's built-in Canvas widget with Marker::Braille

### HD/Pixel Mode (--hd flag)
- Auto-detects sixel or kitty graphics protocol support
- Full RGBA pixel-level rendering
- Renders to off-screen buffer, outputs via escape sequences
- Falls back to braille if no protocol detected

### 3D Pipeline
1. Load atoms from PDB/CIF via pdbtbx
2. Extract backbone (C-alpha atoms) or all atoms based on mode
3. Apply rotation matrices (Rx, Ry, Rz) based on camera state
4. Orthographic projection to 2D screen coordinates
5. Z-buffer for depth sorting
6. Color by selected scheme (secondary structure, chain, element, etc.)
7. Rasterize to braille dots or pixel buffer

## Color Schemes

### Secondary Structure (default)
- Helix: red/magenta (#FF0080)
- Sheet: yellow/gold (#FFC800)
- Coil: green (#00CC00)
- Turn: pale blue (#6080FF)

### Element (CPK)
- C: gray (#909090), N: blue (#3050F8), O: red (#FF0D0D)
- S: yellow (#FFFF30), P: orange (#FF8000), H: white (#FFFFFF)

### Chain
- Each chain gets a distinct color from a curated palette

### B-factor
- Blue (low mobility) → Red (high mobility) gradient

### Rainbow
- N-terminus (blue) → C-terminus (red) by residue position

## Data Flow

```
PDB/CIF file → pdbtbx parser → Protein model (atoms, residues, chains)
                                       ↓
                              Camera (rotation, zoom, pan)
                                       ↓
                              3D → 2D projection + z-buffer
                                       ↓
                        ┌──────────────┴──────────────┐
                   Braille renderer            Pixel renderer
                  (ratatui Canvas)          (sixel/kitty escapes)
                        ↓                          ↓
                   Terminal output            Terminal output
```

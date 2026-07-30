# ProteinView for Oh My Pi

This plugin adds an essential `render_protein` tool and a read-approved
`proteinview_control` tool to Oh My Pi. A request such as:

```text
show me 1UBQ PDB protein
```

downloads the structure from RCSB and opens a persistent live ProteinView panel
directly above the editor. The panel remains part of Oh My Pi's layout while
messages scroll and updates its mounted image in place instead of adding a new
transcript image per frame. An exact 1920×1080 FullHD snapshot is rendered only
when the live panel is unavailable or fails to start.

## Requirements

- Oh My Pi 17 or newer.
- Bun 1.3.14 or newer.
- An Oh My Pi TUI build with mutable `Image.setData(...)` support.
- A ProteinView build that supports the snapshot flags and the v1
  `--panel-server` NDJSON protocol.

Install ProteinView from this checkout:

```bash
cargo install --path .
```

The plugin invokes `proteinview` from `PATH`. To select another executable,
set `PROTEINVIEW_BIN` to its absolute path.

ProteinView's optional Cargo `fetch` feature is not required. The plugin
downloads PDB entries itself so it can validate IDs, enforce a response-size
limit and timeout, and forward Oh My Pi cancellation.

## Install

From the ProteinView repository:

```bash
omp plugin link ./integrations/oh-my-pi
```

Restart Oh My Pi after linking the plugin.

## Inputs

The tool accepts either:

- A four-character RCSB PDB ID beginning with `1`–`9`, such as `1UBQ`.
- A `.pdb`, `.cif`, `.mmcif`, or `.xyz` file inside the current Oh My Pi
  workspace.

Local files are resolved through `realpath`, required to remain within the
workspace, size-capped, and copied into a private temporary directory before
ProteinView runs. PDB downloads use a fixed RCSB URL, a strict ID grammar,
redirect rejection, a 30-second timeout, and a 64 MiB cap.

The model may choose `cartoon`, `backbone`, or `wireframe` rendering and any
ProteinView color scheme exposed by the tool. When no mode is requested,
ProteinView uses cartoon for ordinary structures and automatically falls back
to backbone for very large structures. Snapshot dimensions are fixed at
1920×1080 and are not model-controlled. Because the plugin runs the
ProteinView executable, Oh My Pi classifies it as an `exec`-approval tool.
XYZ inputs preserve ProteinView's Element/Wireframe defaults unless the request
explicitly selects another scheme or mode.

Biological view controls are conversational too. The tool can color by chain,
show pLDDT or B-factor confidence, focus the interface for a named chain,
overlay ProteinView's classified hydrogen-bond/salt-bridge/hydrophobic contact
lines, show or hide ligands and ions, and assign exact RGB colors to polymer
residues. For example, “show the interface for chain A with interactions, hide
ligands, and color residue A:42 red” maps to `interface_chain: "A"`,
`show_interactions: true`, `show_ligands: false`, and:

```json
{
  "residue_colors": [
    {"chain": "A", "residue_number": 42, "color": "FF0000"}
  ]
}
```

Insertion-coded residues are exact: residue `42[A]` uses
`{"chain":"A","residue_number":42,"insertion_code":"A","color":"FF0000"}`;
omitting `insertion_code` selects only the blank insertion code, never every
residue numbered 42. Up to 256 unique residue colors can be supplied. Colors
are six hexadecimal RGB digits and ProteinView applies them above every
regular, element, pLDDT, and interface palette. Ligands are intentionally not
selected by this polymer-residue API.
Interface highlighting uses its own green/orange palette, so it is intentionally
mutually exclusive with regular `color` schemes. Exact residue colors may still
be layered over an interface view.

## Live panel

The extension starts one persistent ProteinView process before any fallback
snapshot is attempted. ProteinView parses the private input before the render
tool removes it, then owns molecular state and atomically replaces a separate
private `frame.png`. Commands use request IDs, strict size-capped NDJSON v1
records, timeouts, serialized writes, bounded stderr, and a validated PNG read
after each successful frame revision.

The panel is full-width and mounted with `placement: "aboveEditor"`. It is
bounded to roughly 45% of terminal rows, reserves room for the editor, and
collapses to two text rows on small terminals. The source framebuffer remains
FullHD (1920×1080); Oh My Pi scales that source into the available terminal
cells.

Press `F8` to focus or release the viewer. While focused:

- Arrow keys rotate; `+`/`-` zoom.
- `r` resets and `f` fits the structure.
- `c` cycles color, `v` cycles visualization, and `[`/`]` select a chain.
- `i` toggles interface coloring, `n` toggles interface contacts, and `l`
  toggles ligands.
- `a` toggles auto-rotation.
- `Esc` or `Tab` returns to chat; `q` closes the panel.

The agent uses `proteinview_control` for the same operations, plus pan, exact
color/mode selection, explicit boolean settings, state inspection, and atomic
`set_residue_colors` / `clear_residue_colors` commands. Presentation commands
are classified as `read` approval because they only change the private live
view.

## Snapshot fallback

When the live panel opens, the tool emits no transcript image. If no interactive
handler is available (for example print, RPC, ACP, or headless execution), or if
live-panel startup fails, the tool renders one exact FullHD PNG and returns it in
`details.images`, not in the model-facing tool content. Oh My Pi therefore
displays the fallback inline without adding its bytes to model context. Normal
image settings still apply:

- `terminal.showImages`
- `tui.maxInlineImageColumns`
- `tui.maxInlineImageRows`
- `tui.maxInlineImages`

When inline graphics are unavailable or disabled, Oh My Pi shows the textual
render summary and its normal image placeholder. The plugin never launches
ProteinView's alternate-screen TUI and never emits raw terminal graphics
escapes; all graphics go through Oh My Pi's image component and budget.

## Test

```bash
cd integrations/oh-my-pi
bun test
```

The tests mock Oh My Pi's custom-tool API, panel process/transport, TUI, and
network access; they do not need a running Oh My Pi session or installed
ProteinView binary.

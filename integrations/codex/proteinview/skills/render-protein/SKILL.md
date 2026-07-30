---
name: render-protein
description: Render and display protein structures with ProteinView FullHD pixel graphics. Use whenever a user asks to show, view, visualize, inspect, or render a Protein Data Bank/PDB identifier (for example, "show me 1UBQ PDB protein") or a local .pdb/.cif/.mmcif structure file.
---

# Render a protein

Use the bundled deterministic helper, then pass its output to Codex's built-in `view_image` tool. This produces ProteinView's FullHD pixel mode at a 1920×1080 source resolution; it is not ProteinView's HD text-cell renderer. The helper also writes private, bounded local launch metadata beside the PNG. A live-viewer-enabled Codex build uses that metadata to open a persistent ProteinView pane while keeping chat visible; other Codex builds continue to show the static PNG.

1. Set `SKILL_DIR` to the absolute directory containing this `SKILL.md`.
2. For a PDB identifier, run:

   ```sh
   python3 "$SKILL_DIR/scripts/render_protein.py" --pdb 1UBQ
   ```

   Substitute the requested four-character ID. Do not fetch from an arbitrary host.

3. For a local structure inside the current workspace, run:

   ```sh
   python3 "$SKILL_DIR/scripts/render_protein.py" --file "/absolute/path/to/model.cif"
   ```

4. Read the JSON object printed by the helper. Call `view_image` with its absolute `path` and `detail: "original"`. In the live-enabled build this also opens the interactive pane; F8 moves focus between it and chat.
5. Tell the user the rendered PDB/file, visualization mode, colors, and that the source snapshot is ProteinView FullHD 1920×1080.

Use `--mode cartoon|backbone|wireframe` and `--color structure|element|chain|plddt|bfactor|rainbow` only when the user asks for a specific presentation. With no mode, ProteinView uses cartoon for ordinary structures and automatically falls back to backbone for very large structures; the default color is structure.

For biologically focused requests:

- “color by chain” → add `--color chain`.
- “color chain A residue 42 red” → add `--residue-color A:42=FF0000`.
- “color insertion residue 42A in chain A cyan” → add `--residue-color 'A:42[A]=00FFFF'`.
- “show the interface for chain A” → add `--interface-chain A`.
- “show interface contacts/interactions for chain A” → add `--interface-chain A --interactions`.
- “hide ligands” → add `--hide-ligands`; ligands and ions are shown by default.

Repeat `--residue-color CHAIN:RES[ICODE]=RRGGBB` to select and color multiple exact polymer residues. Use six-digit RGB hex colors. The helper normalizes lowercase hex and insertion codes to uppercase, rejects duplicates, and records the normalized selections in its JSON result and live-view metadata. A selector without `[ICODE]` means the residue with a blank insertion code, not every residue sharing that sequence number. Exact residue colors take precedence over the global or interface palette for those polymer residues, so they can be combined with `--color`, `--interface-chain`, and interface interactions.

Interface mode highlights focus-chain interface residues in bright green and partner-chain interface residues in bright orange. With `--interactions`, cyan lines denote candidate hydrogen bonds, red salt bridges, yellow hydrophobic contacts, and gray other contacts, using ProteinView's 4.5 Å interface analysis.
The interface palette replaces regular color schemes, so never combine `--interface-chain` with `--color`.

Never launch ProteinView's standalone interactive TUI from inside Codex. The live pane uses ProteinView's bounded headless control server so Codex retains ownership of terminal input and layout. If the helper reports a missing or old binary, surface its actionable `PROTEINVIEW_BIN` message. If live or inline graphics are unavailable, return the saved PNG path as the fallback rather than substituting the HD text renderer.

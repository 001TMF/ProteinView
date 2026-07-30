# ProteinView for Codex

This plugin maps natural-language PDB requests to a cached ProteinView FullHD
1920×1080 render and displays the PNG through Codex's built-in image tool. It
also places private launch metadata beside the image. A live-viewer-enabled
Codex build uses that metadata to keep a controllable ProteinView pane beside
the chat; unmodified Codex versions simply ignore it and retain the PNG.

The ProteinView binary must support `--snapshot` and `--panel-server`. Put it on `PATH` or set
`PROTEINVIEW_BIN` to its absolute path. Compatible Kitty/Ghostty terminals show
the pixel image inline; F8 focuses or releases the live pane. Other terminals
retain the saved PNG path as a fallback.
Snapshots use a private, bounded temporary cache by default. Set
`PROTEINVIEW_CACHE_DIR` to a persistent writable cache root when previews must
survive operating-system temp cleanup. The helper creates and exclusively
prunes its own `proteinview/fullhd-v1` child; it never prunes the supplied root.

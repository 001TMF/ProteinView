use std::sync::mpsc;

use ratatui_image::picker::Picker;

use crate::model::interface::{InterfaceAnalysis, analyze_binding_pockets, analyze_interface};
use crate::model::protein::{Protein, Residue};
use crate::model::residue_selection::ResidueSelection;
use crate::model::selection::ResidueColorOverrides;
use crate::model::sequence::{SeqRow, SequenceLayout, wrap_for_width};
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, ColorSchemeType};
use crate::render::ribbon::{RibbonTriangle, generate_ribbon_mesh};

/// Structures with more residues than this threshold trigger performance
/// optimizations (background interface analysis, backbone default, reduced LOD).
pub const LARGE_STRUCTURE_THRESHOLD: usize = 5000;

/// Upper bound on the FullHD framebuffer, in pixels.
///
/// A graphics-protocol viewport is sized in *device* pixels, so on a HiDPI panel
/// it is several times the area the cell grid suggests, and every per-pixel
/// stage scales with it.  This is a backstop against a framebuffer so large it
/// costs real memory, not a frame-rate control: a still frame is rendered once
/// and then the loop idles, and everything drawn *while* the view moves is
/// already quartered by [`FULLHD_INTERACTIVE_SCALE`].
///
/// Set it generously, because capping is not free.  The terminal scales the
/// result back up, and a cap that barely engages buys a few percent of the
/// pixels in exchange for a non-integer resample of every still frame — worse
/// output for no useful saving.  12 MP clears a 4K viewport and a full-screen
/// HiDPI laptop with room to spare, so the cap only meets 5K and above, where
/// the framebuffer would otherwise run past a hundred megabytes.
pub const FULLHD_MAX_PIXELS: f64 = 12_000_000.0;

/// Resolution multiplier used while the camera is moving.
///
/// Halving each axis quarters every per-pixel cost, and the terminal scales the
/// result back up via the protocol's `c=`/`r=` keys.  Motion hides the
/// softness; the full-resolution frame lands as soon as the camera settles.
pub const FULLHD_INTERACTIVE_SCALE: f64 = 0.5;

/// How long after the last camera change the view still counts as interacting.
///
/// Long enough to cover the gap between key repeats, so held keys never
/// oscillate between resolutions, and short enough that the sharp frame feels
/// immediate once the user stops.
pub const INTERACTION_LINGER: std::time::Duration = std::time::Duration::from_millis(220);

/// Still-frame pixel dimensions of the FullHD framebuffer for a viewport of
/// `vp_cols` by `vp_rows` cells.
///
/// This is the single source of truth for FullHD sizing: both the zoom
/// calculation and the renderer go through it, so the framebuffer and the zoom
/// computed for it can never disagree.  While the camera is moving the renderer
/// scales this down by [`FULLHD_INTERACTIVE_SCALE`]; the still-frame size is
/// what zoom is defined against.
pub fn fullhd_framebuffer_size(
    vp_cols: f64,
    vp_rows: f64,
    font_w: u16,
    font_h: u16,
    is_graphics: bool,
) -> (f64, f64) {
    if !is_graphics {
        // Colored-braille fallback: 2x4 dots per cell.
        return (vp_cols * 2.0, vp_rows * 4.0);
    }

    let native_w = vp_cols * f64::from(font_w);
    let native_h = vp_rows * f64::from(font_h);

    // Cap by area rather than by either axis, so ultrawide and tall windows are
    // treated alike.
    let area = native_w * native_h;
    let scale = if area > FULLHD_MAX_PIXELS {
        (FULLHD_MAX_PIXELS / area).sqrt()
    } else {
        1.0
    };

    ((native_w * scale).max(1.0), (native_h * scale).max(1.0))
}

/// Visualization mode
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VizMode {
    Backbone,
    Cartoon,
    Wireframe,
}

impl VizMode {
    pub fn next(&self) -> Self {
        match self {
            Self::Backbone => Self::Cartoon,
            Self::Cartoon => Self::Wireframe,
            Self::Wireframe => Self::Backbone,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Self::Backbone => "Backbone",
            Self::Cartoon => "Cartoon",
            Self::Wireframe => "Wireframe",
        }
    }
}

/// Rendering mode for the 3D viewport
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RenderMode {
    /// Braille dots - highest text-mode spatial resolution, monochrome per cell
    Braille,
    /// HD-quality colored braille via software rasterizer (Lambert shading,
    /// z-buffer, depth fog).  Fast everywhere including SSH.
    HalfBlock,
    /// Same 2x4 braille grid as [`RenderMode::HalfBlock`], but rasterized into a
    /// supersampled framebuffer and box-filtered back down, with color
    /// quantization applied during the conversion.  Anti-aliased silhouettes and
    /// stable per-cell color, for the same number of characters on the wire.
    HalfBlockPlus,
    /// Full pixel graphics via Sixel/Kitty/iTerm2 - best quality, high bandwidth
    FullHD,
}

impl RenderMode {
    pub fn name(&self) -> &str {
        match self {
            Self::Braille => "Braille",
            Self::HalfBlock => "HD",
            Self::HalfBlockPlus => "HDplus",
            Self::FullHD => "FullHD",
        }
    }
}

/// Whether the terminal session is local or over SSH.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConnectionType {
    Local,
    Ssh,
}

impl ConnectionType {
    /// Detect whether the current session is running over SSH.
    ///
    /// This checks the `SSH_CLIENT`, `SSH_TTY`, and `SSH_CONNECTION`
    /// environment variables. Note that this can produce false positives
    /// in containers, CI environments, or VS Code Remote sessions where
    /// these variables may be inherited. Users can override the default
    /// render mode with `--fullhd` if detection is wrong.
    pub fn detect() -> Self {
        if std::env::var("SSH_CLIENT").is_ok()
            || std::env::var("SSH_TTY").is_ok()
            || std::env::var("SSH_CONNECTION").is_ok()
        {
            Self::Ssh
        } else {
            Self::Local
        }
    }
}

/// Configuration bundle for [`App::new`], replacing individual parameters
/// to avoid too_many_arguments.
pub struct AppConfig {
    pub render_mode: RenderMode,
    pub viz_mode: VizMode,
    pub user_explicit_mode: bool,
    pub color_override: Option<ColorSchemeType>,
    pub residue_colors: ResidueColorOverrides,
}

/// Main application state
pub struct App {
    pub protein: Protein,
    pub camera: Camera,
    pub color_scheme: ColorScheme,
    pub viz_mode: VizMode,
    pub current_chain: usize,
    pub render_mode: RenderMode,
    pub show_help: bool,
    pub show_ligands: bool,
    pub show_interface: bool,
    pub show_interactions: bool,
    pub interface_analysis: InterfaceAnalysis,
    pub should_quit: bool,
    /// Whether the B-factor column likely contains pLDDT confidence scores.
    pub has_plddt: bool,
    /// Cached ribbon mesh — regenerated only when color scheme changes.
    pub mesh_cache: Vec<RibbonTriangle>,
    mesh_dirty: bool,
    /// ratatui-image protocol picker for Sixel/Kitty/iTerm2 graphics.
    pub picker: Picker,
    /// Detected connection type (local vs SSH).
    pub connection_type: ConnectionType,
    /// Temporary warning when user enters FullHD over SSH.
    pub ssh_hd_warning: bool,
    /// Countdown frames to auto-dismiss the SSH HD warning (~90 frames = 3 seconds at 30fps).
    pub ssh_hd_warning_frames: u8,
    /// Set to `true` after a render-mode switch so the main loop can call
    /// `terminal.clear()` before the next draw, forcing ratatui to redraw
    /// every cell and preventing stale content from the previous mode.
    pub needs_clear: bool,
    /// Saved color scheme type to restore when leaving interface mode.
    /// When interface mode is active, we display Interface colors but
    /// preserve the user's chosen scheme so it can be restored on exit.
    saved_color_scheme_type: ColorSchemeType,
    residue_colors: ResidueColorOverrides,
    /// Whether interface analysis has been computed. For large structures
    /// (> LARGE_STRUCTURE_THRESHOLD residues), computation starts on a
    /// background thread at startup and completes before the user needs it.
    /// If the user requests interface mode before computation completes,
    /// the toggle is a no-op until the next frame.
    interface_computed: bool,
    /// Receiver for background interface analysis (large structures only).
    interface_rx: Option<mpsc::Receiver<InterfaceAnalysis>>,
    /// When the user last moved the camera, driving `is_interacting`.
    /// `None` until the first camera change.
    last_camera_change: Option<std::time::Instant>,
    /// Whether the terminal accepted a shared-memory graphics transmission at
    /// startup.  Set once by the probe in `main`; `false` over SSH.
    pub kitty_shm: bool,
    /// Next shared-memory slot to hand the terminal.  A `Cell` because the
    /// viewport renders from `&App` but each frame needs its own object.
    shm_slot: std::cell::Cell<u32>,
    /// Whether the scrollable chain-sequence panel is open.
    pub show_sequence: bool,
    /// Residues picked in the sequence panel.
    pub selection: ResidueSelection,
    /// Whether the selection is drawn as ball-and-stick in the 3D view.
    /// When off, selected residues are still marked with a highlight sphere.
    pub show_ball_stick: bool,
    /// Cursor position in the sequence panel, as `(chain, residue)` indices.
    pub seq_cursor: (usize, usize),
    /// Anchor for shift-extended range selection, cleared by any unshifted move.
    seq_anchor: Option<(usize, usize)>,
    /// First layout row drawn in the panel.
    pub seq_scroll: usize,
    /// Wrapped row layout, rebuilt whenever the panel width changes.
    seq_layout: SequenceLayout,
    /// Sequence rows the panel can currently show, set from the drawn area.
    seq_visible_rows: usize,
    /// Panel height in terminal rows, adjustable with `<` and `>`.
    pub seq_panel_height: u16,
    /// A chain to scroll to once the layout for the current width exists.
    /// Opening the panel happens before the frame that sizes it, so the jump
    /// has to wait for a layout to jump within.
    seq_pending_goto: Option<usize>,
}

impl App {
    pub fn new(
        mut protein: Protein,
        config: AppConfig,
        term_cols: u16,
        term_rows: u16,
        picker: Picker,
    ) -> Self {
        let AppConfig {
            render_mode,
            viz_mode,
            user_explicit_mode,
            color_override,
            residue_colors,
        } = config;
        protein.center();
        // If user explicitly requested pLDDT via CLI, trust that even if
        // the heuristic disagrees.
        let has_plddt = protein.has_plddt() || color_override == Some(ColorSchemeType::Plddt);
        let total_residues = protein.residue_count();
        let radius = protein.bounding_radius().max(1.0);

        let vp_rows = term_rows.saturating_sub(4) as f64;
        let vp_cols = term_cols as f64;
        let (font_w, font_h) = picker.font_size();

        let auto_zoom = match render_mode {
            RenderMode::FullHD => {
                let proto = picker.protocol_type();
                let (px_w, px_h) = if proto != ratatui_image::picker::ProtocolType::Halfblocks
                    && font_w > 0
                    && font_h > 0
                {
                    (vp_cols * font_w as f64, vp_rows * font_h as f64)
                } else {
                    // Fallback to braille-like resolution
                    (vp_cols * 2.0, vp_rows * 4.0)
                };
                0.9 * px_w.min(px_h) / (2.0 * radius)
            }
            RenderMode::HalfBlock | RenderMode::HalfBlockPlus => {
                let px_w = vp_cols * 2.0;
                let px_h = vp_rows * 4.0;
                0.9 * px_w.min(px_h) / (2.0 * radius)
            }
            RenderMode::Braille => {
                let px_w = vp_cols * 2.0;
                let px_h = vp_rows * 4.0;
                0.9 * px_w.min(px_h) / (2.0 * radius)
            }
        };
        let mut camera = Camera::default();
        camera.zoom = auto_zoom;

        let is_large = total_residues > LARGE_STRUCTURE_THRESHOLD;

        // For large structures, start interface analysis on a background thread
        // so it's ready by the time the user presses 'f'.
        let interface_rx = if is_large {
            let bg_protein = protein.clone();
            let (tx, rx) = mpsc::channel();
            std::thread::spawn(move || {
                let mut ia = analyze_interface(&bg_protein, 4.5);
                if !bg_protein.ligands.is_empty() {
                    ia.binding_pockets = Some(analyze_binding_pockets(&bg_protein, 4.5));
                }
                let _ = tx.send(ia);
            });
            // Interface analysis is running in the background — it'll be ready
            // by the time the user presses 'f'.
            Some(rx)
        } else {
            None
        };

        let (interface_analysis, interface_computed) = if is_large {
            let empty = InterfaceAnalysis {
                contacts: Vec::new(),
                interface_residues: std::collections::HashSet::new(),
                chain_interface_counts: vec![0; protein.chains.len()],
                total_interface_residues: 0,
                binding_pockets: None,
                interactions: Vec::new(),
            };
            (empty, false)
        } else {
            let mut ia = analyze_interface(&protein, 4.5);
            if !protein.ligands.is_empty() {
                ia.binding_pockets = Some(analyze_binding_pockets(&protein, 4.5));
            }
            (ia, true)
        };

        // For large structures, default to Backbone mode for instant
        // interactivity — but only if the user didn't explicitly choose a mode.
        let viz_mode = if is_large && !user_explicit_mode && viz_mode == VizMode::Cartoon {
            VizMode::Backbone
        } else {
            viz_mode
        };

        let initial_scheme = color_override.unwrap_or(ColorSchemeType::Structure);
        let color_scheme = ColorScheme::new(initial_scheme, total_residues)
            .with_residue_colors(residue_colors.clone());
        // Only build ribbon mesh eagerly if we're actually in Cartoon mode.
        // For Backbone/Wireframe, defer until the user switches to Cartoon.
        let (mesh_cache, mesh_dirty) = if viz_mode == VizMode::Cartoon {
            (generate_ribbon_mesh(&protein, &color_scheme), false)
        } else {
            (Vec::new(), true)
        };

        let connection_type = ConnectionType::detect();

        let selection = ResidueSelection::new(&protein);
        // The first chain with residues is where the sequence cursor starts;
        // a structure can open with an empty leading chain.
        let cursor_chain = protein
            .chains
            .iter()
            .position(|chain| !chain.residues.is_empty())
            .unwrap_or(0);

        Self {
            protein,
            camera,
            color_scheme,
            viz_mode,
            current_chain: 0,
            render_mode,
            show_help: false,
            show_ligands: true,
            show_interface: false,
            show_interactions: false,
            interface_analysis,
            should_quit: false,
            has_plddt,
            mesh_cache,
            mesh_dirty,
            picker,
            connection_type,
            ssh_hd_warning: false,
            ssh_hd_warning_frames: 0,
            needs_clear: false,
            saved_color_scheme_type: initial_scheme,
            residue_colors,
            interface_computed,
            interface_rx,
            last_camera_change: None,
            kitty_shm: false,
            shm_slot: std::cell::Cell::new(0),
            show_sequence: false,
            selection,
            show_ball_stick: true,
            seq_cursor: (cursor_chain, 0),
            seq_anchor: None,
            seq_scroll: 0,
            seq_layout: SequenceLayout::default(),
            seq_visible_rows: 1,
            seq_panel_height: DEFAULT_SEQUENCE_PANEL_HEIGHT,
            seq_pending_goto: None,
        }
    }

    pub fn cycle_color(&mut self) {
        if self.show_interface {
            // While interface mode is active, cycle the saved scheme so the
            // user's preference is tracked, but keep displaying Interface colors.
            self.saved_color_scheme_type = self.saved_color_scheme_type.next(self.has_plddt);
        } else {
            let next = self.color_scheme.scheme_type.next(self.has_plddt);
            self.color_scheme = ColorScheme::new(next, self.protein.residue_count())
                .with_residue_colors(self.residue_colors.clone());
            self.mesh_dirty = true;
        }
    }

    /// Whether the ribbon mesh will be rebuilt on the next `ribbon_mesh()` call.
    pub fn mesh_is_dirty(&self) -> bool {
        self.mesh_dirty
    }

    /// Whether a background interface analysis is still outstanding.
    pub fn interface_pending(&self) -> bool {
        !self.interface_computed
    }

    /// Poll the background interface analysis thread (non-blocking).
    /// Called each frame so results are absorbed as soon as they're ready.
    pub fn poll_background_interface(&mut self) {
        if self.interface_computed {
            return;
        }
        if let Some(rx) = &self.interface_rx {
            match rx.try_recv() {
                Ok(ia) => {
                    self.interface_analysis = ia;
                    self.interface_computed = true;
                    self.interface_rx = None;
                }
                Err(mpsc::TryRecvError::Empty) => {
                    // Still computing — nothing to do yet.
                }
                Err(mpsc::TryRecvError::Disconnected) => {
                    // Background thread panicked or dropped the sender.
                    // Drop the rx and fall back to synchronous computation.
                    self.interface_rx = None;
                    let mut ia = analyze_interface(&self.protein, 4.5);
                    if !self.protein.ligands.is_empty() {
                        ia.binding_pockets = Some(analyze_binding_pockets(&self.protein, 4.5));
                    }
                    self.interface_analysis = ia;
                    self.interface_computed = true;
                }
            }
        }
    }

    pub fn cycle_viz_mode(&mut self) {
        self.viz_mode = self.viz_mode.next();
    }

    fn rebuild_interface_colors(&mut self) {
        self.color_scheme = ColorScheme::new_interface(
            self.protein.residue_count(),
            self.current_chain,
            &self.interface_analysis,
            &self.protein,
        )
        .with_residue_colors(self.residue_colors.clone());
        self.mesh_dirty = true;
    }

    pub fn toggle_interface(&mut self) {
        self.show_interface = !self.show_interface;
        if self.show_interface {
            // Check if background analysis is ready, otherwise compute synchronously.
            if !self.interface_computed {
                // Determine background thread status without holding a
                // long-lived borrow on self.interface_rx.
                let bg_status = self.interface_rx.as_ref().map(|rx| rx.try_recv());
                match bg_status {
                    Some(Ok(ia)) => {
                        self.interface_analysis = ia;
                        self.interface_computed = true;
                        self.interface_rx = None;
                    }
                    Some(Err(mpsc::TryRecvError::Empty)) => {
                        // Still computing — don't enter interface mode yet.
                        // poll_background_interface() will absorb the result
                        // when ready; the user can press `f` again.
                        self.show_interface = false;
                        return;
                    }
                    Some(Err(mpsc::TryRecvError::Disconnected)) => {
                        // Thread panicked — drop the rx and fall through to
                        // synchronous computation below.
                        self.interface_rx = None;
                    }
                    None => {
                        // No background thread was spawned.
                    }
                }
                // If we still don't have it (no rx, or disconnected), compute synchronously.
                if !self.interface_computed {
                    let mut ia = analyze_interface(&self.protein, 4.5);
                    if !self.protein.ligands.is_empty() {
                        ia.binding_pockets = Some(analyze_binding_pockets(&self.protein, 4.5));
                    }
                    self.interface_analysis = ia;
                    self.interface_computed = true;
                }
            }
            // Save the user's current color scheme before switching to Interface
            self.saved_color_scheme_type = self.color_scheme.scheme_type;
            self.rebuild_interface_colors();
        } else {
            self.show_interactions = false;
            // Restore the user's saved color scheme instead of hardcoding Structure
            self.color_scheme =
                ColorScheme::new(self.saved_color_scheme_type, self.protein.residue_count())
                    .with_residue_colors(self.residue_colors.clone());
            self.mesh_dirty = true;
        }
    }

    pub fn toggle_interactions(&mut self) {
        if self.show_interface {
            self.show_interactions = !self.show_interactions;
        }
    }

    pub fn toggle_ligands(&mut self) {
        self.show_ligands = !self.show_ligands;
    }

    /// Get the cached ribbon mesh, regenerating if dirty.
    pub fn ribbon_mesh(&mut self) -> &[RibbonTriangle] {
        if self.mesh_dirty {
            self.mesh_cache = generate_ribbon_mesh(&self.protein, &self.color_scheme);
            self.mesh_dirty = false;
        }
        &self.mesh_cache
    }

    pub fn next_chain(&mut self) {
        if !self.protein.chains.is_empty() {
            self.current_chain = (self.current_chain + 1) % self.protein.chains.len();
            if self.show_interface {
                self.rebuild_interface_colors();
            }
            if self.show_sequence {
                self.seq_goto_chain(self.current_chain);
            }
        }
    }

    pub fn prev_chain(&mut self) {
        if !self.protein.chains.is_empty() {
            self.current_chain = if self.current_chain == 0 {
                self.protein.chains.len() - 1
            } else {
                self.current_chain - 1
            };
            if self.show_interface {
                self.rebuild_interface_colors();
            }
            if self.show_sequence {
                self.seq_goto_chain(self.current_chain);
            }
        }
    }

    pub fn chain_names(&self) -> Vec<String> {
        self.protein.chains.iter().map(|c| c.id.clone()).collect()
    }

    /// Returns `true` while the view is moving — auto-rotating, or within
    /// [`INTERACTION_LINGER`] of the last manual camera change.
    ///
    /// FullHD renders at reduced resolution whenever this holds, so rotating,
    /// panning and zooming stay smooth; the full-resolution frame is drawn once
    /// the camera settles.
    pub fn is_interacting(&self) -> bool {
        if self.camera.auto_rotate {
            return true;
        }
        self.last_camera_change
            .is_some_and(|at| at.elapsed() < INTERACTION_LINGER)
    }

    /// Claim the next shared-memory slot, cycling through the ring so the
    /// terminal is never still reading the object we are about to replace.
    pub fn next_shm_slot(&self) -> u32 {
        let slot = self.shm_slot.get();
        self.shm_slot
            .set((slot + 1) % crate::render::kitty_shm::SLOTS);
        slot
    }

    /// Record that the user just moved the camera.
    ///
    /// Called by the input loop for every key that rotates, pans, zooms or
    /// resets the view.
    pub fn note_camera_change(&mut self) {
        self.last_camera_change = Some(std::time::Instant::now());
    }

    pub fn tick(&mut self) {
        self.camera.tick();

        // Tick down SSH HD warning
        if self.ssh_hd_warning && self.ssh_hd_warning_frames > 0 {
            self.ssh_hd_warning_frames -= 1;
            if self.ssh_hd_warning_frames == 0 {
                self.ssh_hd_warning = false;
            }
        }
    }

    /// Mark the ribbon mesh cache as dirty, forcing a rebuild on the next frame.
    /// Called when terminal resize occurs or other events invalidate the mesh.
    pub fn mesh_dirty_flag(&mut self) {
        self.mesh_dirty = true;
    }

    /// Recalculate the zoom factor based on current render mode and terminal size.
    /// Call this after changing `render_mode` so the protein fills the viewport
    /// correctly for the new framebuffer dimensions.
    pub fn recalculate_zoom(&mut self, term_cols: u16, term_rows: u16) {
        let radius = self.protein.bounding_radius().max(1.0);
        let vp_rows = term_rows.saturating_sub(4) as f64;
        let vp_cols = term_cols as f64;
        let (font_w, font_h) = self.picker.font_size();

        let (px_w, px_h) = match self.render_mode {
            RenderMode::FullHD => {
                let proto = self.picker.protocol_type();
                let is_graphics = proto != ratatui_image::picker::ProtocolType::Halfblocks
                    && font_w > 0
                    && font_h > 0;
                // Zoom is defined against the still-frame resolution; the
                // renderer scales the camera itself when it drops to the
                // interactive resolution.
                fullhd_framebuffer_size(vp_cols, vp_rows, font_w, font_h, is_graphics)
            }
            RenderMode::HalfBlock | RenderMode::HalfBlockPlus => (vp_cols * 2.0, vp_rows * 4.0),
            RenderMode::Braille => (vp_cols * 2.0, vp_rows * 4.0),
        };
        self.camera.zoom = 0.9 * px_w.min(px_h) / (2.0 * radius);
    }

    /// Cycle lower render tiers: Braille -> HD -> HDplus -> Braille.
    /// From FullHD, steps down to HD (next lower tier).
    /// Bound to `m`.
    pub fn toggle_hd(&mut self, term_cols: u16, term_rows: u16) {
        self.render_mode = match self.render_mode {
            RenderMode::Braille => RenderMode::HalfBlock,
            RenderMode::HalfBlock => RenderMode::HalfBlockPlus,
            RenderMode::HalfBlockPlus => RenderMode::Braille,
            RenderMode::FullHD => RenderMode::HalfBlock,
        };
        // Dismiss any stale SSH warning (no longer in FullHD)
        self.ssh_hd_warning = false;
        self.ssh_hd_warning_frames = 0;
        self.needs_clear = true;
        self.recalculate_zoom(term_cols, term_rows);
    }

    /// Upgrade to FullHD (Sixel/Kitty) or back to HalfBlock.
    /// Bound to `M` (Shift+M).  Warns when entering FullHD over SSH.
    pub fn toggle_fullhd(&mut self, term_cols: u16, term_rows: u16) {
        self.render_mode = match self.render_mode {
            RenderMode::FullHD => RenderMode::HalfBlock,
            _ => RenderMode::FullHD,
        };

        self.needs_clear = true;

        if self.render_mode == RenderMode::FullHD && self.connection_type == ConnectionType::Ssh {
            self.ssh_hd_warning = true;
            self.ssh_hd_warning_frames = 90;
        } else {
            // Leaving FullHD — dismiss any active SSH warning
            self.ssh_hd_warning = false;
            self.ssh_hd_warning_frames = 0;
        }

        self.recalculate_zoom(term_cols, term_rows);
    }
}

/// Panel rows that are not sequence: the top border and the cursor line.
pub const SEQUENCE_PANEL_CHROME: u16 = 2;

/// Default height of the sequence panel, in terminal rows.
pub const DEFAULT_SEQUENCE_PANEL_HEIGHT: u16 = 10;

/// Bounds for `<` / `>` resizing.  The panel is additionally capped at half the
/// terminal height by the layout, so the 3D view never disappears.
pub const MIN_SEQUENCE_PANEL_HEIGHT: u16 = 4;
pub const MAX_SEQUENCE_PANEL_HEIGHT: u16 = 30;

/// Sequence panel: layout, cursor navigation and residue selection.
impl App {
    /// Open or close the panel.  The selection outlives it, so closing the
    /// panel keeps whatever is drawn in 3D.
    pub fn toggle_sequence_panel(&mut self) {
        self.show_sequence = !self.show_sequence;
        if self.show_sequence {
            // Follow the chain the rest of the UI is focused on, once there is
            // a layout to scroll within.
            self.seq_pending_goto = Some(self.current_chain);
        }
    }

    /// Tell the panel how much room it has, rebuilding the wrapped layout when
    /// the width changes.  Called once per frame before input is handled, so
    /// navigation and rendering always share one layout.
    pub fn set_sequence_viewport(&mut self, width: u16, height: u16) {
        let avail = width.saturating_sub(crate::ui::sequence_panel::GUTTER) as usize;
        let wrap = wrap_for_width(avail);
        if self.seq_layout.width != width || self.seq_layout.wrap != wrap {
            self.seq_layout = SequenceLayout::build(&self.protein, wrap, width);
        }
        self.seq_visible_rows = height.saturating_sub(SEQUENCE_PANEL_CHROME).max(1) as usize;
        if let Some(chain) = self.seq_pending_goto.take() {
            self.seq_goto_chain(chain);
        }
        self.scroll_cursor_into_view();
    }

    pub fn sequence_layout(&self) -> &SequenceLayout {
        &self.seq_layout
    }

    /// Grow or shrink the panel by `delta` rows, within the fixed bounds.
    pub fn resize_sequence_panel(&mut self, delta: i16) {
        let height = i32::from(self.seq_panel_height) + i32::from(delta);
        self.seq_panel_height = height.clamp(
            i32::from(MIN_SEQUENCE_PANEL_HEIGHT),
            i32::from(MAX_SEQUENCE_PANEL_HEIGHT),
        ) as u16;
    }

    /// The residue under the cursor, if the structure has one.
    pub fn seq_cursor_residue(&self) -> Option<(&crate::model::protein::Chain, &Residue)> {
        let chain = self.protein.chains.get(self.seq_cursor.0)?;
        let residue = chain.residues.get(self.seq_cursor.1)?;
        Some((chain, residue))
    }

    fn seq_chain_len(&self, chain: usize) -> usize {
        self.protein
            .chains
            .get(chain)
            .map_or(0, |chain| chain.residues.len())
    }

    /// Move the cursor to an exact residue, extending the selection when
    /// `extend` is set.
    ///
    /// Extension is anchored at the cursor position the first shifted move
    /// started from and only ever adds residues, so a shift-arrow can never
    /// silently drop part of an existing selection.
    fn seq_move_to(&mut self, chain: usize, residue: usize, extend: bool) {
        if self.seq_chain_len(chain) == 0 {
            return;
        }
        let previous = self.seq_cursor;
        let residue = residue.min(self.seq_chain_len(chain) - 1);
        self.seq_cursor = (chain, residue);
        self.current_chain = chain;
        // Interface coloring is keyed on the focus chain, so moving the cursor
        // into another chain has to refresh it exactly as `[` / `]` does.
        if previous.0 != chain && self.show_interface {
            self.rebuild_interface_colors();
        }

        if extend {
            let anchor = *self.seq_anchor.get_or_insert(previous);
            if anchor.0 == chain {
                self.selection.set_range(chain, anchor.1, residue, true);
            } else {
                // A range that crosses chains is not a range; restart the
                // anchor in the new chain rather than selecting everything in
                // between.
                self.selection.set(chain, residue, true);
                self.seq_anchor = Some((chain, residue));
            }
        } else {
            self.seq_anchor = None;
        }

        self.scroll_cursor_into_view();
    }

    /// Flat index of a residue across all chains, used for horizontal movement
    /// that runs off the end of a chain into the next one.
    fn seq_flat_index(&self, chain: usize, residue: usize) -> usize {
        self.protein.chains[..chain.min(self.protein.chains.len())]
            .iter()
            .map(|chain| chain.residues.len())
            .sum::<usize>()
            + residue
    }

    fn seq_from_flat(&self, mut index: usize) -> Option<(usize, usize)> {
        for (chain_index, chain) in self.protein.chains.iter().enumerate() {
            if index < chain.residues.len() {
                return Some((chain_index, index));
            }
            index -= chain.residues.len();
        }
        None
    }

    /// Move by `delta` residues, crossing chain boundaries at the ends.
    pub fn seq_move_horizontal(&mut self, delta: isize, extend: bool) {
        let total: usize = self
            .protein
            .chains
            .iter()
            .map(|chain| chain.residues.len())
            .sum();
        if total == 0 {
            return;
        }
        let flat = self.seq_flat_index(self.seq_cursor.0, self.seq_cursor.1) as isize;
        let target = (flat + delta).clamp(0, total as isize - 1) as usize;
        if let Some((chain, residue)) = self.seq_from_flat(target) {
            self.seq_move_to(chain, residue, extend);
        }
    }

    /// Move `delta` layout rows, keeping the column and skipping chain headers.
    pub fn seq_move_vertical(&mut self, delta: isize, extend: bool) {
        let Some((row, column)) = self.seq_layout.locate(self.seq_cursor.0, self.seq_cursor.1)
        else {
            return;
        };
        let rows = self.seq_layout.row_count() as isize;
        if rows == 0 {
            return;
        }
        let target = (row as isize + delta).clamp(0, rows - 1);
        let step = if delta >= 0 { 1 } else { -1 };

        // Headers carry no residue, so walk past them — first onwards in the
        // direction of travel, then backwards if that ran off the end.
        let landing = seek_residue_row(&self.seq_layout, target, step, rows)
            .or_else(|| seek_residue_row(&self.seq_layout, target, -step, rows));
        let Some(landing) = landing else {
            return;
        };
        if let Some((chain, residue)) = self.seq_layout.residue_at(landing, column) {
            self.seq_move_to(chain, residue, extend);
        }
    }

    /// Jump to the first or last residue of the cursor's chain.
    pub fn seq_move_to_chain_edge(&mut self, end: bool, extend: bool) {
        let chain = self.seq_cursor.0;
        let len = self.seq_chain_len(chain);
        if len == 0 {
            return;
        }
        self.seq_move_to(chain, if end { len - 1 } else { 0 }, extend);
    }

    /// Move by one screenful.
    pub fn seq_page(&mut self, forward: bool, extend: bool) {
        let rows = self.seq_visible_rows.max(1) as isize;
        self.seq_move_vertical(if forward { rows } else { -rows }, extend);
    }

    /// Put the cursor on the first residue of `chain` and show its header.
    pub fn seq_goto_chain(&mut self, chain: usize) {
        if self.protein.chains.is_empty() {
            return;
        }
        let chain = chain.min(self.protein.chains.len() - 1);
        if self.seq_chain_len(chain) == 0 {
            // Nothing to put a cursor on; still scroll the header into view.
            if let Some(header) = self.seq_layout.header_row(chain) {
                self.seq_scroll = self.clamp_scroll(header);
            }
            self.current_chain = chain;
            return;
        }
        self.seq_move_to(chain, 0, false);
        if let Some(header) = self.seq_layout.header_row(chain) {
            self.seq_scroll = self.clamp_scroll(header);
        }
    }

    /// Toggle the residue under the cursor.
    pub fn seq_toggle_selection(&mut self) {
        let (chain, residue) = self.seq_cursor;
        if self.seq_chain_len(chain) == 0 {
            return;
        }
        self.selection.toggle(chain, residue);
        // A later shift-arrow extends from here.
        self.seq_anchor = Some((chain, residue));
    }

    /// Select the cursor's whole chain, or clear it if it is already fully in.
    pub fn seq_toggle_chain_selection(&mut self) {
        let chain = self.seq_cursor.0;
        let len = self.seq_chain_len(chain);
        if len == 0 {
            return;
        }
        let fully_selected = (0..len).all(|residue| self.selection.contains(chain, residue));
        self.selection.set_chain(chain, !fully_selected);
        self.seq_anchor = None;
    }

    pub fn clear_selection(&mut self) {
        self.selection.clear();
        self.seq_anchor = None;
    }

    pub fn toggle_ball_stick(&mut self) {
        self.show_ball_stick = !self.show_ball_stick;
    }

    /// Centre the view on the selection, or on the cursor residue when nothing
    /// is selected.  Zoom is left alone: finding the residue is the hard part,
    /// and the user still owns the magnification.
    pub fn focus_on_selection(&mut self) -> bool {
        let target = self.selection.centroid(&self.protein).or_else(|| {
            let (_, residue) = self.seq_cursor_residue()?;
            let atoms = &residue.atoms;
            if atoms.is_empty() {
                return None;
            }
            let n = atoms.len() as f64;
            Some([
                atoms.iter().map(|a| a.x).sum::<f64>() / n,
                atoms.iter().map(|a| a.y).sum::<f64>() / n,
                atoms.iter().map(|a| a.z).sum::<f64>() / n,
            ])
        });
        let Some([x, y, z]) = target else {
            return false;
        };
        // `project` already includes the current pan, so subtracting the
        // projected offset lands the target exactly on the view centre for any
        // rotation or zoom.
        let projected = self.camera.project(x, y, z);
        self.camera.pan_x -= projected.x;
        self.camera.pan_y -= projected.y;
        self.note_camera_change();
        true
    }

    fn clamp_scroll(&self, scroll: usize) -> usize {
        let max = self
            .seq_layout
            .row_count()
            .saturating_sub(self.seq_visible_rows.max(1));
        scroll.min(max)
    }

    fn scroll_cursor_into_view(&mut self) {
        let Some((row, _)) = self.seq_layout.locate(self.seq_cursor.0, self.seq_cursor.1) else {
            return;
        };
        let visible = self.seq_visible_rows.max(1);
        if row < self.seq_scroll {
            self.seq_scroll = row;
        } else if row >= self.seq_scroll + visible {
            self.seq_scroll = row + 1 - visible;
        }
        self.seq_scroll = self.clamp_scroll(self.seq_scroll);
    }
}

/// First residue row at or after `from`, walking in `step` direction.
fn seek_residue_row(
    layout: &SequenceLayout,
    from: isize,
    step: isize,
    rows: isize,
) -> Option<usize> {
    let mut row = from;
    while (0..rows).contains(&row) {
        if !matches!(layout.rows[row as usize], SeqRow::Header(_)) {
            return Some(row as usize);
        }
        row += step;
    }
    None
}

#[cfg(test)]
mod sequence_navigation_tests {
    use super::*;
    use crate::model::protein::{Atom, Chain, MoleculeType, SecondaryStructure};
    use crate::model::selection::ResidueColorOverrides;
    use ratatui_image::picker::Picker;

    fn chain(id: &str, count: usize) -> Chain {
        Chain {
            id: id.to_string(),
            molecule_type: MoleculeType::Protein,
            residues: (0..count)
                .map(|i| Residue {
                    name: "ALA".to_string(),
                    seq_num: i as i32 + 1,
                    insertion_code: None,
                    atoms: vec![Atom {
                        name: "CA".to_string(),
                        element: "C".to_string(),
                        x: i as f64,
                        y: 0.0,
                        z: 0.0,
                        b_factor: 10.0,
                        is_backbone: true,
                        is_hetero: false,
                    }],
                    secondary_structure: SecondaryStructure::Coil,
                })
                .collect(),
        }
    }

    /// Chains of 25, 0 and 8 residues: the empty one in the middle is what
    /// makes navigation interesting.
    fn app() -> App {
        let protein = Protein {
            name: "nav".to_string(),
            chains: vec![chain("A", 25), chain("B", 0), chain("C", 8)],
            ligands: Vec::new(),
        };
        let mut app = App::new(
            protein,
            AppConfig {
                render_mode: RenderMode::Braille,
                viz_mode: VizMode::Backbone,
                user_explicit_mode: true,
                color_override: None,
                residue_colors: ResidueColorOverrides::default(),
            },
            80,
            40,
            Picker::halfblocks(),
        );
        app.show_sequence = true;
        // Gutter plus three groups of ten: ten residues per row.
        app.set_sequence_viewport(7 + 10, 8);
        app
    }

    #[test]
    fn horizontal_movement_crosses_into_the_next_non_empty_chain() {
        let mut app = app();
        app.seq_move_to_chain_edge(true, false);
        assert_eq!(app.seq_cursor, (0, 24));
        // Past the end of chain A, straight over the empty chain B.
        app.seq_move_horizontal(1, false);
        assert_eq!(app.seq_cursor, (2, 0));
        assert_eq!(app.current_chain, 2);
        // And back again.
        app.seq_move_horizontal(-1, false);
        assert_eq!(app.seq_cursor, (0, 24));
    }

    #[test]
    fn vertical_movement_skips_chain_headers() {
        let mut app = app();
        app.seq_move_to(0, 22, false);
        // Row below the last row of chain A is chain C's header; the cursor
        // must land on residues, not on it.
        app.seq_move_vertical(1, false);
        assert_eq!(app.seq_cursor.0, 2, "expected to land in chain C");
        app.seq_move_vertical(-1, false);
        assert_eq!(app.seq_cursor.0, 0, "expected to land back in chain A");
    }

    #[test]
    fn shift_extension_selects_a_range_and_restarts_across_chains() {
        let mut app = app();
        app.seq_move_to(0, 4, false);
        app.seq_toggle_selection();
        app.seq_move_horizontal(3, true);
        assert_eq!(app.selection.count(), 4);
        assert!(app.selection.contains(0, 7));

        // Extending into another chain must not select everything between.
        app.seq_move_to(2, 3, true);
        assert_eq!(app.selection.count(), 5);
        assert!(app.selection.contains(2, 3));
        assert!(!app.selection.contains(2, 0));
    }

    #[test]
    fn the_cursor_stays_inside_the_visible_rows() {
        let mut app = app();
        let visible = app.seq_visible_rows;
        app.seq_move_to_chain_edge(true, false);
        let (row, _) = app
            .sequence_layout()
            .locate(app.seq_cursor.0, app.seq_cursor.1)
            .unwrap();
        assert!(
            (app.seq_scroll..app.seq_scroll + visible).contains(&row),
            "row {row} outside scroll window {}..{}",
            app.seq_scroll,
            app.seq_scroll + visible
        );
    }

    #[test]
    fn empty_chains_never_take_the_cursor() {
        let mut app = app();
        app.seq_goto_chain(1);
        assert_eq!(app.current_chain, 1, "the header still gets focus");
        assert_ne!(app.seq_cursor.0, 1, "but no residue cursor lands there");
    }

    #[test]
    fn focusing_the_selection_centres_it() {
        let mut app = app();
        app.selection.set_range(0, 0, 4, true);
        assert!(app.focus_on_selection());
        let centroid = app.selection.centroid(&app.protein).unwrap();
        let projected = app.camera.project(centroid[0], centroid[1], centroid[2]);
        assert!(
            projected.x.abs() < 1e-9 && projected.y.abs() < 1e-9,
            "selection should project to the view centre, got {projected:?}"
        );
    }

    #[test]
    fn the_panel_never_squeezes_out_the_viewport() {
        let mut app = app();
        for _ in 0..100 {
            app.resize_sequence_panel(1);
        }
        assert_eq!(app.seq_panel_height, MAX_SEQUENCE_PANEL_HEIGHT);
        // On a 20-row terminal the layout still leaves the rest of the UI room.
        assert_eq!(
            crate::ui::sequence_panel::height_for(&app, 20),
            13,
            "panel should give back the seven rows the rest of the UI needs"
        );
        for _ in 0..100 {
            app.resize_sequence_panel(-1);
        }
        assert_eq!(app.seq_panel_height, MIN_SEQUENCE_PANEL_HEIGHT);
    }
}

#[cfg(test)]
mod fullhd_sizing_tests {
    use super::*;

    /// Measured: a full-screen kitty at font_size 14 on a 2560x1600 Retina
    /// panel reports 20x43 device-pixel cells over a 144x36 viewport.
    const RETINA: (f64, f64, u16, u16) = (144.0, 36.0, 20, 43);

    #[test]
    fn native_resolution_is_used_below_the_cap() {
        let (cols, rows, fw, fh) = RETINA;
        let (w, h) = fullhd_framebuffer_size(cols, rows, fw, fh, true);
        assert_eq!((w, h), (cols * f64::from(fw), rows * f64::from(fh)));
        assert!(
            w * h <= FULLHD_MAX_PIXELS,
            "{w}x{h} should be under the cap"
        );
    }

    /// Regression: the cap was first set at 4 MP, which a real full-screen
    /// HiDPI laptop viewport (4.46 MP) tripped by 10% — paying a non-integer
    /// upscale of every still frame to save almost nothing.
    #[test]
    fn a_full_screen_hidpi_laptop_is_not_capped() {
        let (cols, rows, fw, fh) = RETINA;
        let native = cols * f64::from(fw) * rows * f64::from(fh);
        assert!(
            native < FULLHD_MAX_PIXELS,
            "{native} px viewport should render natively, cap is {FULLHD_MAX_PIXELS}"
        );
    }

    /// A 4K viewport should also pass through untouched.
    #[test]
    fn a_4k_viewport_is_not_capped() {
        let (w, h) = fullhd_framebuffer_size(3840.0, 2160.0, 1, 1, true);
        assert_eq!((w, h), (3840.0, 2160.0));
    }

    #[test]
    fn oversized_viewports_are_capped_by_area_and_keep_their_aspect() {
        // A 6K display: well past the cap.
        let (native_w, native_h) = (6016.0, 3384.0);
        let (w, h) = fullhd_framebuffer_size(6016.0, 3384.0, 1, 1, true);

        assert!(
            w * h <= FULLHD_MAX_PIXELS * 1.001,
            "{w}x{h} = {} px exceeds the cap",
            w * h
        );
        let aspect_error = (w / h) - (native_w / native_h);
        assert!(
            aspect_error.abs() < 1e-6,
            "aspect drifted by {aspect_error}"
        );
    }

    #[test]
    fn braille_fallback_ignores_font_size() {
        let (w, h) = fullhd_framebuffer_size(100.0, 40.0, 17, 35, false);
        assert_eq!((w, h), (200.0, 160.0));
    }

    #[test]
    fn dimensions_never_collapse_to_zero() {
        let (w, h) = fullhd_framebuffer_size(0.0, 0.0, 17, 35, true);
        assert!(w >= 1.0 && h >= 1.0, "got {w}x{h}");
    }

    #[test]
    fn interactive_scale_quarters_the_pixel_count() {
        let (cols, rows, fw, fh) = RETINA;
        let (still_w, still_h) = fullhd_framebuffer_size(cols, rows, fw, fh, true);
        let (moving_w, moving_h) = (
            still_w * FULLHD_INTERACTIVE_SCALE,
            still_h * FULLHD_INTERACTIVE_SCALE,
        );
        let ratio = (still_w * still_h) / (moving_w * moving_h);
        assert!(
            (ratio - 4.0).abs() < 1e-9,
            "expected 4x fewer pixels, got {ratio}"
        );
    }
}

use std::sync::mpsc;

use ratatui_image::picker::Picker;

use crate::model::interface::{InterfaceAnalysis, analyze_binding_pockets, analyze_interface};
use crate::model::protein::Protein;
use crate::model::selection::ResidueColorOverrides;
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, ColorSchemeType};
use crate::render::ribbon::{RibbonTriangle, generate_ribbon_mesh};

/// Structures with more residues than this threshold trigger performance
/// optimizations (background interface analysis, backbone default, reduced LOD).
pub const LARGE_STRUCTURE_THRESHOLD: usize = 5000;

/// Upper bound on the FullHD framebuffer, in pixels.
///
/// A graphics-protocol viewport is sized in *device* pixels, so on a HiDPI
/// panel it is four times the area the cell grid suggests, and every per-pixel
/// stage scales with it.  This caps the still-frame resolution on very large or
/// very dense displays; below the cap the render stays at native resolution, so
/// a normal window is unaffected.  4 MP covers a full-screen Retina laptop.
pub const FULLHD_MAX_PIXELS: f64 = 4_000_000.0;

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

#[cfg(test)]
mod fullhd_sizing_tests {
    use super::*;

    /// A full-screen kitty at font_size 14 on a 2560x1600 Retina panel.
    const RETINA: (f64, f64, u16, u16) = (150.0, 41.0, 17, 35);

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

    #[test]
    fn oversized_viewports_are_capped_by_area_and_keep_their_aspect() {
        // A 5K display: well past the cap.
        let (native_w, native_h) = (5120.0, 2880.0);
        let (w, h) = fullhd_framebuffer_size(5120.0, 2880.0, 1, 1, true);

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

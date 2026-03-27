use ratatui_image::picker::Picker;

use crate::model::interface::{analyze_interface, analyze_binding_pockets, InterfaceAnalysis};
use crate::model::protein::Protein;
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, ColorSchemeType};
use crate::render::ribbon::{generate_ribbon_mesh, RibbonTriangle};

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
    /// Full pixel graphics via Sixel/Kitty/iTerm2 - best quality, high bandwidth
    FullHD,
}

impl RenderMode {
    pub fn name(&self) -> &str {
        match self {
            Self::Braille => "Braille",
            Self::HalfBlock => "HD",
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
}

impl App {
    pub fn new(
        mut protein: Protein,
        render_mode: RenderMode,
        term_cols: u16,
        term_rows: u16,
        picker: Picker,
        color_override: Option<ColorSchemeType>,
        viz_mode: VizMode,
    ) -> Self {
        protein.center();
        // If user explicitly requested pLDDT via CLI, trust that even if
        // the heuristic disagrees.
        let has_plddt = protein.has_plddt()
            || color_override == Some(ColorSchemeType::Plddt);
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
            RenderMode::HalfBlock => {
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

        // Pre-compute interface analysis (4.5A cutoff)
        let interface_analysis = {
            let mut ia = analyze_interface(&protein, 4.5);
            if !protein.ligands.is_empty() {
                ia.binding_pockets = Some(analyze_binding_pockets(&protein, 4.5));
            }
            ia
        };

        let initial_scheme = color_override.unwrap_or(ColorSchemeType::Structure);
        let color_scheme = ColorScheme::new(initial_scheme, total_residues);
        let mesh_cache = generate_ribbon_mesh(&protein, &color_scheme);

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
            mesh_dirty: false,
            picker,
            connection_type,
            ssh_hd_warning: false,
            ssh_hd_warning_frames: 0,
            needs_clear: false,
            saved_color_scheme_type: initial_scheme,
        }
    }

    pub fn cycle_color(&mut self) {
        if self.show_interface {
            // While interface mode is active, cycle the saved scheme so the
            // user's preference is tracked, but keep displaying Interface colors.
            self.saved_color_scheme_type = self.saved_color_scheme_type.next(self.has_plddt);
        } else {
            let next = self.color_scheme.scheme_type.next(self.has_plddt);
            self.color_scheme = ColorScheme::new(next, self.protein.residue_count());
            self.mesh_dirty = true;
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
        );
        self.mesh_dirty = true;
    }

    pub fn toggle_interface(&mut self) {
        self.show_interface = !self.show_interface;
        if self.show_interface {
            // Save the user's current color scheme before switching to Interface
            self.saved_color_scheme_type = self.color_scheme.scheme_type;
            self.rebuild_interface_colors();
        } else {
            self.show_interactions = false;
            // Restore the user's saved color scheme instead of hardcoding Structure
            self.color_scheme = ColorScheme::new(
                self.saved_color_scheme_type,
                self.protein.residue_count(),
            );
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
                if proto != ratatui_image::picker::ProtocolType::Halfblocks
                    && font_w > 0
                    && font_h > 0
                {
                    (vp_cols * font_w as f64, vp_rows * font_h as f64)
                } else {
                    (vp_cols * 2.0, vp_rows * 4.0)
                }
            }
            RenderMode::HalfBlock => (vp_cols * 2.0, vp_rows * 4.0),
            RenderMode::Braille => (vp_cols * 2.0, vp_rows * 4.0),
        };
        self.camera.zoom = 0.9 * px_w.min(px_h) / (2.0 * radius);
    }

    /// Toggle between Braille and HalfBlock (HD) mode.
    /// Bound to `m`.
    pub fn toggle_hd(&mut self, term_cols: u16, term_rows: u16) {
        self.render_mode = match self.render_mode {
            RenderMode::Braille => RenderMode::HalfBlock,
            _ => RenderMode::Braille,
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

        if self.render_mode == RenderMode::FullHD
            && self.connection_type == ConnectionType::Ssh
        {
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

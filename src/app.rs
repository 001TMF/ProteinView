use ratatui_image::picker::Picker;

use crate::model::interface::{InterfaceAnalysis, analyze_interface};
use crate::model::protein::Protein;
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, ColorSchemeType};

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

    pub fn from_cli(mode: &str) -> Self {
        match mode.to_ascii_lowercase().as_str() {
            "backbone" => Self::Backbone,
            "wireframe" => Self::Wireframe,
            "cartoon" => Self::Cartoon,
            _ => Self::Cartoon,
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
    pub hd_mode: bool,
    pub show_help: bool,
    pub show_interface: bool,
    pub interface_analysis: InterfaceAnalysis,
    pub should_quit: bool,
    /// ratatui-image protocol picker for Sixel/Kitty/iTerm2 graphics.
    pub picker: Picker,
}

impl App {
    pub fn new(
        mut protein: Protein,
        hd_mode: bool,
        viz_mode: VizMode,
        term_cols: u16,
        term_rows: u16,
        picker: Picker,
    ) -> Self {
        protein.center();
        let total_residues = protein.residue_count();
        let radius = protein.bounding_radius().max(1.0);
        // Dynamic zoom based on actual terminal size.
        // When a true graphics protocol (Sixel/Kitty/iTerm2) is available,
        // we render at full pixel resolution (cols*font_w x rows*font_h)
        // instead of the half-block resolution (cols x rows*2).  Scale the
        // zoom factor accordingly so the protein fills the viewport.
        let vp_rows = term_rows.saturating_sub(4) as f64;
        let vp_cols = term_cols as f64;
        let (font_w, font_h) = picker.font_size();
        let has_graphics = hd_mode
            && picker.protocol_type() != ratatui_image::picker::ProtocolType::Halfblocks
            && font_w > 0
            && font_h > 0;
        let auto_zoom = if hd_mode {
            let (px_w, px_h) = if has_graphics {
                (vp_cols * font_w as f64, vp_rows * font_h as f64)
            } else {
                // Must match viewport.rs braille dimensions: cols*2, rows*4
                (vp_cols * 2.0, vp_rows * 4.0)
            };
            0.9 * px_w.min(px_h) / (2.0 * radius)
        } else {
            0.9 * (vp_cols * 2.0).min(vp_rows * 4.0) / (2.0 * radius)
        };
        let mut camera = Camera::default();
        camera.zoom = auto_zoom;

        // Pre-compute interface analysis (4.5A cutoff)
        let interface_analysis = analyze_interface(&protein, 4.5);

        let color_scheme = ColorScheme::new(ColorSchemeType::Structure, total_residues);
        Self {
            protein,
            camera,
            color_scheme,
            viz_mode,
            current_chain: 0,
            hd_mode,
            show_help: false,
            show_interface: false,
            interface_analysis,
            should_quit: false,
            picker,
        }
    }

    pub fn cycle_color(&mut self) {
        let next = self.color_scheme.scheme_type.next();
        self.color_scheme = ColorScheme::new(next, self.protein.residue_count());
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
    }

    pub fn toggle_interface(&mut self) {
        self.show_interface = !self.show_interface;
        if self.show_interface {
            self.rebuild_interface_colors();
        } else {
            self.color_scheme =
                ColorScheme::new(ColorSchemeType::Structure, self.protein.residue_count());
        }
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
    }

    /// Recalculate the zoom factor based on current HD mode and terminal size.
    /// Call this after toggling `hd_mode` so the protein fills the viewport
    /// correctly for the new framebuffer dimensions.
    pub fn recalculate_zoom(&mut self, term_cols: u16, term_rows: u16) {
        let radius = self.protein.bounding_radius().max(1.0);
        let vp_rows = term_rows.saturating_sub(4) as f64;
        let vp_cols = term_cols as f64;
        let (font_w, font_h) = self.picker.font_size();
        let has_graphics = self.hd_mode
            && self.picker.protocol_type() != ratatui_image::picker::ProtocolType::Halfblocks
            && font_w > 0
            && font_h > 0;
        if self.hd_mode {
            let (px_w, px_h) = if has_graphics {
                (vp_cols * font_w as f64, vp_rows * font_h as f64)
            } else {
                // Must match viewport.rs braille dimensions: cols*2, rows*4
                (vp_cols * 2.0, vp_rows * 4.0)
            };
            self.camera.zoom = 0.9 * px_w.min(px_h) / (2.0 * radius);
        } else {
            self.camera.zoom = 0.9 * (vp_cols * 2.0).min(vp_rows * 4.0) / (2.0 * radius);
        }
    }
}

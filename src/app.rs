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
    fn auto_zoom_for_radius(
        &self,
        radius: f64,
        hd_mode: bool,
        term_cols: u16,
        term_rows: u16,
    ) -> f64 {
        let radius = radius.max(1.0);
        let vp_rows = term_rows.saturating_sub(4) as f64;
        let vp_cols = term_cols as f64;
        let (font_w, font_h) = self.picker.font_size();
        let has_graphics = hd_mode
            && self.picker.protocol_type() != ratatui_image::picker::ProtocolType::Halfblocks
            && font_w > 0
            && font_h > 0;

        if hd_mode {
            let (px_w, px_h) = if has_graphics {
                (vp_cols * font_w as f64, vp_rows * font_h as f64)
            } else {
                // Must match viewport.rs braille dimensions: cols*2, rows*4
                (vp_cols * 2.0, vp_rows * 4.0)
            };
            0.9 * px_w.min(px_h) / (2.0 * radius)
        } else {
            0.9 * (vp_cols * 2.0).min(vp_rows * 4.0) / (2.0 * radius)
        }
    }

    fn auto_zoom_for_mode(&self, hd_mode: bool, term_cols: u16, term_rows: u16) -> f64 {
        self.auto_zoom_for_radius(
            self.protein.bounding_radius(),
            hd_mode,
            term_cols,
            term_rows,
        )
    }

    pub fn new(
        mut protein: Protein,
        hd_mode: bool,
        initial_color: ColorSchemeType,
        viz_mode: VizMode,
        term_cols: u16,
        term_rows: u16,
        picker: Picker,
    ) -> Self {
        protein.center();
        let total_residues = protein.residue_count();
        // Pre-compute interface analysis (4.5A cutoff)
        let interface_analysis = analyze_interface(&protein, 4.5);

        let color_scheme = ColorScheme::new(initial_color, total_residues);
        let mut app = Self {
            protein,
            camera: Camera::default(),
            color_scheme,
            viz_mode,
            current_chain: 0,
            hd_mode,
            show_help: false,
            show_interface: false,
            interface_analysis,
            should_quit: false,
            picker,
        };
        app.camera.zoom = app.auto_zoom_for_mode(app.hd_mode, term_cols, term_rows);
        app
    }

    fn atomic_mass(element: &str) -> f64 {
        match element.trim().to_ascii_uppercase().as_str() {
            "H" => 1.008,
            "C" => 12.011,
            "N" => 14.007,
            "O" => 15.999,
            "P" => 30.974,
            "S" => 32.06,
            _ => 12.0,
        }
    }

    fn interface_pivot_for_chain(&self, chain_idx: usize) -> Option<[f64; 3]> {
        let mut mx = 0.0;
        let mut my = 0.0;
        let mut mz = 0.0;
        let mut m_total = 0.0;

        for &(ci, ri) in &self.interface_analysis.interface_residues {
            if ci != chain_idx {
                continue;
            }
            let residue = self.protein.chains.get(ci)?.residues.get(ri)?;
            for atom in &residue.atoms {
                let mass = Self::atomic_mass(&atom.element);
                mx += atom.x * mass;
                my += atom.y * mass;
                mz += atom.z * mass;
                m_total += mass;
            }
        }

        if m_total > 0.0 {
            Some([mx / m_total, my / m_total, mz / m_total])
        } else {
            None
        }
    }

    fn chain_center_of_mass(&self, chain_idx: usize) -> Option<[f64; 3]> {
        let chain = self.protein.chains.get(chain_idx)?;
        let mut mx = 0.0;
        let mut my = 0.0;
        let mut mz = 0.0;
        let mut m_total = 0.0;

        for residue in &chain.residues {
            for atom in &residue.atoms {
                let mass = Self::atomic_mass(&atom.element);
                mx += atom.x * mass;
                my += atom.y * mass;
                mz += atom.z * mass;
                m_total += mass;
            }
        }

        if m_total > 0.0 {
            Some([mx / m_total, my / m_total, mz / m_total])
        } else {
            None
        }
    }

    fn chain_bounding_radius_from_pivot(&self, chain_idx: usize, pivot: [f64; 3]) -> Option<f64> {
        let chain = self.protein.chains.get(chain_idx)?;

        let mut max_backbone = 0.0f64;
        let mut found_backbone = false;
        let mut max_any = 0.0f64;
        let mut found_any = false;

        for residue in &chain.residues {
            for atom in &residue.atoms {
                let dx = atom.x - pivot[0];
                let dy = atom.y - pivot[1];
                let dz = atom.z - pivot[2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                max_any = max_any.max(d);
                found_any = true;
                if atom.is_backbone {
                    max_backbone = max_backbone.max(d);
                    found_backbone = true;
                }
            }
        }

        if found_backbone {
            Some(max_backbone)
        } else if found_any {
            Some(max_any)
        } else {
            None
        }
    }

    fn refresh_interface_pivot(&mut self) {
        if !self.show_interface {
            self.camera.set_pivot([0.0, 0.0, 0.0]);
            return;
        }
        if let Some(pivot) = self.interface_pivot_for_chain(self.current_chain) {
            self.camera.set_pivot(pivot);
        } else {
            self.camera.set_pivot([0.0, 0.0, 0.0]);
        }
    }

    pub fn cycle_color(&mut self) {
        let next = self.color_scheme.scheme_type.next(self.protein.has_plddt());
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
        self.refresh_interface_pivot();
    }

    pub fn next_chain(&mut self) {
        if !self.protein.chains.is_empty() {
            self.current_chain = (self.current_chain + 1) % self.protein.chains.len();
            if self.show_interface {
                self.rebuild_interface_colors();
                self.refresh_interface_pivot();
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
                self.refresh_interface_pivot();
            }
        }
    }

    /// Global focus helper (bound to `/`): advance chain, highlight focused
    /// chain, and fit it to the viewport using its COM and radius.
    pub fn focus_next_chain(&mut self, term_cols: u16, term_rows: u16) {
        if self.protein.chains.is_empty() {
            return;
        }

        self.current_chain = (self.current_chain + 1) % self.protein.chains.len();
        self.color_scheme = ColorScheme::new_focus_chain(
            self.protein.residue_count(),
            self.current_chain,
            &self.protein,
        );

        let pivot = self
            .chain_center_of_mass(self.current_chain)
            .unwrap_or([0.0, 0.0, 0.0]);
        self.camera.set_pivot(pivot);
        self.camera.pan_x = 0.0;
        self.camera.pan_y = 0.0;

        if let Some(radius) = self.chain_bounding_radius_from_pivot(self.current_chain, pivot) {
            self.camera.zoom =
                self.auto_zoom_for_radius(radius, self.hd_mode, term_cols, term_rows);
        }
    }

    pub fn chain_names(&self) -> Vec<String> {
        self.protein.chains.iter().map(|c| c.id.clone()).collect()
    }

    pub fn tick(&mut self) {
        self.camera.tick();
    }

    pub fn reset_camera(&mut self) {
        self.camera.reset();
        self.refresh_interface_pivot();
    }

    /// Toggle HD mode while preserving the current framing.
    /// Keeps the user's current zoom and pan proportional to the framebuffer
    /// resolution used by each render backend.
    pub fn toggle_hd_mode_preserve_view(&mut self, term_cols: u16, term_rows: u16) {
        let old_auto_zoom = self.auto_zoom_for_mode(self.hd_mode, term_cols, term_rows);
        self.hd_mode = !self.hd_mode;
        let new_auto_zoom = self.auto_zoom_for_mode(self.hd_mode, term_cols, term_rows);

        if old_auto_zoom > f64::EPSILON && new_auto_zoom > f64::EPSILON {
            let scale = new_auto_zoom / old_auto_zoom;
            self.camera.zoom *= scale;
            self.camera.pan_x *= scale;
            self.camera.pan_y *= scale;
        }
    }
}

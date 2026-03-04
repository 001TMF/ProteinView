use crate::model::protein::Protein;
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, ColorSchemeType};

/// Visualization mode
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VizMode {
    Backbone,
    BallAndStick,
    Wireframe,
}

impl VizMode {
    pub fn next(&self) -> Self {
        match self {
            Self::Backbone => Self::BallAndStick,
            Self::BallAndStick => Self::Wireframe,
            Self::Wireframe => Self::Backbone,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Self::Backbone => "Backbone",
            Self::BallAndStick => "Ball+Stick",
            Self::Wireframe => "Wireframe",
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
    pub should_quit: bool,
}

impl App {
    pub fn new(mut protein: Protein, hd_mode: bool) -> Self {
        protein.center();
        let total_residues = protein.residue_count();
        Self {
            protein,
            camera: Camera::default(),
            color_scheme: ColorScheme::new(ColorSchemeType::Structure, total_residues),
            viz_mode: VizMode::Backbone,
            current_chain: 0,
            hd_mode,
            show_help: false,
            should_quit: false,
        }
    }

    pub fn cycle_color(&mut self) {
        let next = self.color_scheme.scheme_type.next();
        self.color_scheme = ColorScheme::new(next, self.protein.residue_count());
    }

    pub fn cycle_viz_mode(&mut self) {
        self.viz_mode = self.viz_mode.next();
    }

    pub fn next_chain(&mut self) {
        if !self.protein.chains.is_empty() {
            self.current_chain = (self.current_chain + 1) % self.protein.chains.len();
        }
    }

    pub fn prev_chain(&mut self) {
        if !self.protein.chains.is_empty() {
            self.current_chain = if self.current_chain == 0 {
                self.protein.chains.len() - 1
            } else {
                self.current_chain - 1
            };
        }
    }

    pub fn tick(&mut self) {
        self.camera.tick();
    }
}

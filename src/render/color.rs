use ratatui::style::Color;
use crate::model::protein::{Chain, Residue, SecondaryStructure};

/// Available color schemes
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ColorSchemeType {
    Structure,
    Chain,
    Element,
    BFactor,
    Rainbow,
}

impl ColorSchemeType {
    pub fn next(&self) -> Self {
        match self {
            Self::Structure => Self::Chain,
            Self::Chain => Self::Element,
            Self::Element => Self::BFactor,
            Self::BFactor => Self::Rainbow,
            Self::Rainbow => Self::Structure,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Self::Structure => "Structure",
            Self::Chain => "Chain",
            Self::Element => "Element",
            Self::BFactor => "B-Factor",
            Self::Rainbow => "Rainbow",
        }
    }
}

/// Color scheme for rendering
#[derive(Debug, Clone)]
pub struct ColorScheme {
    pub scheme_type: ColorSchemeType,
    total_residues: usize,
}

impl ColorScheme {
    pub fn new(scheme_type: ColorSchemeType, total_residues: usize) -> Self {
        Self { scheme_type, total_residues }
    }

    /// Get color for a residue based on current scheme
    pub fn residue_color(&self, residue: &Residue, chain: &Chain) -> Color {
        match self.scheme_type {
            ColorSchemeType::Structure => self.structure_color(residue),
            ColorSchemeType::Chain => self.chain_color(chain),
            ColorSchemeType::Element => Color::Rgb(144, 144, 144), // Default to carbon gray
            ColorSchemeType::BFactor => self.bfactor_color(residue),
            ColorSchemeType::Rainbow => self.rainbow_color(residue),
        }
    }

    fn structure_color(&self, residue: &Residue) -> Color {
        match residue.secondary_structure {
            SecondaryStructure::Helix => Color::Rgb(255, 0, 128),   // Magenta-red
            SecondaryStructure::Sheet => Color::Rgb(255, 200, 0),   // Gold
            SecondaryStructure::Turn  => Color::Rgb(96, 128, 255),  // Pale blue
            SecondaryStructure::Coil  => Color::Rgb(0, 204, 0),     // Green
        }
    }

    fn chain_color(&self, chain: &Chain) -> Color {
        let chain_colors = [
            Color::Rgb(0, 180, 255),   // Cyan
            Color::Rgb(255, 100, 0),   // Orange
            Color::Rgb(0, 220, 100),   // Green
            Color::Rgb(255, 50, 150),  // Pink
            Color::Rgb(180, 100, 255), // Purple
            Color::Rgb(255, 220, 0),   // Yellow
            Color::Rgb(0, 200, 200),   // Teal
            Color::Rgb(255, 150, 150), // Salmon
        ];
        let idx = chain.id.bytes().next().unwrap_or(b'A') as usize % chain_colors.len();
        chain_colors[idx]
    }

    fn bfactor_color(&self, residue: &Residue) -> Color {
        // Average B-factor of residue atoms, map to blue-red gradient
        let avg_b: f64 = if residue.atoms.is_empty() {
            0.0
        } else {
            residue.atoms.iter().map(|a| a.b_factor).sum::<f64>() / residue.atoms.len() as f64
        };
        // Normalize: typical B-factors range 5-80
        let t = ((avg_b - 5.0) / 75.0).clamp(0.0, 1.0);
        let r = (t * 255.0) as u8;
        let b = ((1.0 - t) * 255.0) as u8;
        Color::Rgb(r, 0, b)
    }

    fn rainbow_color(&self, residue: &Residue) -> Color {
        if self.total_residues == 0 { return Color::White; }
        let t = residue.seq_num as f64 / self.total_residues as f64;
        // HSV to RGB with hue 0-300 (blue to red)
        let hue = (1.0 - t) * 300.0;
        let (r, g, b) = hsv_to_rgb(hue, 1.0, 1.0);
        Color::Rgb(r, g, b)
    }
}

/// Convert HSV to RGB (h: 0-360, s: 0-1, v: 0-1)
fn hsv_to_rgb(h: f64, s: f64, v: f64) -> (u8, u8, u8) {
    let c = v * s;
    let x = c * (1.0 - ((h / 60.0) % 2.0 - 1.0).abs());
    let m = v - c;
    let (r, g, b) = match h as u32 {
        0..=59 => (c, x, 0.0),
        60..=119 => (x, c, 0.0),
        120..=179 => (0.0, c, x),
        180..=239 => (0.0, x, c),
        240..=299 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    (((r + m) * 255.0) as u8, ((g + m) * 255.0) as u8, ((b + m) * 255.0) as u8)
}

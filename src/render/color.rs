use std::collections::HashSet;

use crate::model::interface::InterfaceAnalysis;
use crate::model::protein::{Atom, Chain, Residue, SecondaryStructure};
use ratatui::style::Color;

/// Available color schemes
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ColorSchemeType {
    Structure,
    Plddt,
    Chain,
    Element,
    BFactor,
    Rainbow,
    Interface,
}

impl ColorSchemeType {
    pub fn next(&self, has_plddt: bool) -> Self {
        match self {
            Self::Structure => {
                if has_plddt {
                    Self::Plddt
                } else {
                    Self::Chain
                }
            }
            Self::Plddt => Self::Chain,
            Self::Chain => Self::Element,
            Self::Element => Self::BFactor,
            Self::BFactor => Self::Rainbow,
            Self::Rainbow => Self::Structure,
            // Interface is toggled separately, skip it in the cycle
            Self::Interface => Self::Structure,
        }
    }

    pub fn from_cli(name: &str, has_plddt: bool) -> Self {
        match name.to_ascii_lowercase().as_str() {
            "structure" => Self::Structure,
            "plddt" => {
                if has_plddt {
                    Self::Plddt
                } else {
                    Self::Structure
                }
            }
            "chain" => Self::Chain,
            "element" => Self::Element,
            "bfactor" | "b-factor" => Self::BFactor,
            "rainbow" => Self::Rainbow,
            "interface" => Self::Interface,
            _ => Self::Structure,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Self::Structure => "Structure",
            Self::Plddt => "pLDDT",
            Self::Chain => "Chain",
            Self::Element => "Element",
            Self::BFactor => "B-Factor",
            Self::Rainbow => "Rainbow",
            Self::Interface => "Interface",
        }
    }
}

/// Color scheme for rendering
#[derive(Debug, Clone)]
pub struct ColorScheme {
    pub scheme_type: ColorSchemeType,
    total_residues: usize,
    /// For Interface mode: chain ID of the "focus" (antibody) chain
    focus_chain_id: String,
    /// For Interface mode: set of (chain_id, seq_num) at the interface
    interface_residues_by_id: HashSet<(String, i32)>,
}

impl ColorScheme {
    pub fn new(scheme_type: ColorSchemeType, total_residues: usize) -> Self {
        Self {
            scheme_type,
            total_residues,
            focus_chain_id: String::new(),
            interface_residues_by_id: HashSet::new(),
        }
    }

    pub fn new_interface(
        total_residues: usize,
        focus_chain: usize,
        analysis: &InterfaceAnalysis,
        protein: &crate::model::protein::Protein,
    ) -> Self {
        let focus_chain_id = protein
            .chains
            .get(focus_chain)
            .map(|c| c.id.clone())
            .unwrap_or_default();
        Self {
            scheme_type: ColorSchemeType::Interface,
            total_residues,
            focus_chain_id,
            interface_residues_by_id: analysis.interface_residues_by_id_with_protein(protein),
        }
    }

    /// Get color for a residue based on current scheme
    pub fn residue_color(&self, residue: &Residue, chain: &Chain) -> Color {
        match self.scheme_type {
            ColorSchemeType::Structure => self.structure_color(residue),
            ColorSchemeType::Plddt => self.plddt_color(residue),
            ColorSchemeType::Chain => self.chain_color(chain),
            ColorSchemeType::Element => Color::Rgb(144, 144, 144),
            ColorSchemeType::BFactor => self.bfactor_color(residue),
            ColorSchemeType::Rainbow => self.rainbow_color(residue),
            ColorSchemeType::Interface => self.interface_color(residue, chain),
        }
    }

    /// Get color for an individual atom, respecting the current scheme.
    pub fn atom_color(&self, atom: &Atom, residue: &Residue, chain: &Chain) -> Color {
        match self.scheme_type {
            ColorSchemeType::Element => Self::element_color(atom),
            ColorSchemeType::Interface => {
                if self.is_interface_residue(residue, chain) {
                    match atom.element.trim() {
                        // Keep heteroatom cues visible on highlighted interface residues.
                        "N" => Color::Rgb(48, 80, 248),
                        "O" => Color::Rgb(255, 13, 13),
                        _ => self.interface_color(residue, chain),
                    }
                } else {
                    self.interface_color(residue, chain)
                }
            }
            _ => self.residue_color(residue, chain),
        }
    }

    /// Interface color scheme:
    /// - Focus chain (antibody): green tones
    ///   - Interface residues: bright green
    ///   - Non-interface: dim green
    /// - Other chains (antigen): orange tones
    ///   - Interface residues: bright orange
    ///   - Non-interface: dim gray-brown
    fn interface_color(&self, residue: &Residue, chain: &Chain) -> Color {
        let is_contact = self.is_interface_residue(residue, chain);
        let is_focus = chain.id == self.focus_chain_id;

        match (is_focus, is_contact) {
            (true, true) => Color::Rgb(0, 255, 100), // Bright green — antibody interface
            (true, false) => Color::Rgb(40, 100, 60), // Dim green — antibody non-interface
            (false, true) => Color::Rgb(255, 165, 0), // Bright orange — antigen interface
            (false, false) => Color::Rgb(100, 80, 60), // Dim brown — antigen non-interface
        }
    }

    fn is_interface_residue(&self, residue: &Residue, chain: &Chain) -> bool {
        self.interface_residues_by_id
            .contains(&(chain.id.clone(), residue.seq_num))
    }

    /// CPK-style element coloring
    fn element_color(atom: &Atom) -> Color {
        match atom.element.trim() {
            "C" => Color::Rgb(144, 144, 144),
            "N" => Color::Rgb(48, 80, 248),
            "O" => Color::Rgb(255, 13, 13),
            "S" => Color::Rgb(255, 255, 48),
            "H" => Color::Rgb(255, 255, 255),
            "P" => Color::Rgb(255, 128, 0),
            "FE" | "Fe" => Color::Rgb(224, 102, 51),
            _ => Color::Rgb(200, 200, 200),
        }
    }

    fn structure_color(&self, residue: &Residue) -> Color {
        // Nucleotide residues get base-type coloring
        if let Some(color) = nucleotide_base_color(&residue.name) {
            return color;
        }

        match residue.secondary_structure {
            SecondaryStructure::Helix => Color::Rgb(255, 0, 128),
            SecondaryStructure::Sheet => Color::Rgb(255, 200, 0),
            SecondaryStructure::Turn => Color::Rgb(96, 128, 255),
            SecondaryStructure::Coil => Color::Rgb(0, 204, 0),
        }
    }

    fn chain_color(&self, chain: &Chain) -> Color {
        let chain_colors = [
            Color::Rgb(0, 180, 255),
            Color::Rgb(255, 100, 0),
            Color::Rgb(0, 220, 100),
            Color::Rgb(255, 50, 150),
            Color::Rgb(180, 100, 255),
            Color::Rgb(255, 220, 0),
            Color::Rgb(0, 200, 200),
            Color::Rgb(255, 150, 150),
        ];
        let idx = chain.id.bytes().next().unwrap_or(b'A') as usize % chain_colors.len();
        chain_colors[idx]
    }

    fn bfactor_color(&self, residue: &Residue) -> Color {
        let avg_b: f64 = if residue.atoms.is_empty() {
            0.0
        } else {
            residue.atoms.iter().map(|a| a.b_factor).sum::<f64>() / residue.atoms.len() as f64
        };
        let t = ((avg_b - 5.0) / 75.0).clamp(0.0, 1.0);
        let r = (t * 255.0) as u8;
        let b = ((1.0 - t) * 255.0) as u8;
        Color::Rgb(r, 0, b)
    }

    fn plddt_color(&self, residue: &Residue) -> Color {
        let avg_plddt: f64 = if residue.atoms.is_empty() {
            0.0
        } else {
            residue.atoms.iter().map(|a| a.b_factor).sum::<f64>() / residue.atoms.len() as f64
        };

        if avg_plddt >= 90.0 {
            Color::Rgb(0, 83, 214)
        } else if avg_plddt >= 70.0 {
            Color::Rgb(101, 203, 243)
        } else if avg_plddt >= 50.0 {
            Color::Rgb(255, 219, 19)
        } else {
            Color::Rgb(255, 125, 69)
        }
    }

    fn rainbow_color(&self, residue: &Residue) -> Color {
        if self.total_residues == 0 {
            return Color::White;
        }
        let t = residue.seq_num as f64 / self.total_residues as f64;
        let hue = (1.0 - t) * 300.0;
        let (r, g, b) = hsv_to_rgb(hue, 1.0, 1.0);
        Color::Rgb(r, g, b)
    }
}

/// Returns a base-type color for nucleotide residues, or `None` for non-nucleotides.
fn nucleotide_base_color(name: &str) -> Option<Color> {
    match name {
        "A" | "DA" | "AMP" => Some(Color::Rgb(220, 60, 60)), // Adenine — red
        "U" | "UMP" => Some(Color::Rgb(60, 60, 220)),        // Uracil — blue
        "T" | "DT" => Some(Color::Rgb(60, 60, 220)),         // Thymine — blue
        "G" | "DG" | "GMP" => Some(Color::Rgb(60, 180, 60)), // Guanine — green
        "C" | "DC" | "CMP" => Some(Color::Rgb(220, 200, 40)), // Cytosine — yellow
        "I" | "DI" => Some(Color::Rgb(150, 100, 180)),       // Inosine — purple
        _ => None,
    }
}

/// Convert a ratatui `Color` to an `[u8; 3]` RGB triple.
///
/// Returns `[180, 180, 180]` (light gray) for non-RGB color variants.
pub fn color_to_rgb(color: Color) -> [u8; 3] {
    match color {
        Color::Rgb(r, g, b) => [r, g, b],
        _ => [180, 180, 180],
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
    (
        ((r + m) * 255.0) as u8,
        ((g + m) * 255.0) as u8,
        ((b + m) * 255.0) as u8,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::protein::{Atom, Residue, SecondaryStructure, is_nucleotide};

    /// Build a minimal residue for testing color assignment.
    fn make_residue(name: &str, ss: SecondaryStructure) -> Residue {
        Residue {
            name: name.to_string(),
            seq_num: 1,
            atoms: vec![],
            secondary_structure: ss,
        }
    }

    // ---- is_nucleotide helper ----

    #[test]
    fn is_nucleotide_rna_bases() {
        for name in &["A", "U", "G", "C"] {
            assert!(
                is_nucleotide(name),
                "{name} should be recognized as nucleotide"
            );
        }
    }

    #[test]
    fn is_nucleotide_dna_bases() {
        for name in &["DA", "DT", "DG", "DC"] {
            assert!(
                is_nucleotide(name),
                "{name} should be recognized as nucleotide"
            );
        }
    }

    #[test]
    fn is_nucleotide_modified_forms() {
        for name in &["AMP", "UMP", "GMP", "CMP"] {
            assert!(
                is_nucleotide(name),
                "{name} should be recognized as nucleotide"
            );
        }
    }

    #[test]
    fn is_nucleotide_rejects_amino_acids() {
        for name in &["ALA", "GLY", "CYS", "THR", "TRP", "LEU"] {
            assert!(
                !is_nucleotide(name),
                "{name} should NOT be recognized as nucleotide"
            );
        }
    }

    // ---- structure_color: nucleotide residues ----

    #[test]
    fn structure_color_adenine_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(220, 60, 60);

        for name in &["A", "DA", "AMP"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(
                scheme.structure_color(&r),
                expected,
                "Adenine variant {name}"
            );
        }
    }

    #[test]
    fn structure_color_uracil_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(60, 60, 220);

        for name in &["U", "UMP"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(
                scheme.structure_color(&r),
                expected,
                "Uracil variant {name}"
            );
        }
    }

    #[test]
    fn structure_color_thymine_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(60, 60, 220);

        for name in &["T", "DT"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(
                scheme.structure_color(&r),
                expected,
                "Thymine variant {name}"
            );
        }
    }

    #[test]
    fn structure_color_guanine_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(60, 180, 60);

        for name in &["G", "DG", "GMP"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(
                scheme.structure_color(&r),
                expected,
                "Guanine variant {name}"
            );
        }
    }

    #[test]
    fn structure_color_cytosine_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(220, 200, 40);

        for name in &["C", "DC", "CMP"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(
                scheme.structure_color(&r),
                expected,
                "Cytosine variant {name}"
            );
        }
    }

    #[test]
    fn structure_color_inosine_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(150, 100, 180);

        for name in &["I", "DI"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(
                scheme.structure_color(&r),
                expected,
                "Inosine variant {name}"
            );
        }
    }

    #[test]
    fn is_nucleotide_inosine() {
        assert!(is_nucleotide("I"), "I should be recognized as nucleotide");
        assert!(is_nucleotide("DI"), "DI should be recognized as nucleotide");
    }

    // ---- structure_color: protein residues still get secondary-structure colors ----

    #[test]
    fn structure_color_protein_helix() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let r = make_residue("ALA", SecondaryStructure::Helix);
        assert_eq!(scheme.structure_color(&r), Color::Rgb(255, 0, 128));
    }

    #[test]
    fn structure_color_protein_sheet() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let r = make_residue("GLY", SecondaryStructure::Sheet);
        assert_eq!(scheme.structure_color(&r), Color::Rgb(255, 200, 0));
    }

    #[test]
    fn structure_color_protein_turn() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let r = make_residue("CYS", SecondaryStructure::Turn);
        assert_eq!(scheme.structure_color(&r), Color::Rgb(96, 128, 255));
    }

    #[test]
    fn structure_color_protein_coil() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let r = make_residue("LEU", SecondaryStructure::Coil);
        assert_eq!(scheme.structure_color(&r), Color::Rgb(0, 204, 0));
    }

    #[test]
    fn plddt_color_uses_standard_af_bands() {
        let scheme = ColorScheme::new(ColorSchemeType::Plddt, 100);

        let mk = |score: f64| Residue {
            name: "ALA".to_string(),
            seq_num: 1,
            atoms: vec![Atom {
                name: "CA".to_string(),
                element: "C".to_string(),
                x: 0.0,
                y: 0.0,
                z: 0.0,
                b_factor: score,
                is_backbone: true,
            }],
            secondary_structure: SecondaryStructure::Coil,
        };

        assert_eq!(scheme.plddt_color(&mk(95.0)), Color::Rgb(0, 83, 214));
        assert_eq!(scheme.plddt_color(&mk(80.0)), Color::Rgb(101, 203, 243));
        assert_eq!(scheme.plddt_color(&mk(60.0)), Color::Rgb(255, 219, 19));
        assert_eq!(scheme.plddt_color(&mk(40.0)), Color::Rgb(255, 125, 69));
    }
}

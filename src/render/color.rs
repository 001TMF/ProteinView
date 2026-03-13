use std::collections::HashSet;

use ratatui::style::Color;
use crate::model::interface::InterfaceAnalysis;
use crate::model::protein::{Atom, Chain, Ligand, LigandType, Residue, SecondaryStructure};

/// Available color schemes
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ColorSchemeType {
    Structure,
    Chain,
    Element,
    BFactor,
    Rainbow,
    Interface,
    Plddt,
}

impl ColorSchemeType {
    pub fn next(&self, has_plddt: bool) -> Self {
        match self {
            Self::Structure => Self::Element,
            Self::Element => Self::Chain,
            Self::Chain => Self::BFactor,
            Self::BFactor => Self::Rainbow,
            Self::Rainbow => {
                if has_plddt {
                    Self::Plddt
                } else {
                    Self::Structure
                }
            }
            Self::Plddt => Self::Structure,
            // Interface is toggled separately via 'f' key
            Self::Interface => Self::Structure,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Self::Structure => "Structure",
            Self::Chain => "Chain",
            Self::Element => "Element",
            Self::BFactor => "B-Factor",
            Self::Rainbow => "Rainbow",
            Self::Interface => "Interface",
            Self::Plddt => "pLDDT",
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
        let focus_chain_id = protein.chains.get(focus_chain)
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
            ColorSchemeType::Chain => self.chain_color(chain),
            ColorSchemeType::Element => Color::Rgb(144, 144, 144),
            ColorSchemeType::BFactor => self.bfactor_color(residue),
            ColorSchemeType::Rainbow => self.rainbow_color(residue),
            ColorSchemeType::Interface => self.interface_color(residue, chain),
            ColorSchemeType::Plddt => self.plddt_residue_color(residue),
        }
    }

    /// Get color for an individual atom, respecting the current scheme.
    pub fn atom_color(&self, atom: &Atom, residue: &Residue, chain: &Chain) -> Color {
        match self.scheme_type {
            ColorSchemeType::Element => Self::element_color(atom),
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
        let is_contact = self.interface_residues_by_id
            .contains(&(chain.id.clone(), residue.seq_num));
        let is_focus = chain.id == self.focus_chain_id;

        match (is_focus, is_contact) {
            (true, true) => Color::Rgb(0, 255, 100),    // Bright green — antibody interface
            (true, false) => Color::Rgb(40, 100, 60),   // Dim green — antibody non-interface
            (false, true) => Color::Rgb(255, 165, 0),   // Bright orange — antigen interface
            (false, false) => Color::Rgb(100, 80, 60),  // Dim brown — antigen non-interface
        }
    }

    /// CPK-style element coloring
    pub fn element_color(atom: &Atom) -> Color {
        match atom.element.trim() {
            "C" => Color::Rgb(144, 144, 144),
            "N" => Color::Rgb(48, 80, 248),
            "O" => Color::Rgb(255, 13, 13),
            "S" => Color::Rgb(255, 255, 48),
            "H" => Color::Rgb(255, 255, 255),
            "P" => Color::Rgb(255, 128, 0),
            "FE" | "Fe" => Color::Rgb(224, 102, 51),
            "MG" | "Mg" => Color::Rgb(0, 180, 0),        // Magnesium — green
            "ZN" | "Zn" => Color::Rgb(125, 128, 176),    // Zinc — blue-gray
            "CA" | "Ca" => Color::Rgb(61, 255, 0),       // Calcium — green
            "MN" | "Mn" => Color::Rgb(156, 122, 199),    // Manganese — purple
            "CO" | "Co" => Color::Rgb(240, 144, 160),    // Cobalt — pink
            "CU" | "Cu" => Color::Rgb(200, 128, 51),     // Copper — brown-orange
            "NI" | "Ni" => Color::Rgb(80, 208, 80),      // Nickel — green
            "CL" | "Cl" => Color::Rgb(31, 240, 31),      // Chlorine — green
            "BR" | "Br" => Color::Rgb(166, 41, 41),      // Bromine — dark red
            _ => Color::Rgb(200, 200, 200),
        }
    }

    /// Get base color for a ligand based on current scheme.
    pub fn ligand_color(&self, ligand: &Ligand) -> Color {
        match self.scheme_type {
            ColorSchemeType::Structure => match ligand.ligand_type {
                LigandType::Ligand => Color::Rgb(255, 0, 255),   // magenta for ligands
                LigandType::Ion => Color::Rgb(0, 255, 255),      // cyan for ions
            },
            ColorSchemeType::Element => Color::Rgb(144, 144, 144), // overridden per-atom
            ColorSchemeType::BFactor => {
                let avg_b = if ligand.atoms.is_empty() {
                    0.0
                } else {
                    ligand.atoms.iter().map(|a| a.b_factor).sum::<f64>() / ligand.atoms.len() as f64
                };
                let t = ((avg_b - 5.0) / 75.0).clamp(0.0, 1.0);
                let r = (t * 255.0) as u8;
                let b = ((1.0 - t) * 255.0) as u8;
                Color::Rgb(r, 0, b)
            },
            ColorSchemeType::Chain => {
                // Match parent chain's color using chain_id
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
                let idx = ligand.chain_id.bytes().next().unwrap_or(b'A') as usize % chain_colors.len();
                chain_colors[idx]
            },
            ColorSchemeType::Rainbow => Color::Rgb(255, 0, 255),
            ColorSchemeType::Interface => Color::Rgb(255, 255, 255), // bright white to stand out
            // pLDDT mode: fall back to Structure-mode colors for ligands
            ColorSchemeType::Plddt => match ligand.ligand_type {
                LigandType::Ligand => Color::Rgb(255, 0, 255),   // magenta for ligands
                LigandType::Ion => Color::Rgb(0, 255, 255),      // cyan for ions
            },
        }
    }

    /// Get color for a specific atom within a ligand.
    pub fn ligand_atom_color(&self, atom: &Atom, ligand: &Ligand) -> Color {
        match self.scheme_type {
            ColorSchemeType::Element => Self::element_color(atom),
            // pLDDT mode: fall back to Element-mode (CPK) colors for ligand atoms
            ColorSchemeType::Plddt => Self::element_color(atom),
            _ => self.ligand_color(ligand),
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
            SecondaryStructure::Turn  => Color::Rgb(96, 128, 255),
            SecondaryStructure::Coil  => Color::Rgb(0, 204, 0),
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

    fn rainbow_color(&self, residue: &Residue) -> Color {
        if self.total_residues == 0 { return Color::White; }
        let t = residue.seq_num as f64 / self.total_residues as f64;
        let hue = (1.0 - t) * 300.0;
        let (r, g, b) = hsv_to_rgb(hue, 1.0, 1.0);
        Color::Rgb(r, g, b)
    }

    /// AlphaFold pLDDT confidence color for a raw B-factor value.
    ///
    /// Uses the standard AlphaFold color bands:
    /// - >= 90: dark blue  (very high confidence)
    /// - >= 70: light blue (confident)
    /// - >= 50: yellow     (low confidence)
    /// - <  50: orange     (very low confidence)
    pub fn plddt_color(b_factor: f64) -> Color {
        if b_factor >= 90.0 {
            Color::Rgb(0, 83, 214)
        } else if b_factor >= 70.0 {
            Color::Rgb(101, 203, 243)
        } else if b_factor >= 50.0 {
            Color::Rgb(255, 219, 19)
        } else {
            Color::Rgb(255, 125, 69)
        }
    }

    /// pLDDT color for a residue, using the average B-factor of its atoms.
    fn plddt_residue_color(&self, residue: &Residue) -> Color {
        let avg_b = if residue.atoms.is_empty() {
            0.0
        } else {
            residue.atoms.iter().map(|a| a.b_factor).sum::<f64>() / residue.atoms.len() as f64
        };
        Self::plddt_color(avg_b)
    }
}

/// Returns a base-type color for nucleotide residues, or `None` for non-nucleotides.
fn nucleotide_base_color(name: &str) -> Option<Color> {
    match name {
        "A" | "DA" | "AMP" => Some(Color::Rgb(220, 60, 60)),   // Adenine — red
        "U" | "UMP"        => Some(Color::Rgb(60, 60, 220)),   // Uracil — blue
        "T" | "DT"         => Some(Color::Rgb(60, 60, 220)),   // Thymine — blue
        "G" | "DG" | "GMP" => Some(Color::Rgb(60, 180, 60)),   // Guanine — green
        "C" | "DC" | "CMP" => Some(Color::Rgb(220, 200, 40)),  // Cytosine — yellow
        "I" | "DI"         => Some(Color::Rgb(150, 100, 180)), // Inosine — purple
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
    (((r + m) * 255.0) as u8, ((g + m) * 255.0) as u8, ((b + m) * 255.0) as u8)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::protein::{is_nucleotide, LigandType, Residue, SecondaryStructure};

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
            assert!(is_nucleotide(name), "{name} should be recognized as nucleotide");
        }
    }

    #[test]
    fn is_nucleotide_dna_bases() {
        for name in &["DA", "DT", "DG", "DC"] {
            assert!(is_nucleotide(name), "{name} should be recognized as nucleotide");
        }
    }

    #[test]
    fn is_nucleotide_modified_forms() {
        for name in &["AMP", "UMP", "GMP", "CMP"] {
            assert!(is_nucleotide(name), "{name} should be recognized as nucleotide");
        }
    }

    #[test]
    fn is_nucleotide_rejects_amino_acids() {
        for name in &["ALA", "GLY", "CYS", "THR", "TRP", "LEU"] {
            assert!(!is_nucleotide(name), "{name} should NOT be recognized as nucleotide");
        }
    }

    // ---- structure_color: nucleotide residues ----

    #[test]
    fn structure_color_adenine_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(220, 60, 60);

        for name in &["A", "DA", "AMP"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(scheme.structure_color(&r), expected, "Adenine variant {name}");
        }
    }

    #[test]
    fn structure_color_uracil_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(60, 60, 220);

        for name in &["U", "UMP"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(scheme.structure_color(&r), expected, "Uracil variant {name}");
        }
    }

    #[test]
    fn structure_color_thymine_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(60, 60, 220);

        for name in &["T", "DT"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(scheme.structure_color(&r), expected, "Thymine variant {name}");
        }
    }

    #[test]
    fn structure_color_guanine_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(60, 180, 60);

        for name in &["G", "DG", "GMP"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(scheme.structure_color(&r), expected, "Guanine variant {name}");
        }
    }

    #[test]
    fn structure_color_cytosine_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(220, 200, 40);

        for name in &["C", "DC", "CMP"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(scheme.structure_color(&r), expected, "Cytosine variant {name}");
        }
    }

    #[test]
    fn structure_color_inosine_variants() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let expected = Color::Rgb(150, 100, 180);

        for name in &["I", "DI"] {
            let r = make_residue(name, SecondaryStructure::Coil);
            assert_eq!(scheme.structure_color(&r), expected, "Inosine variant {name}");
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

    // ---- ligand coloring ----

    #[test]
    fn ligand_color_structure_magenta() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let ligand = crate::model::protein::Ligand {
            name: "HEM".to_string(),
            chain_id: "A".to_string(),
            seq_num: 1,
            atoms: vec![],
            ligand_type: LigandType::Ligand,
        };
        assert_eq!(scheme.ligand_color(&ligand), Color::Rgb(255, 0, 255));
    }

    #[test]
    fn ligand_color_ion_cyan() {
        let scheme = ColorScheme::new(ColorSchemeType::Structure, 100);
        let ligand = crate::model::protein::Ligand {
            name: "ZN".to_string(),
            chain_id: "A".to_string(),
            seq_num: 1,
            atoms: vec![],
            ligand_type: LigandType::Ion,
        };
        assert_eq!(scheme.ligand_color(&ligand), Color::Rgb(0, 255, 255));
    }

    #[test]
    fn ligand_atom_color_element_mode() {
        let scheme = ColorScheme::new(ColorSchemeType::Element, 100);
        let atom = Atom {
            name: "FE".to_string(),
            element: "Fe".to_string(),
            x: 0.0, y: 0.0, z: 0.0,
            b_factor: 0.0,
            is_backbone: false,
            is_hetero: true,
        };
        let ligand = crate::model::protein::Ligand {
            name: "HEM".to_string(),
            chain_id: "A".to_string(),
            seq_num: 1,
            atoms: vec![],
            ligand_type: LigandType::Ligand,
        };
        // In Element mode, should use element-based CPK color for Fe
        assert_eq!(scheme.ligand_atom_color(&atom, &ligand), Color::Rgb(224, 102, 51));
    }

    #[test]
    fn element_color_zinc() {
        let atom = Atom {
            name: "ZN".to_string(),
            element: "Zn".to_string(),
            x: 0.0, y: 0.0, z: 0.0,
            b_factor: 0.0,
            is_backbone: false,
            is_hetero: true,
        };
        assert_eq!(ColorScheme::element_color(&atom), Color::Rgb(125, 128, 176));
    }

    #[test]
    fn element_color_magnesium() {
        let atom = Atom {
            name: "MG".to_string(),
            element: "Mg".to_string(),
            x: 0.0, y: 0.0, z: 0.0,
            b_factor: 0.0,
            is_backbone: false,
            is_hetero: true,
        };
        assert_eq!(ColorScheme::element_color(&atom), Color::Rgb(0, 180, 0));
    }

    // ---- pLDDT coloring ----

    #[test]
    fn test_plddt_color_bands() {
        // Very high confidence (>= 90): dark blue
        assert_eq!(ColorScheme::plddt_color(90.0), Color::Rgb(0, 83, 214));
        assert_eq!(ColorScheme::plddt_color(95.0), Color::Rgb(0, 83, 214));
        assert_eq!(ColorScheme::plddt_color(100.0), Color::Rgb(0, 83, 214));

        // Confident (>= 70, < 90): light blue
        assert_eq!(ColorScheme::plddt_color(70.0), Color::Rgb(101, 203, 243));
        assert_eq!(ColorScheme::plddt_color(80.0), Color::Rgb(101, 203, 243));
        assert_eq!(ColorScheme::plddt_color(89.9), Color::Rgb(101, 203, 243));

        // Low confidence (>= 50, < 70): yellow
        assert_eq!(ColorScheme::plddt_color(50.0), Color::Rgb(255, 219, 19));
        assert_eq!(ColorScheme::plddt_color(60.0), Color::Rgb(255, 219, 19));
        assert_eq!(ColorScheme::plddt_color(69.9), Color::Rgb(255, 219, 19));

        // Very low confidence (< 50): orange
        assert_eq!(ColorScheme::plddt_color(49.9), Color::Rgb(255, 125, 69));
        assert_eq!(ColorScheme::plddt_color(30.0), Color::Rgb(255, 125, 69));
        assert_eq!(ColorScheme::plddt_color(0.0), Color::Rgb(255, 125, 69));
    }

    #[test]
    fn test_plddt_next_cycling() {
        // With pLDDT: Structure -> Element -> Chain -> BFactor -> Rainbow -> Plddt -> Structure
        assert_eq!(ColorSchemeType::Structure.next(true), ColorSchemeType::Element);
        assert_eq!(ColorSchemeType::Element.next(true), ColorSchemeType::Chain);
        assert_eq!(ColorSchemeType::Chain.next(true), ColorSchemeType::BFactor);
        assert_eq!(ColorSchemeType::BFactor.next(true), ColorSchemeType::Rainbow);
        assert_eq!(ColorSchemeType::Rainbow.next(true), ColorSchemeType::Plddt);
        assert_eq!(ColorSchemeType::Plddt.next(true), ColorSchemeType::Structure);

        // Without pLDDT: Structure -> Element -> Chain -> BFactor -> Rainbow -> Structure
        assert_eq!(ColorSchemeType::Structure.next(false), ColorSchemeType::Element);
        assert_eq!(ColorSchemeType::Element.next(false), ColorSchemeType::Chain);
        assert_eq!(ColorSchemeType::Chain.next(false), ColorSchemeType::BFactor);
        assert_eq!(ColorSchemeType::BFactor.next(false), ColorSchemeType::Rainbow);
        assert_eq!(ColorSchemeType::Rainbow.next(false), ColorSchemeType::Structure);

        // Interface always cycles to Structure regardless of pLDDT flag
        assert_eq!(ColorSchemeType::Interface.next(true), ColorSchemeType::Structure);
        assert_eq!(ColorSchemeType::Interface.next(false), ColorSchemeType::Structure);
    }

    #[test]
    fn test_plddt_ligand_fallback() {
        let scheme = ColorScheme::new(ColorSchemeType::Plddt, 100);

        // Ligands should fall back to Structure-mode colors (magenta)
        let ligand = crate::model::protein::Ligand {
            name: "HEM".to_string(),
            chain_id: "A".to_string(),
            seq_num: 1,
            atoms: vec![],
            ligand_type: LigandType::Ligand,
        };
        assert_eq!(scheme.ligand_color(&ligand), Color::Rgb(255, 0, 255));

        // Ions should fall back to Structure-mode colors (cyan)
        let ion = crate::model::protein::Ligand {
            name: "ZN".to_string(),
            chain_id: "A".to_string(),
            seq_num: 1,
            atoms: vec![],
            ligand_type: LigandType::Ion,
        };
        assert_eq!(scheme.ligand_color(&ion), Color::Rgb(0, 255, 255));

        // Ligand atoms should fall back to Element-mode (CPK) colors
        let fe_atom = Atom {
            name: "FE".to_string(),
            element: "Fe".to_string(),
            x: 0.0, y: 0.0, z: 0.0,
            b_factor: 0.0,
            is_backbone: false,
            is_hetero: true,
        };
        assert_eq!(scheme.ligand_atom_color(&fe_atom, &ligand), Color::Rgb(224, 102, 51));
    }
}

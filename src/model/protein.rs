/// Classification of the polymer type for a chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(clippy::upper_case_acronyms)]
pub enum MoleculeType {
    Protein,
    RNA,
    DNA,
    #[allow(dead_code)]
    SmallMolecule,
}

/// Distinguishes multi-atom ligands from single-atom ions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LigandType {
    Ligand,
    Ion,
}

/// Standard RNA residue names.
pub const RNA_RESIDUES: &[&str] = &["A", "U", "G", "C", "I", "AMP", "UMP", "GMP", "CMP"];

/// Standard DNA residue names.
pub const DNA_RESIDUES: &[&str] = &["DA", "DT", "DG", "DC", "DI", "T"];

/// Residue names used for crystallographic water molecules.
pub const WATER_NAMES: &[&str] = &["HOH", "WAT", "DOD", "H2O", "OH2"];

/// Common single-atom ions found as HETATM in PDB files.
pub const COMMON_IONS: &[&str] = &[
    "ZN", "MG", "CA", "FE", "MN", "CO", "CU", "NI", "CD",
    "NA", "K", "CL", "BR", "I", "F",
    "HG", "PT", "AU", "AG", "PB",
];

/// Returns true if the residue name is a nucleotide (RNA or DNA).
#[allow(dead_code)]
pub fn is_nucleotide(name: &str) -> bool {
    RNA_RESIDUES.contains(&name) || DNA_RESIDUES.contains(&name)
}

/// Returns true if the residue name is a purine base (A, G, I and their variants).
pub fn is_purine(name: &str) -> bool {
    matches!(name, "A" | "DA" | "AMP" | "G" | "DG" | "GMP" | "I" | "DI")
}

/// A complete protein structure
#[derive(Debug, Clone)]
pub struct Protein {
    pub name: String,
    pub chains: Vec<Chain>,
    pub ligands: Vec<Ligand>,
}

/// A polypeptide chain
#[derive(Debug, Clone)]
pub struct Chain {
    pub id: String,
    pub residues: Vec<Residue>,
    pub molecule_type: MoleculeType,
}

/// A small molecule (ligand, cofactor, or ion) from HETATM records.
#[derive(Debug, Clone)]
pub struct Ligand {
    pub name: String,
    pub chain_id: String,
    pub seq_num: i32,
    pub atoms: Vec<Atom>,
    pub ligand_type: LigandType,
}

/// An amino acid residue
#[derive(Debug, Clone)]
pub struct Residue {
    pub name: String,
    pub seq_num: i32,
    pub atoms: Vec<Atom>,
    pub secondary_structure: SecondaryStructure,
}

/// An individual atom
#[derive(Debug, Clone)]
pub struct Atom {
    pub name: String,
    pub element: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub b_factor: f64,
    pub is_backbone: bool,
    #[allow(dead_code)]
    pub is_hetero: bool,
}

/// Secondary structure classification
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SecondaryStructure {
    Helix,
    Sheet,
    #[allow(dead_code)]
    Turn,
    Coil,
}

impl Protein {
    /// Get all backbone atoms (C-alpha for proteins, C4' for nucleic acids)
    /// for backbone trace rendering.
    pub fn backbone_atoms(&self) -> Vec<(&Atom, &Residue, &Chain)> {
        let mut cas = Vec::new();
        for chain in &self.chains {
            for residue in &chain.residues {
                for atom in &residue.atoms {
                    if atom.is_backbone {
                        cas.push((atom, residue, chain));
                    }
                }
            }
        }
        cas
    }

    /// Get total atom count
    pub fn atom_count(&self) -> usize {
        self.chains.iter()
            .flat_map(|c| &c.residues)
            .flat_map(|r| &r.atoms)
            .count()
            + self.ligands.iter().flat_map(|l| &l.atoms).count()
    }

    /// Get total residue count
    pub fn residue_count(&self) -> usize {
        self.chains.iter()
            .flat_map(|c| &c.residues)
            .count()
    }

    /// Get the bounding radius from origin (call after centering)
    pub fn bounding_radius(&self) -> f64 {
        let chain_atoms = self.chains.iter()
            .flat_map(|c| &c.residues)
            .flat_map(|r| &r.atoms)
            .filter(|a| a.is_backbone);
        let ligand_atoms = self.ligands.iter()
            .flat_map(|l| &l.atoms);
        chain_atoms.map(|a| (a.x * a.x + a.y * a.y + a.z * a.z).sqrt())
            .chain(ligand_atoms.map(|a| (a.x * a.x + a.y * a.y + a.z * a.z).sqrt()))
            .fold(0.0f64, f64::max)
    }

    /// Center the protein at the origin
    pub fn center(&mut self) {
        let chain_atoms: Vec<&Atom> = self.chains.iter()
            .flat_map(|c| &c.residues)
            .flat_map(|r| &r.atoms)
            .collect();
        let ligand_atoms: Vec<&Atom> = self.ligands.iter()
            .flat_map(|l| &l.atoms)
            .collect();

        let total = chain_atoms.len() + ligand_atoms.len();
        if total == 0 { return; }

        let n = total as f64;
        let cx: f64 = (chain_atoms.iter().map(|a| a.x).sum::<f64>()
            + ligand_atoms.iter().map(|a| a.x).sum::<f64>()) / n;
        let cy: f64 = (chain_atoms.iter().map(|a| a.y).sum::<f64>()
            + ligand_atoms.iter().map(|a| a.y).sum::<f64>()) / n;
        let cz: f64 = (chain_atoms.iter().map(|a| a.z).sum::<f64>()
            + ligand_atoms.iter().map(|a| a.z).sum::<f64>()) / n;

        for chain in &mut self.chains {
            for residue in &mut chain.residues {
                for atom in &mut residue.atoms {
                    atom.x -= cx;
                    atom.y -= cy;
                    atom.z -= cz;
                }
            }
        }
        for ligand in &mut self.ligands {
            for atom in &mut ligand.atoms {
                atom.x -= cx;
                atom.y -= cy;
                atom.z -= cz;
            }
        }
    }

    /// Get total number of ligands (including ions)
    pub fn ligand_count(&self) -> usize {
        self.ligands.len()
    }

    /// Get total number of atoms across all ligands
    #[allow(dead_code)]
    pub fn ligand_atom_count(&self) -> usize {
        self.ligands.iter().flat_map(|l| &l.atoms).count()
    }

    /// Heuristically detect whether the B-factor column stores pLDDT scores.
    ///
    /// AlphaFold/ModelCIF outputs store confidence values in [0, 100], with
    /// most atoms above 50 and many above 70.  Classic experimental B-factors
    /// are usually much lower on average even when they overlap numerically.
    ///
    /// Only polymer chain atoms are considered (ligands are excluded).
    pub fn has_plddt(&self) -> bool {
        let mut total = 0usize;
        let mut in_range = 0usize;
        let mut high_conf = 0usize;
        let mut sum = 0.0f64;

        for chain in &self.chains {
            if chain.molecule_type == MoleculeType::SmallMolecule {
                continue;
            }
            for residue in &chain.residues {
                for atom in &residue.atoms {
                    total += 1;
                    let value = atom.b_factor;
                    sum += value;
                    if (0.0..=100.0).contains(&value) {
                        in_range += 1;
                    }
                    if value >= 70.0 {
                        high_conf += 1;
                    }
                }
            }
        }

        if total == 0 {
            return false;
        }

        let mean = sum / total as f64;
        let in_range_fraction = in_range as f64 / total as f64;
        let high_conf_fraction = high_conf as f64 / total as f64;

        // All three conditions must be true:
        // 1) >= 95% of B-factors in [0, 100]
        // 2) Mean B-factor >= 50
        // 3) >= 25% of atoms have B-factor >= 70
        in_range_fraction >= 0.95 && mean >= 50.0 && high_conf_fraction >= 0.25
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom_with_bfactor(b: f64) -> Atom {
        Atom {
            name: "CA".to_string(),
            element: "C".to_string(),
            x: 0.0, y: 0.0, z: 0.0,
            b_factor: b,
            is_backbone: true,
            is_hetero: false,
        }
    }

    fn protein_from_bfactors(values: &[f64]) -> Protein {
        Protein {
            name: "test".to_string(),
            chains: vec![Chain {
                id: "A".to_string(),
                molecule_type: MoleculeType::Protein,
                residues: values.iter().enumerate().map(|(i, &b)| Residue {
                    name: "ALA".to_string(),
                    seq_num: i as i32 + 1,
                    atoms: vec![atom_with_bfactor(b)],
                    secondary_structure: SecondaryStructure::Coil,
                }).collect(),
            }],
            ligands: vec![],
        }
    }

    #[test]
    fn test_has_plddt_alphafold_like() {
        // Typical AlphaFold pLDDT scores: mostly high, all in [0,100]
        let protein = protein_from_bfactors(&[
            95.0, 92.0, 88.0, 76.0, 67.0, 54.0, 91.0, 85.0, 73.0, 80.0,
        ]);
        assert!(protein.has_plddt());
    }

    #[test]
    fn test_has_plddt_crystallographic() {
        // Typical crystallographic B-factors: low values, wide range
        let protein = protein_from_bfactors(&[
            12.0, 18.0, 22.0, 30.0, 16.0, 25.0, 8.0, 14.0, 20.0, 35.0,
        ]);
        assert!(!protein.has_plddt());
    }

    #[test]
    fn test_has_plddt_empty_protein() {
        let protein = Protein {
            name: "empty".to_string(),
            chains: vec![],
            ligands: vec![],
        };
        assert!(!protein.has_plddt());
    }

    #[test]
    fn test_has_plddt_borderline_mean_below_threshold() {
        // Mean = 49.0, just below 50.0 threshold — should reject
        let protein = protein_from_bfactors(&[
            45.0, 42.0, 65.0, 48.0, 40.0, 55.0, 50.0, 38.0, 60.0, 47.0,
        ]);
        assert!(!protein.has_plddt());
    }

    #[test]
    fn test_has_plddt_negative_bfactors() {
        // NMR structures can have negative B-factors — outside [0,100]
        let protein = protein_from_bfactors(&[
            -5.0, 80.0, 75.0, 90.0, 85.0, -2.0, 78.0, 92.0, 70.0, 88.0,
        ]);
        assert!(!protein.has_plddt());
    }
}

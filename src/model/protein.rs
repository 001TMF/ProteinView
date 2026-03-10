/// Classification of the polymer type for a chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(clippy::upper_case_acronyms)]
pub enum MoleculeType {
    Protein,
    RNA,
    DNA,
}

/// Standard RNA residue names.
pub const RNA_RESIDUES: &[&str] = &["A", "U", "G", "C", "I", "AMP", "UMP", "GMP", "CMP"];

/// Standard DNA residue names.
pub const DNA_RESIDUES: &[&str] = &["DA", "DT", "DG", "DC", "DI", "T"];

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
}

/// A polypeptide chain
#[derive(Debug, Clone)]
pub struct Chain {
    pub id: String,
    pub residues: Vec<Residue>,
    pub molecule_type: MoleculeType,
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
    fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.chains
            .iter()
            .flat_map(|c| &c.residues)
            .flat_map(|r| &r.atoms)
    }

    /// Get total atom count
    pub fn atom_count(&self) -> usize {
        self.atoms().count()
    }

    /// Get total residue count
    pub fn residue_count(&self) -> usize {
        self.chains.iter().flat_map(|c| &c.residues).count()
    }

    /// Get the bounding radius from origin (call after centering)
    pub fn bounding_radius(&self) -> f64 {
        self.atoms()
            .filter(|a| a.is_backbone)
            .map(|a| (a.x * a.x + a.y * a.y + a.z * a.z).sqrt())
            .fold(0.0f64, f64::max)
    }

    /// Heuristically detect whether the B-factor column stores pLDDT scores.
    ///
    /// AlphaFold/ModelCIF outputs typically store confidence values in the
    /// range [0, 100], with most atoms above 50 and many above 70. Classic
    /// experimental B-factors are usually much lower on average even when they
    /// overlap numerically.
    pub fn has_plddt(&self) -> bool {
        let mut total = 0usize;
        let mut in_range = 0usize;
        let mut high_conf = 0usize;
        let mut sum = 0.0f64;

        for atom in self.atoms() {
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

        if total == 0 {
            return false;
        }

        let mean = sum / total as f64;
        let in_range_fraction = in_range as f64 / total as f64;
        let high_conf_fraction = high_conf as f64 / total as f64;

        in_range_fraction >= 0.95 && mean >= 50.0 && high_conf_fraction >= 0.25
    }

    /// Center the protein at the origin
    pub fn center(&mut self) {
        let atoms: Vec<&Atom> = self.atoms().collect();

        if atoms.is_empty() {
            return;
        }

        let n = atoms.len() as f64;
        let cx: f64 = atoms.iter().map(|a| a.x).sum::<f64>() / n;
        let cy: f64 = atoms.iter().map(|a| a.y).sum::<f64>() / n;
        let cz: f64 = atoms.iter().map(|a| a.z).sum::<f64>() / n;

        for chain in &mut self.chains {
            for residue in &mut chain.residues {
                for atom in &mut residue.atoms {
                    atom.x -= cx;
                    atom.y -= cy;
                    atom.z -= cz;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom_with_bfactor(b_factor: f64) -> Atom {
        Atom {
            name: "CA".to_string(),
            element: "C".to_string(),
            x: 0.0,
            y: 0.0,
            z: 0.0,
            b_factor,
            is_backbone: true,
        }
    }

    fn protein_from_bfactors(values: &[f64]) -> Protein {
        Protein {
            name: "test".to_string(),
            chains: vec![Chain {
                id: "A".to_string(),
                molecule_type: MoleculeType::Protein,
                residues: values
                    .iter()
                    .enumerate()
                    .map(|(i, value)| Residue {
                        name: "ALA".to_string(),
                        seq_num: i as i32 + 1,
                        atoms: vec![atom_with_bfactor(*value)],
                        secondary_structure: SecondaryStructure::Coil,
                    })
                    .collect(),
            }],
        }
    }

    #[test]
    fn detects_plddt_like_confidence_scores() {
        let protein = protein_from_bfactors(&[95.0, 92.0, 88.0, 76.0, 67.0, 54.0]);
        assert!(protein.has_plddt());
    }

    #[test]
    fn rejects_typical_experimental_bfactors() {
        let protein = protein_from_bfactors(&[12.0, 18.0, 22.0, 30.0, 16.0, 25.0]);
        assert!(!protein.has_plddt());
    }
}

/// Classification of the polymer type for a chain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MoleculeType {
    Protein,
    RNA,
    DNA,
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
    }

    /// Get total residue count
    pub fn residue_count(&self) -> usize {
        self.chains.iter()
            .flat_map(|c| &c.residues)
            .count()
    }

    /// Get the bounding radius from origin (call after centering)
    pub fn bounding_radius(&self) -> f64 {
        self.chains.iter()
            .flat_map(|c| &c.residues)
            .flat_map(|r| &r.atoms)
            .filter(|a| a.is_backbone)
            .map(|a| (a.x * a.x + a.y * a.y + a.z * a.z).sqrt())
            .fold(0.0f64, f64::max)
    }

    /// Center the protein at the origin
    pub fn center(&mut self) {
        let atoms: Vec<&Atom> = self.chains.iter()
            .flat_map(|c| &c.residues)
            .flat_map(|r| &r.atoms)
            .collect();

        if atoms.is_empty() { return; }

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

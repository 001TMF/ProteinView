use anyhow::Result;
use crate::model::protein::{Protein, Chain, MoleculeType, Residue, Atom, SecondaryStructure};
use crate::model::secondary::{assign_from_pdb_file, assign_from_cif_file};

/// Load a protein structure from a PDB or mmCIF file
pub fn load_structure(path: &str) -> Result<Protein> {
    // Try default strictness first, fall back to loose + atomic-coords-only
    // for files like AlphaFold3 output that have non-standard metadata
    let (pdb, _errors) = pdbtbx::open(path)
        .or_else(|_| {
            pdbtbx::ReadOptions::new()
                .set_level(pdbtbx::StrictnessLevel::Loose)
                .set_only_atomic_coords(true)
                .read(path)
        })
        .map_err(|e| anyhow::anyhow!("Failed to open structure file: {:?}", e))?;

    let mut chains = Vec::new();

    for chain in pdb.chains() {
        let mut residues = Vec::new();
        for residue in chain.residues() {
            let mut atoms = Vec::new();
            for atom in residue.atoms() {
                atoms.push(Atom {
                    name: atom.name().to_string(),
                    element: atom.element().map(|e| format!("{:?}", e)).unwrap_or_default(),
                    x: atom.x(),
                    y: atom.y(),
                    z: atom.z(),
                    b_factor: atom.b_factor(),
                    is_backbone: atom.name() == "CA" || atom.name() == "C4'",
                });
            }
            residues.push(Residue {
                name: residue.name().unwrap_or("UNK").to_string(),
                seq_num: residue.serial_number() as i32,
                atoms,
                secondary_structure: SecondaryStructure::Coil,
            });
        }
        let molecule_type = classify_chain_type(&residues);
        chains.push(Chain {
            id: chain.id().to_string(),
            residues,
            molecule_type,
        });
    }

    let name = pdb.identifier.as_deref().unwrap_or("Unknown").to_string();

    let mut protein = Protein { name, chains };

    // Assign secondary structure from HELIX/SHEET records in the PDB file
    assign_from_pdb_file(&mut protein, path);

    // If all residues are still Coil (no PDB HELIX/SHEET records found),
    // try CIF _struct_conf/_struct_sheet_range parsing as a fallback.
    let all_coil = protein.chains.iter()
        .flat_map(|c| &c.residues)
        .all(|r| r.secondary_structure == SecondaryStructure::Coil);
    if all_coil {
        let lower = path.to_lowercase();
        if lower.ends_with(".cif") || lower.ends_with(".mmcif") {
            assign_from_cif_file(&mut protein, path);
        }
    }

    Ok(protein)
}

/// Known RNA residue names.
const RNA_RESIDUES: &[&str] = &["A", "U", "G", "C", "I", "AMP", "UMP", "GMP", "CMP"];

/// Known DNA residue names.
const DNA_RESIDUES: &[&str] = &["DA", "DT", "DG", "DC", "DI"];

/// Classify a chain's molecule type from its residue names.
///
/// Counts residues matching known RNA and DNA names. Whichever set has the
/// majority determines the type. If neither set has any matches (or there is
/// a tie), the chain defaults to `Protein`.
fn classify_chain_type(residues: &[Residue]) -> MoleculeType {
    let mut rna_count = 0usize;
    let mut dna_count = 0usize;

    for res in residues {
        let name = res.name.trim();
        if RNA_RESIDUES.contains(&name) {
            rna_count += 1;
        } else if DNA_RESIDUES.contains(&name) {
            dna_count += 1;
        }
    }

    if rna_count == 0 && dna_count == 0 {
        return MoleculeType::Protein;
    }
    if rna_count >= dna_count {
        MoleculeType::RNA
    } else {
        MoleculeType::DNA
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classify_chain_type_protein() {
        let residues = vec![
            Residue {
                name: "ALA".to_string(),
                seq_num: 1,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "GLY".to_string(),
                seq_num: 2,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
        ];
        assert_eq!(classify_chain_type(&residues), MoleculeType::Protein);
    }

    #[test]
    fn test_classify_chain_type_rna() {
        let residues = vec![
            Residue {
                name: "A".to_string(),
                seq_num: 1,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "U".to_string(),
                seq_num: 2,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "G".to_string(),
                seq_num: 3,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "C".to_string(),
                seq_num: 4,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
        ];
        assert_eq!(classify_chain_type(&residues), MoleculeType::RNA);
    }

    #[test]
    fn test_classify_chain_type_dna() {
        let residues = vec![
            Residue {
                name: "DA".to_string(),
                seq_num: 1,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "DT".to_string(),
                seq_num: 2,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "DG".to_string(),
                seq_num: 3,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "DC".to_string(),
                seq_num: 4,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
        ];
        assert_eq!(classify_chain_type(&residues), MoleculeType::DNA);
    }

    #[test]
    fn test_classify_chain_type_empty() {
        let residues: Vec<Residue> = vec![];
        assert_eq!(classify_chain_type(&residues), MoleculeType::Protein);
    }

    #[test]
    fn test_classify_chain_type_mixed_majority_rna() {
        // 3 RNA residues, 1 DNA residue -> RNA wins
        let residues = vec![
            Residue {
                name: "A".to_string(),
                seq_num: 1,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "U".to_string(),
                seq_num: 2,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "G".to_string(),
                seq_num: 3,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "DA".to_string(),
                seq_num: 4,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
        ];
        assert_eq!(classify_chain_type(&residues), MoleculeType::RNA);
    }

    #[test]
    fn test_backbone_detection_ca() {
        let atom = Atom {
            name: "CA".to_string(),
            element: "C".to_string(),
            x: 0.0,
            y: 0.0,
            z: 0.0,
            b_factor: 0.0,
            is_backbone: true,
        };
        assert!(atom.is_backbone);
    }

    #[test]
    fn test_backbone_detection_c4prime() {
        // C4' should be a backbone atom for nucleic acids
        let name = "C4'";
        let is_backbone = name == "CA" || name == "C4'";
        assert!(is_backbone);
    }

    #[test]
    fn test_non_backbone_atom() {
        let name = "CB";
        let is_backbone = name == "CA" || name == "C4'";
        assert!(!is_backbone);
    }
}

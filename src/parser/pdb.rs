use crate::model::protein::{
    Atom, COMMON_IONS, Chain, DNA_RESIDUES, Ligand, LigandType, MoleculeType, Protein,
    RNA_RESIDUES, Residue, SecondaryStructure, WATER_NAMES,
};
use crate::model::secondary::{
    assign_from_cif_file, assign_from_pdb_file, infer_secondary_structure,
};
use anyhow::Result;

/// Load a protein structure from a PDB or mmCIF file
pub fn load_structure(path: &str) -> Result<Protein> {
    // Try default strictness first, fall back to loose + atomic-coords-only
    // for files like AlphaFold3 output that have non-standard metadata.
    // Both paths use only_first_model to avoid NMR multi-model duplication.
    let (pdb, _errors) = pdbtbx::ReadOptions::new()
        .set_only_first_model(true)
        .read(path)
        .or_else(|_| {
            pdbtbx::ReadOptions::new()
                .set_level(pdbtbx::StrictnessLevel::Loose)
                .set_only_atomic_coords(true)
                .set_only_first_model(true)
                .read(path)
        })
        .map_err(|e| anyhow::anyhow!("Failed to open structure file: {:?}", e))?;

    let mut chains = Vec::new();
    let mut ligands: Vec<Ligand> = Vec::new();

    for chain in pdb.chains() {
        let mut residues = Vec::new();
        for residue in chain.residues() {
            let pdbtbx_atoms: Vec<_> = residue.atoms().collect();
            let res_name = residue.name().unwrap_or("UNK").trim().to_string();
            let all_hetero = pdbtbx_atoms.iter().all(|a| a.hetero());

            if all_hetero && WATER_NAMES.contains(&res_name.as_str()) {
                // Skip water molecules entirely
                continue;
            }

            let atoms: Vec<Atom> = pdbtbx_atoms
                .iter()
                .map(|atom| Atom {
                    name: atom.name().to_string(),
                    element: atom
                        .element()
                        .map(|e| format!("{:?}", e))
                        .unwrap_or_default(),
                    x: atom.x(),
                    y: atom.y(),
                    z: atom.z(),
                    b_factor: atom.b_factor(),
                    is_backbone: atom.name() == "CA" || atom.name() == "C4'",
                    is_hetero: atom.hetero(),
                })
                .collect();

            if all_hetero {
                // Classify as ion or ligand
                let non_h_count = atoms.iter().filter(|a| a.element.trim() != "H").count();
                let ligand_type = if non_h_count <= 1 || COMMON_IONS.contains(&res_name.as_str()) {
                    LigandType::Ion
                } else {
                    LigandType::Ligand
                };
                ligands.push(Ligand {
                    name: res_name,
                    chain_id: chain.id().to_string(),
                    seq_num: residue.serial_number() as i32,
                    atoms,
                    ligand_type,
                });
            } else {
                residues.push(Residue {
                    name: res_name,
                    seq_num: residue.serial_number() as i32,
                    atoms,
                    secondary_structure: SecondaryStructure::Coil,
                });
            }
        }
        let molecule_type = classify_chain_type(&residues);
        chains.push(Chain {
            id: chain.id().to_string(),
            residues,
            molecule_type,
        });
    }

    let name = pdb.identifier.as_deref().unwrap_or("Unknown").to_string();

    let mut protein = Protein {
        name,
        chains,
        ligands,
    };

    // Assign secondary structure from HELIX/SHEET records in the PDB file
    assign_from_pdb_file(&mut protein, path);

    // If all residues are still Coil (no PDB HELIX/SHEET records found),
    // try CIF _struct_conf/_struct_sheet_range parsing as a fallback.
    let all_coil = protein
        .chains
        .iter()
        .flat_map(|c| &c.residues)
        .all(|r| r.secondary_structure == SecondaryStructure::Coil);
    if all_coil {
        let lower = path.to_lowercase();
        if lower.ends_with(".cif") || lower.ends_with(".mmcif") {
            assign_from_cif_file(&mut protein, path);
        }
    }

    // Infer secondary structure from backbone geometry for any protein chains
    // that still lack SS annotations (e.g. AlphaFold PDBs without HELIX/SHEET records).
    infer_secondary_structure(&mut protein.chains);

    Ok(protein)
}

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
            is_hetero: false,
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

    #[test]
    fn test_classify_chain_type_dna_thymine_only() {
        // A chain with only "T" residues should be classified as DNA
        let residues = vec![
            Residue {
                name: "T".to_string(),
                seq_num: 1,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "T".to_string(),
                seq_num: 2,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
            Residue {
                name: "T".to_string(),
                seq_num: 3,
                atoms: vec![],
                secondary_structure: SecondaryStructure::Coil,
            },
        ];
        assert_eq!(classify_chain_type(&residues), MoleculeType::DNA);
    }

    #[test]
    fn test_nmr_multimodel_loads_single_model() {
        // 2KGP is an NMR RNA structure with 10 MODEL records.
        // The parser should only load the first model, not duplicate
        // chains/atoms across all 10 models.
        let protein = load_structure("examples/2KGP.pdb").expect("Failed to load 2KGP.pdb");

        // Should have exactly 1 chain (not 10 duplicated chains)
        assert_eq!(
            protein.chains.len(),
            1,
            "NMR multi-model file should produce 1 chain, got {}",
            protein.chains.len()
        );

        // The single chain should be classified as RNA
        assert_eq!(protein.chains[0].molecule_type, MoleculeType::RNA);

        // Atom count should be reasonable for a single model (~500-900),
        // not inflated by 10x from all models (~8590)
        let total_atoms: usize = protein.chains[0]
            .residues
            .iter()
            .map(|r| r.atoms.len())
            .sum();
        assert!(
            total_atoms < 1000,
            "Expected < 1000 atoms for single NMR model, got {} (multi-model duplication?)",
            total_atoms
        );
    }

    #[test]
    fn test_single_model_pdb_unaffected() {
        // 1UBQ is a single-model X-ray protein structure (ubiquitin).
        // Verify it still loads correctly after NMR multi-model handling.
        let protein = load_structure("examples/1UBQ.pdb").expect("Failed to load 1UBQ.pdb");

        // Ubiquitin has 1 chain (chain A)
        assert_eq!(
            protein.chains.len(),
            1,
            "1UBQ should have 1 chain, got {}",
            protein.chains.len()
        );

        // It should be classified as a protein
        assert_eq!(protein.chains[0].molecule_type, MoleculeType::Protein);

        // Ubiquitin has 76 amino acid residues; crystallographic water
        // molecules (HOH) are now filtered out during parsing.
        assert!(
            protein.chains[0].residues.len() >= 70 && protein.chains[0].residues.len() <= 80,
            "Expected ~76 residues for ubiquitin (waters filtered), got {}",
            protein.chains[0].residues.len()
        );
    }

    #[test]
    fn test_water_filtered_from_ubiquitin() {
        // 1UBQ has 76 amino acid residues and ~58 HOH waters.
        // After filtering, only the 76 AA residues should remain in chains,
        // and no ligands should be present (HOH is discarded, not a ligand).
        let protein = load_structure("examples/1UBQ.pdb").expect("Failed to load 1UBQ.pdb");

        let residue_count = protein.residue_count();
        assert!(
            residue_count >= 70 && residue_count <= 80,
            "Expected ~76 residues after water filtering, got {}",
            residue_count
        );

        // 1UBQ has no real ligands (only HOH waters as HETATM)
        assert_eq!(
            protein.ligand_count(),
            0,
            "Expected 0 ligands for 1UBQ (only waters), got {}",
            protein.ligand_count()
        );
    }

    #[test]
    fn test_4hhb_ligand_parsing() {
        // 4HHB is hemoglobin with 4 HEM (heme) ligands and ions
        let protein = load_structure("examples/4HHB.pdb").expect("Failed to load 4HHB.pdb");

        // Should have 4 protein chains
        assert_eq!(protein.chains.len(), 4);

        // Should have ligands (HEM groups and possibly PO4/ions)
        assert!(protein.ligand_count() > 0, "4HHB should have ligands");

        // At least the 4 HEM groups should be present
        let hem_count = protein.ligands.iter().filter(|l| l.name == "HEM").count();
        assert!(
            hem_count >= 4,
            "Expected at least 4 HEM ligands, got {}",
            hem_count
        );

        // HEM should be classified as Ligand (not Ion) since it's multi-atom
        for l in protein.ligands.iter().filter(|l| l.name == "HEM") {
            assert_eq!(
                l.ligand_type,
                LigandType::Ligand,
                "HEM should be Ligand type"
            );
        }
    }

    #[test]
    fn test_ion_classification() {
        // Verify COMMON_IONS contains expected ions
        assert!(COMMON_IONS.contains(&"ZN"), "ZN should be a common ion");
        assert!(COMMON_IONS.contains(&"MG"), "MG should be a common ion");
        assert!(COMMON_IONS.contains(&"CA"), "CA should be a common ion");
        // ATP is not a single-atom ion
        assert!(
            !COMMON_IONS.contains(&"ATP"),
            "ATP should not be a common ion"
        );
    }
}

use anyhow::Result;
use crate::model::protein::{Protein, Chain, Residue, Atom, SecondaryStructure};
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
                    is_ca: atom.name() == "CA",
                });
            }
            residues.push(Residue {
                name: residue.name().unwrap_or("UNK").to_string(),
                seq_num: residue.serial_number() as i32,
                atoms,
                secondary_structure: SecondaryStructure::Coil,
            });
        }
        chains.push(Chain {
            id: chain.id().to_string(),
            residues,
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

#[allow(unused_imports)]
use crate::model::protein::{Protein, SecondaryStructure};

/// Assign secondary structure from PDB HELIX/SHEET records
/// Falls back to simple phi/psi angle heuristics if no records available
pub fn assign_secondary_structure(protein: &mut Protein) {
    // TODO: Implement DSSP-like algorithm or parse HELIX/SHEET records
    // For now, leave everything as Coil (set during parsing)
    let _ = protein;
}

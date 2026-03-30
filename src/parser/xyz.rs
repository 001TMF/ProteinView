use anyhow::{bail, Context, Result};
use crate::model::protein::{Protein, Chain, MoleculeType, Residue, Atom, SecondaryStructure};

/// Load a molecular structure from an XYZ file.
///
/// XYZ format:
///   Line 1: atom count (integer)
///   Line 2: comment / title
///   Lines 3+: Element  x  y  z
pub fn load_xyz(path: &str) -> Result<Protein> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read XYZ file: {}", path))?;

    let mut lines = content.lines();

    // Line 1: atom count
    let count_line = lines.next().context("XYZ file is empty")?;
    let atom_count: usize = count_line.trim().parse()
        .with_context(|| format!("First line must be an atom count, got: '{}'", count_line.trim()))?;

    // Line 2: comment / title
    let name = lines.next().unwrap_or("").trim().to_string();
    let name = if name.is_empty() {
        std::path::Path::new(path)
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "Unknown".to_string())
    } else {
        name
    };

    // Atom lines
    let mut atoms = Vec::with_capacity(atom_count);
    for (i, line) in lines.enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            bail!("Line {}: expected 'Element x y z', got: '{}'", i + 3, line);
        }
        let element = parts[0].to_string();
        let x: f64 = parts[1].parse()
            .with_context(|| format!("Line {}: invalid x coordinate '{}'", i + 3, parts[1]))?;
        let y: f64 = parts[2].parse()
            .with_context(|| format!("Line {}: invalid y coordinate '{}'", i + 3, parts[2]))?;
        let z: f64 = parts[3].parse()
            .with_context(|| format!("Line {}: invalid z coordinate '{}'", i + 3, parts[3]))?;

        atoms.push(Atom {
            name: element.clone(),
            element,
            x,
            y,
            z,
            b_factor: 0.0,
            is_backbone: false,
            is_hetero: false,
        });
    }

    if atoms.len() != atom_count {
        bail!(
            "XYZ header declares {} atoms but file contains {}",
            atom_count,
            atoms.len()
        );
    }

    let residue = Residue {
        name: "MOL".to_string(),
        seq_num: 1,
        atoms,
        secondary_structure: SecondaryStructure::Coil,
    };

    let chain = Chain {
        id: "A".to_string(),
        residues: vec![residue],
        molecule_type: MoleculeType::SmallMolecule,
    };

    Ok(Protein {
        name,
        chains: vec![chain],
        ligands: vec![],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp_xyz(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    #[test]
    fn test_load_water() {
        let f = write_temp_xyz(
            "3\nWater molecule\nO  0.000  0.000  0.117\nH  0.000  0.757 -0.469\nH  0.000 -0.757 -0.469\n",
        );
        let protein = load_xyz(f.path().to_str().unwrap()).unwrap();
        assert_eq!(protein.name, "Water molecule");
        assert_eq!(protein.chains.len(), 1);
        assert_eq!(protein.chains[0].molecule_type, MoleculeType::SmallMolecule);
        assert_eq!(protein.chains[0].residues.len(), 1);
        assert_eq!(protein.chains[0].residues[0].atoms.len(), 3);
        assert_eq!(protein.chains[0].residues[0].atoms[0].element, "O");
    }

    #[test]
    fn test_wrong_atom_count() {
        let f = write_temp_xyz("5\nBad count\nO 0 0 0\nH 1 0 0\n");
        let result = load_xyz(f.path().to_str().unwrap());
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("declares 5 atoms but file contains 2"));
    }

    #[test]
    fn test_bad_coordinate() {
        let f = write_temp_xyz("1\nBad\nO 0 abc 0\n");
        let result = load_xyz(f.path().to_str().unwrap());
        assert!(result.is_err());
    }

    #[test]
    fn test_empty_comment_uses_filename() {
        let f = write_temp_xyz("1\n\nC 0 0 0\n");
        let protein = load_xyz(f.path().to_str().unwrap()).unwrap();
        // Name should be derived from temp file name, not empty
        assert!(!protein.name.is_empty());
    }
}

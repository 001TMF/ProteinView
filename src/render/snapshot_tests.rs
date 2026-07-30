use std::fs;

use image::GenericImageView;
use tempfile::tempdir;

use super::*;
use crate::model::protein::{Atom, Chain, MoleculeType, Protein, Residue, SecondaryStructure};
use crate::model::selection::ResidueColorOverrides;

fn fixture_protein() -> Protein {
    let atom = |name: &str, x: f64, y: f64| Atom {
        name: name.to_string(),
        element: "C".to_string(),
        x,
        y,
        z: 0.0,
        b_factor: 20.0,
        is_backbone: true,
        is_hetero: false,
    };
    Protein {
        name: "fixture".to_string(),
        chains: vec![Chain {
            id: "A".to_string(),
            residues: vec![
                Residue {
                    name: "ALA".to_string(),
                    seq_num: 1,
                    insertion_code: None,
                    atoms: vec![atom("CA", -5.0, -2.0)],
                    secondary_structure: SecondaryStructure::Coil,
                },
                Residue {
                    name: "GLY".to_string(),
                    seq_num: 2,
                    insertion_code: None,
                    atoms: vec![atom("CA", 5.0, 2.0)],
                    secondary_structure: SecondaryStructure::Coil,
                },
            ],
            molecule_type: MoleculeType::Protein,
        }],
        ligands: Vec::new(),
    }
}

#[test]
fn writes_a_fullhd_png_without_terminal_setup() {
    let dir = tempdir().unwrap();
    let output = dir.path().join("protein.png");
    save_png(
        fixture_protein(),
        &output,
        SnapshotOptions {
            width: 320,
            height: 180,
            color_override: None,
            residue_colors: ResidueColorOverrides::default(),
            viz_mode: VizMode::Backbone,
            user_explicit_mode: true,
            show_ligands: true,
            interface_chain: None,
            show_interactions: false,
        },
    )
    .unwrap();

    let bytes = fs::read(&output).unwrap();
    assert!(bytes.starts_with(b"\x89PNG\r\n\x1a\n"));
    assert_eq!(
        image::load_from_memory(&bytes).unwrap().dimensions(),
        (320, 180)
    );
}

#[test]
fn rejects_oversized_snapshot_dimensions() {
    let dir = tempdir().unwrap();
    let error = save_png(
        fixture_protein(),
        &dir.path().join("oversized.png"),
        SnapshotOptions {
            width: 4096,
            height: 4096,
            color_override: None,
            residue_colors: ResidueColorOverrides::default(),
            viz_mode: VizMode::Backbone,
            user_explicit_mode: true,
            show_ligands: true,
            interface_chain: None,
            show_interactions: false,
        },
    )
    .unwrap_err();

    assert!(error.to_string().contains("pixel safety limit"));
}

#[test]
fn rejects_an_unknown_interface_chain() {
    let dir = tempdir().unwrap();
    let error = save_png(
        fixture_protein(),
        &dir.path().join("interface.png"),
        SnapshotOptions {
            width: 320,
            height: 180,
            color_override: None,
            residue_colors: ResidueColorOverrides::default(),
            viz_mode: VizMode::Backbone,
            user_explicit_mode: true,
            show_ligands: true,
            interface_chain: Some("Z".to_string()),
            show_interactions: true,
        },
    )
    .unwrap_err();

    assert!(error.to_string().contains("available chains: A"));
}

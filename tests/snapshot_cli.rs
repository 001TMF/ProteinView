use std::process::{Command, Stdio};

use image::GenericImageView;
use tempfile::tempdir;

#[test]
fn snapshot_cli_is_terminal_independent_and_writes_exact_dimensions() {
    let output_dir = tempdir().unwrap();
    let output_path = output_dir.path().join("1ubq-fullhd.png");
    let fixture = format!("{}/examples/1UBQ.pdb", env!("CARGO_MANIFEST_DIR"));

    let result = Command::new(env!("CARGO_BIN_EXE_proteinview"))
        .args([
            &fixture,
            "--snapshot",
            output_path.to_str().unwrap(),
            "--snapshot-width",
            "320",
            "--snapshot-height",
            "180",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "snapshot failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(!result.stdout.contains(&0x1b));
    assert!(!result.stderr.contains(&0x1b));

    let image = image::open(&output_path).unwrap();
    assert_eq!(image.dimensions(), (320, 180));
    assert!(image.to_rgba8().pixels().any(|pixel| pixel[3] > 0));
}

#[test]
fn snapshot_cli_accepts_repeatable_exact_residue_colors() {
    let output_dir = tempdir().unwrap();
    let output_path = output_dir.path().join("1ubq-residues.png");
    let fixture = format!("{}/examples/1UBQ.pdb", env!("CARGO_MANIFEST_DIR"));

    let result = Command::new(env!("CARGO_BIN_EXE_proteinview"))
        .args([
            &fixture,
            "--snapshot",
            output_path.to_str().unwrap(),
            "--snapshot-width",
            "320",
            "--snapshot-height",
            "180",
            "--mode",
            "backbone",
            "--residue-color",
            "A:1=FF0000",
            "--residue-color",
            "A:2=00FF00",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "snapshot failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(image::open(output_path).unwrap().dimensions(), (320, 180));
}

#[test]
fn snapshot_cli_rejects_unknown_residue_before_writing_output() {
    let output_dir = tempdir().unwrap();
    let output_path = output_dir.path().join("unknown.png");
    let fixture = format!("{}/examples/1UBQ.pdb", env!("CARGO_MANIFEST_DIR"));

    let result = Command::new(env!("CARGO_BIN_EXE_proteinview"))
        .args([
            &fixture,
            "--snapshot",
            output_path.to_str().unwrap(),
            "--residue-color",
            "Z:999=FF0000",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();

    assert!(!result.status.success());
    assert!(!output_path.exists());
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("available chains"),
        "unexpected error: {}",
        String::from_utf8_lossy(&result.stderr)
    );
}

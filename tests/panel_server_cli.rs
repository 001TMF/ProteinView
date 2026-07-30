use std::io::Write;
use std::process::{Command, Stdio};

use image::GenericImageView;
use serde_json::Value;
use tempfile::tempdir;

#[test]
fn panel_server_is_terminal_independent_and_persists_until_shutdown() {
    let output_dir = tempdir().unwrap();
    let output_path = output_dir.path().join("live-panel.png");
    let fixture = format!("{}/examples/1UBQ.pdb", env!("CARGO_MANIFEST_DIR"));
    let mut child = Command::new(env!("CARGO_BIN_EXE_proteinview"))
        .args([
            &fixture,
            "--panel-server",
            "--output",
            output_path.to_str().unwrap(),
            "--panel-width",
            "160",
            "--panel-height",
            "96",
            "--mode",
            "backbone",
            "--residue-color",
            "A:1=FF0000",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();

    child
        .stdin
        .take()
        .unwrap()
        .write_all(
            br#"{"id":1,"command":"rotate","axis":"y","delta":0.2}
{"id":2,"command":"resize","width":320,"height":180}
{"id":3,"command":"get_state"}
{"id":4,"command":"shutdown"}
"#,
        )
        .unwrap();
    let result = child.wait_with_output().unwrap();

    assert!(
        result.status.success(),
        "panel server failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(!result.stdout.contains(&0x1b));
    assert!(!result.stderr.contains(&0x1b));

    let records = String::from_utf8(result.stdout)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .collect::<Vec<_>>();
    assert_eq!(records[0]["type"], "ready");
    assert_eq!(records[0]["revision"], 1);
    assert_eq!(records[1]["revision"], 2);
    assert_eq!(records[2]["revision"], 3);
    assert_eq!(records[3]["revision"], 3);
    assert_eq!(records[4]["revision"], 3);
    assert_eq!(records[3]["state"]["camera"]["rot_y"], 0.2);
    assert_eq!(
        records[3]["state"]["presentation"]["residue_colors"][0]["chain"],
        "A"
    );
    assert_eq!(
        records[3]["state"]["presentation"]["residue_colors"][0]["residue_number"],
        1
    );
    assert_eq!(
        records[3]["state"]["presentation"]["residue_colors"][0]["color"],
        "FF0000"
    );

    let image = image::open(&output_path).unwrap();
    assert_eq!(image.dimensions(), (320, 180));
    assert!(image.to_rgba8().pixels().any(|pixel| pixel[3] > 0));
}

#[test]
fn xyz_distinguishes_default_element_from_explicit_structure_color() {
    let output_dir = tempdir().unwrap();
    let fixture = output_dir.path().join("molecule.xyz");
    std::fs::write(
        &fixture,
        "2\nminimal molecule\nC 0.0 0.0 0.0\nO 1.2 0.0 0.0\n",
    )
    .unwrap();

    for (label, explicit_color, expected_color) in [
        ("default", false, "element"),
        ("explicit", true, "structure"),
    ] {
        let output_path = output_dir.path().join(format!("{label}.png"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_proteinview"));
        command.args([
            fixture.to_str().unwrap(),
            "--panel-server",
            "--output",
            output_path.to_str().unwrap(),
            "--panel-width",
            "64",
            "--panel-height",
            "64",
        ]);
        if explicit_color {
            command.args(["--color", "structure"]);
        }
        let mut child = command
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child
            .stdin
            .take()
            .unwrap()
            .write_all(b"{\"id\":1,\"command\":\"shutdown\"}\n")
            .unwrap();
        let result = child.wait_with_output().unwrap();
        assert!(
            result.status.success(),
            "{label} XYZ panel failed: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let ready: Value = serde_json::from_str(
            String::from_utf8(result.stdout)
                .unwrap()
                .lines()
                .next()
                .unwrap(),
        )
        .unwrap();
        assert_eq!(
            ready["state"]["presentation"]["color"], expected_color,
            "{label} XYZ color"
        );
    }
}

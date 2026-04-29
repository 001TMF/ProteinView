//! Integration test: build the proteinview binary, invoke it with --batch,
//! verify the output directory is populated.

use std::process::Command;

#[test]
fn batch_mode_end_to_end() {
    let tmp = tempfile::tempdir().expect("tempdir");
    let out_dir = tmp.path().join("frames");
    std::fs::create_dir_all(&out_dir).expect("mkdir");

    // Write a tiny config pointing at the bundled example PDB
    let cfg_path = tmp.path().join("config.json");
    let cfg = format!(
        r#"{{
            "input": "examples/1UBQ.pdb",
            "output_dir": "{}",
            "frames": 3,
            "width": 200,
            "height": 150,
            "render_mode": "fullhd",
            "waypoints": [
                {{"t": 0.0, "rot_x": 0.0, "rot_y": 0.0, "rot_z": 0.0, "zoom": 1.0}},
                {{"t": 1.0, "rot_x": 0.0, "rot_y": 1.5, "rot_z": 0.0, "zoom": 1.0}}
            ]
        }}"#,
        out_dir.display()
    );
    std::fs::write(&cfg_path, cfg).expect("write config");

    let binary = env!("CARGO_BIN_EXE_proteinview");
    let status = Command::new(binary)
        .arg("--batch")
        .arg(cfg_path.to_str().unwrap())
        .status()
        .expect("run proteinview");
    assert!(status.success(), "binary exited non-zero");

    for i in 1..=3 {
        let p = out_dir.join(format!("frame_{:04}.png", i));
        assert!(p.exists(), "missing frame {}", p.display());
    }
}

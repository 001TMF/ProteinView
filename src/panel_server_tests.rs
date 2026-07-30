use std::io::Cursor;
use std::time::Duration;

use image::GenericImageView;
use serde_json::Value;
use tempfile::tempdir;

use super::*;
use crate::model::protein::{Atom, Chain, MoleculeType, Protein, Residue, SecondaryStructure};
use crate::model::selection::ResidueColorOverrides;

fn fixture_protein() -> Protein {
    let residue = |name: &str, seq_num: i32, x: f64, y: f64| Residue {
        name: name.to_string(),
        seq_num,
        insertion_code: None,
        atoms: vec![Atom {
            name: "CA".to_string(),
            element: "C".to_string(),
            x,
            y,
            z: 0.0,
            b_factor: 20.0,
            is_backbone: true,
            is_hetero: false,
        }],
        secondary_structure: SecondaryStructure::Coil,
    };
    Protein {
        name: "panel-fixture".to_string(),
        chains: vec![
            Chain {
                id: "A".to_string(),
                residues: vec![residue("ALA", 1, -4.0, -1.0), residue("GLY", 2, -1.0, 1.0)],
                molecule_type: MoleculeType::Protein,
            },
            Chain {
                id: "B".to_string(),
                residues: vec![residue("VAL", 1, 1.0, -1.0), residue("LEU", 2, 4.0, 1.0)],
                molecule_type: MoleculeType::Protein,
            },
        ],
        ligands: Vec::new(),
    }
}

fn options(width: u32, height: u32) -> PanelServerOptions {
    PanelServerOptions {
        width,
        height,
        color_override: None,
        residue_colors: ResidueColorOverrides::default(),
        viz_mode: VizMode::Backbone,
        user_explicit_mode: true,
    }
}

fn run_requests(requests: &str) -> (Vec<Value>, tempfile::TempDir, PathBuf) {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("panel.png");
    let mut output = Vec::new();
    serve(
        fixture_protein(),
        &output_path,
        options(160, 96),
        Cursor::new(requests.as_bytes()),
        &mut output,
    )
    .unwrap();
    let records = String::from_utf8(output)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    (records, dir, output_path)
}

fn response_with_id(records: &[Value], id: i64) -> &Value {
    records
        .iter()
        .find(|record| record.get("id").and_then(Value::as_i64) == Some(id))
        .unwrap()
}

#[test]
fn ready_state_bounds_untrusted_structure_name() {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("panel.png");
    let mut protein = fixture_protein();
    protein.name = format!("{}\n{}", "\u{1b}".repeat(32), "X".repeat(2_000));
    let mut output = Vec::new();

    serve(
        protein,
        &output_path,
        options(160, 96),
        Cursor::new(b"{\"id\":1,\"command\":\"shutdown\"}\n"),
        &mut output,
    )
    .unwrap();

    let ready: Value =
        serde_json::from_slice(output.split(|byte| *byte == b'\n').next().unwrap()).unwrap();
    let name = ready["state"]["structure"]["name"].as_str().unwrap();
    assert!(name.len() <= MAX_STATE_NAME_BYTES);
    assert!(!name.chars().any(char::is_control));
}

#[test]
fn oversized_required_chain_metadata_is_rejected_before_rendering() {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("panel.png");
    let mut protein = fixture_protein();
    protein.chains = (0..2_200)
        .map(|index| Chain {
            id: format!("CHAIN-{index:026}"),
            residues: Vec::new(),
            molecule_type: MoleculeType::Protein,
        })
        .collect();
    let mut output = Vec::new();

    let error = serve(
        protein,
        &output_path,
        options(160, 96),
        Cursor::new([]),
        &mut output,
    )
    .unwrap_err();

    assert!(error.to_string().contains("panel response exceeds"));
    assert!(output.is_empty());
    assert!(!output_path.exists());
}

#[test]
fn near_limit_ids_and_commands_return_bounded_recoverable_errors() {
    let oversized_id = "i".repeat(60_000);
    let unknown_command = "x".repeat(60_000);
    let requests = format!(
        "{}\n{}\n{}\n{}\n",
        json!({"id": oversized_id, "command": "get_state"}),
        json!({"id": 7, "command": unknown_command}),
        json!({"id": 8, "command": "get_state"}),
        json!({"id": 9, "command": "shutdown"}),
    );
    let (records, _dir, _output_path) = run_requests(&requests);

    assert_eq!(records[1]["error"]["code"], "invalid_request");
    assert_eq!(records[1]["id"], Value::Null);
    assert_eq!(records[2]["error"]["code"], "unknown_command");
    assert!(records[2]["error"]["message"].as_str().unwrap().len() <= MAX_ERROR_MESSAGE_BYTES);
    assert_eq!(response_with_id(&records, 8)["ok"], true);
    assert!(
        records
            .iter()
            .all(|record| serde_json::to_vec(record).unwrap().len() < MAX_RESPONSE_BYTES)
    );
}

#[test]
fn oversized_residue_state_is_rejected_before_replacing_the_frame() {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("panel.png");
    let mut protein = fixture_protein();
    let residue_template = protein.chains[0].residues[0].clone();
    protein.chains = (0..1_000)
        .map(|chain_index| {
            let id = format!("CHAIN-{chain_index:026}");
            let residues = if chain_index == 0 {
                (0..MAX_RESIDUE_COLORS)
                    .map(|residue_number| {
                        let mut residue = residue_template.clone();
                        residue.seq_num = residue_number as i32;
                        residue.name = "LONG-MODIFIED-RESIDUE-NAME".to_string();
                        residue
                    })
                    .collect()
            } else {
                Vec::new()
            };
            Chain {
                id,
                residues,
                molecule_type: MoleculeType::Protein,
            }
        })
        .collect();
    let chain = protein.chains[0].id.clone();
    let residues = (0..MAX_RESIDUE_COLORS)
        .map(|residue_number| {
            json!({
                "chain": chain,
                "residue_number": residue_number,
                "color": "FF0000",
            })
        })
        .collect::<Vec<_>>();
    let requests = format!(
        "{}\n{}\n{}\n",
        json!({"id": 1, "command": "set_residue_colors", "residues": residues}),
        json!({"id": 2, "command": "get_state"}),
        json!({"id": 3, "command": "shutdown"}),
    );
    assert!(requests.lines().all(|line| line.len() <= MAX_REQUEST_BYTES));
    let mut output = Vec::new();

    serve(
        protein,
        &output_path,
        options(160, 96),
        Cursor::new(requests),
        &mut output,
    )
    .unwrap();
    let records = String::from_utf8(output)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .collect::<Vec<_>>();

    assert_eq!(
        response_with_id(&records, 1)["error"]["code"],
        "response_too_large"
    );
    assert_eq!(response_with_id(&records, 1)["revision"], 1);
    assert_eq!(response_with_id(&records, 2)["revision"], 1);
    assert_eq!(
        response_with_id(&records, 2)["state"]["presentation"]["residue_colors"],
        json!([])
    );
}

#[test]
fn protocol_persists_camera_state_and_revisions() {
    let (records, _dir, output_path) = run_requests(
        r#"{"id":1,"command":"rotate","axis":"x","delta":0.25}
{"id":2,"command":"pan","dx":4.0,"dy":-3.0}
{"id":3,"command":"zoom","factor":1.5}
{"id":4,"command":"get_state"}
{"id":5,"command":"shutdown"}
"#,
    );

    assert_eq!(records[0]["type"], "ready");
    assert_eq!(records[0]["revision"], 1);
    assert_eq!(response_with_id(&records, 1)["revision"], 2);
    assert_eq!(response_with_id(&records, 2)["revision"], 3);
    assert_eq!(response_with_id(&records, 3)["revision"], 4);
    assert_eq!(response_with_id(&records, 4)["revision"], 4);
    assert_eq!(response_with_id(&records, 5)["revision"], 4);

    let state = &response_with_id(&records, 4)["state"];
    assert_eq!(state["camera"]["rot_x"], 0.25);
    assert_eq!(state["camera"]["pan_x"], 4.0);
    assert_eq!(state["camera"]["pan_y"], -3.0);
    assert!(state["camera"]["zoom"].as_f64().unwrap() > 1.0);
    assert_eq!(image::open(output_path).unwrap().dimensions(), (160, 96));
}

#[test]
fn presentation_commands_update_authoritative_state() {
    let (records, _dir, _output_path) = run_requests(
        r#"{"id":1,"command":"set_color","color":"chain"}
{"id":2,"command":"set_viz","mode":"wireframe"}
{"id":3,"command":"select_chain","direction":"next"}
{"id":4,"command":"set_interface","enabled":true}
{"id":5,"command":"set_interactions","enabled":true}
{"id":6,"command":"set_ligands","enabled":false}
{"id":7,"command":"get_state"}
{"id":8,"command":"shutdown"}
"#,
    );

    let state = &response_with_id(&records, 7)["state"];
    assert_eq!(state["presentation"]["color"], "chain");
    assert_eq!(state["presentation"]["effective_color"], "interface");
    assert_eq!(state["presentation"]["viz_mode"], "wireframe");
    assert_eq!(state["presentation"]["current_chain_index"], 1);
    assert_eq!(state["presentation"]["current_chain_id"], "B");
    assert_eq!(state["presentation"]["interface"], true);
    assert_eq!(state["presentation"]["interactions"], true);
    assert_eq!(state["presentation"]["ligands"], false);
}

#[test]
fn set_chain_selects_one_exact_chain_in_one_render() {
    let (records, _dir, _output_path) = run_requests(
        r#"{"id":1,"command":"set_chain","chain":"B"}
{"id":2,"command":"get_state"}
{"id":3,"command":"set_chain","chain":"missing"}
{"id":4,"command":"get_state"}
{"id":5,"command":"shutdown"}
"#,
    );

    assert_eq!(response_with_id(&records, 1)["revision"], 2);
    assert_eq!(
        response_with_id(&records, 2)["state"]["presentation"]["current_chain_id"],
        "B"
    );
    assert_eq!(response_with_id(&records, 3)["ok"], false);
    assert_eq!(
        response_with_id(&records, 3)["error"]["code"],
        "invalid_params"
    );
    assert_eq!(response_with_id(&records, 4)["revision"], 2);
    assert_eq!(
        response_with_id(&records, 4)["state"]["presentation"]["current_chain_id"],
        "B"
    );
}

#[test]
fn residue_colors_update_atomically_and_are_normalized_in_state() {
    let (records, _dir, _output_path) = run_requests(
        r#"{"id":1,"command":"set_residue_colors","residues":[{"chain":"B","residue_number":2,"color":"00FF7F"},{"chain":"A","residue_number":1,"insertion_code":null,"color":"FF0000"}]}
{"id":2,"command":"get_state"}
{"id":3,"command":"set_residue_colors","residues":[]}
{"id":4,"command":"get_state"}
{"id":5,"command":"shutdown"}
"#,
    );

    assert_eq!(response_with_id(&records, 1)["revision"], 2);
    assert_eq!(
        response_with_id(&records, 2)["state"]["presentation"]["residue_colors"],
        json!([
            {
                "chain": "B",
                "residue_number": 2,
                "insertion_code": null,
                "residue_name": "LEU",
                "color": "00FF7F",
            },
            {
                "chain": "A",
                "residue_number": 1,
                "insertion_code": null,
                "residue_name": "ALA",
                "color": "FF0000",
            }
        ])
    );
    assert_eq!(response_with_id(&records, 3)["revision"], 3);
    assert_eq!(
        response_with_id(&records, 4)["state"]["presentation"]["residue_colors"],
        json!([])
    );
}

#[test]
fn invalid_residue_color_batch_preserves_previous_state_and_revision() {
    let (records, _dir, _output_path) = run_requests(
        r#"{"id":1,"command":"set_residue_colors","residues":[{"chain":"A","residue_number":1,"color":"FF0000"}]}
{"id":2,"command":"set_residue_colors","residues":[{"chain":"A","residue_number":1,"color":"00ff00"}]}
{"id":3,"command":"set_residue_colors","residues":[{"chain":"A","residue_number":1,"color":"00FF00"},{"chain":"Z","residue_number":9,"color":"0000FF"}]}
{"id":4,"command":"get_state"}
{"id":5,"command":"shutdown"}
"#,
    );

    assert_eq!(response_with_id(&records, 1)["revision"], 2);
    assert_eq!(response_with_id(&records, 2)["ok"], false);
    assert_eq!(
        response_with_id(&records, 2)["error"]["code"],
        "invalid_params"
    );
    assert_eq!(response_with_id(&records, 2)["revision"], 2);
    assert_eq!(response_with_id(&records, 3)["ok"], false);
    assert_eq!(
        response_with_id(&records, 3)["error"]["code"],
        "invalid_params"
    );
    assert_eq!(response_with_id(&records, 3)["revision"], 2);
    assert_eq!(
        response_with_id(&records, 4)["state"]["presentation"]["residue_colors"],
        json!([{
            "chain": "A",
            "residue_number": 1,
            "insertion_code": null,
            "residue_name": "ALA",
            "color": "FF0000",
        }])
    );
}

#[test]
fn invalid_requests_are_recoverable_and_do_not_advance_revision() {
    let (records, _dir, _output_path) = run_requests(
        r#"
{
{"id":1,"command":"not_a_command"}
{"id":2,"command":"resize","width":4096,"height":4096}
{"id":3,"command":"zoom","factor":0.0}
{"id":4,"command":"set_interactions","enabled":true}
{"id":5,"command":"get_state"}
{"id":6,"command":"shutdown"}
"#,
    );

    let error_codes = records
        .iter()
        .filter_map(|record| record.pointer("/error/code").and_then(Value::as_str))
        .collect::<Vec<_>>();
    assert_eq!(
        error_codes,
        vec![
            "invalid_json",
            "unknown_command",
            "invalid_params",
            "invalid_params",
            "invalid_state",
        ]
    );
    assert_eq!(response_with_id(&records, 5)["revision"], 1);
    assert_eq!(
        response_with_id(&records, 5)["state"]["viewport"],
        json!({"width": 160, "height": 96})
    );
}

#[test]
fn camera_overflow_is_rejected_without_changing_state_or_revision() {
    let (records, _dir, _output_path) = run_requests(
        r#"{"id":1,"command":"rotate","axis":"x","delta":1e308}
{"id":2,"command":"rotate","axis":"x","delta":1e308}
{"id":3,"command":"pan","dx":1e308,"dy":-1e308}
{"id":4,"command":"pan","dx":1e308,"dy":-1e308}
{"id":5,"command":"get_state"}
{"id":6,"command":"shutdown"}
"#,
    );

    for id in [2, 4] {
        let response = response_with_id(&records, id);
        assert_eq!(response["ok"], false, "command {id} should reject overflow");
        assert_eq!(
            response["error"]["code"], "invalid_params",
            "command {id} should report invalid_params"
        );
    }

    let state = &response_with_id(&records, 5)["state"];
    assert_eq!(response_with_id(&records, 5)["revision"], 3);
    assert_eq!(state["camera"]["rot_x"], 1e308);
    assert_eq!(state["camera"]["pan_x"], 1e308);
    assert_eq!(state["camera"]["pan_y"], -1e308);

    let dir = tempdir().unwrap();
    let output_path = dir.path().join("panel.png");
    let mut session = PanelSession::new(fixture_protein(), &output_path, options(160, 96)).unwrap();
    session.settings.camera.zoom = f64::MAX;
    let error = apply_command(
        &mut session,
        PanelCommand::Zoom {
            factor: None,
            direction: Some(ZoomDirection::In),
        },
    )
    .err()
    .unwrap();
    assert_eq!(error.code, "invalid_params");
    assert_eq!(session.settings.camera.zoom, f64::MAX);

    session.settings.camera.zoom = f64::from_bits(1);
    let error = apply_command(
        &mut session,
        PanelCommand::Zoom {
            factor: Some(f64::from_bits(1)),
            direction: None,
        },
    )
    .err()
    .unwrap();
    assert_eq!(error.code, "invalid_params");
    assert_eq!(session.settings.camera.zoom, f64::from_bits(1));
}

#[test]
fn resize_atomically_writes_exact_png_dimensions() {
    let (records, _dir, output_path) = run_requests(
        r#"{"id":1,"command":"resize","width":320,"height":180}
{"id":2,"command":"shutdown"}
"#,
    );

    let bytes = std::fs::read(&output_path).unwrap();
    assert!(bytes.starts_with(b"\x89PNG\r\n\x1a\n"));
    assert_eq!(
        image::load_from_memory(&bytes).unwrap().dimensions(),
        (320, 180)
    );
    assert_eq!(
        response_with_id(&records, 1)["frame"],
        json!({
            "path": output_path.to_string_lossy(),
            "mime_type": "image/png",
            "width": 320,
            "height": 180,
        })
    );
}

#[test]
fn reset_reproduces_the_initial_camera_frame() {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("panel.png");
    let mut session = PanelSession::new(fixture_protein(), &output_path, options(160, 96)).unwrap();
    session.render_next_revision().unwrap();
    let baseline = std::fs::read(&output_path).unwrap();

    apply_command(
        &mut session,
        PanelCommand::Rotate {
            axis: Axis::Y,
            delta: 0.5,
        },
    )
    .unwrap();
    session.render_next_revision().unwrap();
    let rotated = std::fs::read(&output_path).unwrap();
    assert_ne!(rotated, baseline);

    apply_command(&mut session, PanelCommand::Reset).unwrap();
    session.render_next_revision().unwrap();
    assert_eq!(std::fs::read(output_path).unwrap(), baseline);
}

#[test]
fn capped_reader_discards_an_oversized_line_and_recovers() {
    let mut input = vec![b'x'; MAX_REQUEST_BYTES + 1];
    input.extend_from_slice(b"\n{\"id\":1,\"command\":\"get_state\"}\n");
    let mut reader = Cursor::new(input);

    assert!(matches!(
        read_capped_line(&mut reader).unwrap(),
        CappedLine::TooLarge
    ));
    match read_capped_line(&mut reader).unwrap() {
        CappedLine::Line(line) => {
            assert_eq!(
                serde_json::from_slice::<Value>(&line).unwrap()["command"],
                "get_state"
            );
        }
        CappedLine::Eof | CappedLine::TooLarge => panic!("expected the next valid request"),
    }
}

#[test]
fn floating_point_request_ids_are_rejected() {
    let (records, _dir, _output_path) = run_requests(
        r#"{"id":1.5,"command":"get_state"}
{"id":2,"command":"shutdown"}
"#,
    );

    assert_eq!(records[1]["ok"], false);
    assert_eq!(records[1]["error"]["code"], "invalid_request");
    assert_eq!(records[1]["revision"], 1);
    assert_eq!(response_with_id(&records, 2)["ok"], true);
}

#[test]
fn render_advances_an_enabled_auto_rotate_camera() {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("panel.png");
    let mut session = PanelSession::new(fixture_protein(), &output_path, options(160, 96)).unwrap();
    session.render_next_revision().unwrap();
    apply_command(&mut session, PanelCommand::SetAutoRotate { enabled: true }).unwrap();
    session.render_next_revision().unwrap();
    let before = session.settings.camera.rot_y;

    std::thread::sleep(Duration::from_millis(5));
    apply_command(&mut session, PanelCommand::Render).unwrap();
    session.render_next_revision().unwrap();

    assert!(session.settings.camera.rot_y < before);
}

#[test]
fn dimension_limits_match_snapshot_safety_caps() {
    validate_dimensions(64, 64).unwrap();
    validate_dimensions(4096, 2048).unwrap();
    assert!(validate_dimensions(63, 64).is_err());
    assert!(validate_dimensions(4097, 64).is_err());
    assert!(validate_dimensions(4096, 2049).is_err());
}

#[test]
fn capped_reader_handles_an_oversized_final_line_at_eof() {
    let mut reader = Cursor::new(vec![b'x'; MAX_REQUEST_BYTES + 1]);

    assert!(matches!(
        read_capped_line(&mut reader).unwrap(),
        CappedLine::TooLarge
    ));
    assert!(matches!(
        read_capped_line(&mut reader).unwrap(),
        CappedLine::Eof
    ));
}

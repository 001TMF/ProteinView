//! Batch (headless) rendering mode for ProteinView.
//!
//! Reads a JSON config describing input PDB, camera waypoints, color scheme,
//! hotspot residues, and output dir; renders a PNG sequence without entering
//! the TUI.

use std::collections::HashMap;

use image::RgbImage;
use serde::{Deserialize, Serialize};

/// Top-level batch render configuration.
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct BatchConfig {
    /// Path to input PDB / mmCIF / XYZ file.
    pub input: String,

    /// Directory where `frame_NNNN.png` sequence will be written.
    pub output_dir: String,

    /// Total number of frames to render.
    pub frames: u32,

    /// Output image dimensions (pixels).
    pub width: u32,
    pub height: u32,

    /// Render mode: "fullhd", "hd", or "braille".  Batch mode always
    /// produces RGB PNGs regardless — this string controls the rendering
    /// path used internally (FullHD = pixel-perfect rasterization).
    #[serde(default = "default_render_mode")]
    pub render_mode: String,

    /// Color scheme: "structure", "element", "chain", "plddt", "bfactor", "rainbow".
    #[serde(default = "default_color")]
    pub color: String,

    /// Visualization mode: "cartoon", "backbone", "wireframe".
    #[serde(default = "default_viz")]
    pub viz: String,

    /// Camera waypoints — at least 1.  Linearly interpolated across frames.
    pub waypoints: Vec<Waypoint>,

    /// Optional hotspot residues to highlight.
    #[serde(default)]
    pub hotspots: Vec<HotspotSpec>,
}

/// One camera position in the timeline.  `t` is a normalized time in [0.0, 1.0]
/// over the full frame range.
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Waypoint {
    pub t: f64,
    pub rot_x: f64,
    pub rot_y: f64,
    pub rot_z: f64,
    pub zoom: f64,
    #[serde(default)]
    pub pan_x: f64,
    #[serde(default)]
    pub pan_y: f64,
}

/// One residue (or set of residues) to highlight in the render.
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct HotspotSpec {
    /// Chain ID, e.g. "A".
    pub chain: String,
    /// 1-based residue numbers per the PDB.
    pub residues: Vec<i32>,
    /// RGB highlight color, 0..255 each.
    pub color: [u8; 3],
}

use crate::model::protein::Protein;
use crate::render::camera::Camera;

/// Returns an interpolated `Camera` at normalized time `t` along the waypoint
/// timeline.
///
/// Waypoints are sorted by `t` before evaluation. `t` is clamped to [first, last].
/// Interpolation is linear between consecutive waypoints.
pub fn camera_at(waypoints: &[Waypoint], t: f64) -> Camera {
    assert!(!waypoints.is_empty(), "waypoints must be non-empty");
    let mut sorted: Vec<&Waypoint> = waypoints.iter().collect();
    sorted.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));

    let first = sorted.first().unwrap();
    let last = sorted.last().unwrap();
    if t <= first.t {
        return waypoint_to_camera(first);
    }
    if t >= last.t {
        return waypoint_to_camera(last);
    }

    // Find bracketing waypoints
    let mut prev = sorted[0];
    let mut next = sorted[0];
    for w in sorted.windows(2) {
        if w[0].t <= t && t <= w[1].t {
            prev = w[0];
            next = w[1];
            break;
        }
    }

    let span = next.t - prev.t;
    let alpha = if span.abs() < 1e-12 { 0.0 } else { (t - prev.t) / span };
    let lerp = |a: f64, b: f64| a + alpha * (b - a);

    let mut cam = Camera::default();
    cam.rot_x = lerp(prev.rot_x, next.rot_x);
    cam.rot_y = lerp(prev.rot_y, next.rot_y);
    cam.rot_z = lerp(prev.rot_z, next.rot_z);
    cam.zoom = lerp(prev.zoom, next.zoom);
    cam.pan_x = lerp(prev.pan_x, next.pan_x);
    cam.pan_y = lerp(prev.pan_y, next.pan_y);
    cam
}

fn waypoint_to_camera(w: &Waypoint) -> Camera {
    let mut cam = Camera::default();
    cam.rot_x = w.rot_x;
    cam.rot_y = w.rot_y;
    cam.rot_z = w.rot_z;
    cam.zoom = w.zoom;
    cam.pan_x = w.pan_x;
    cam.pan_y = w.pan_y;
    cam
}

/// Render one frame to an RGB image.
///
/// Reuses the existing `render::draw_protein` → `render::hd::render_hd_framebuffer`
/// rasterization path so batch output is pixel-identical to what the FullHD
/// TUI mode would produce for the same camera and color settings.
///
/// Hotspot overrides in `cfg.hotspots` are currently stored in a `HighlightMap`
/// (chain_id → residue_seq_nums → RGB).  Full integration with the color scheme
/// (overriding individual atom colors) requires a more invasive change to the
/// `ColorScheme` type; for the MVP batch pipeline the map is computed here and
/// can be plumbed in future iterations (Task 6+).
pub fn render_frame(
    protein: &Protein,
    camera: &Camera,
    cfg: &BatchConfig,
) -> anyhow::Result<RgbImage> {
    use crate::render::color::ColorScheme;

    // Center the protein at the origin, mirroring App::new().
    // We clone here so the caller's copy is unchanged.
    let mut centered = protein.clone();
    centered.center();

    let color_type = crate::render::color::parse_color_scheme(&cfg.color);
    let viz_mode = crate::app::VizMode::parse(&cfg.viz);

    let total_residues = centered.residue_count();

    // Build highlight map and convert [u8; 3] → ratatui Color::Rgb for the
    // ColorScheme override table.
    let highlight_map = build_highlight_map(&cfg.hotspots);
    let highlight_colors: std::collections::HashMap<(String, i32), ratatui::style::Color> =
        highlight_map
            .into_iter()
            .map(|(k, [r, g, b])| (k, ratatui::style::Color::Rgb(r, g, b)))
            .collect();

    let color_scheme = ColorScheme::new(color_type, total_residues)
        .with_highlights(highlight_colors);

    // Auto-fit zoom: always compute a base zoom so the protein fills the
    // output frame (mirroring App::new(): base = 0.9 * min(w,h) / (2*r)).
    // The waypoint's `zoom` field is then applied as a **multiplier** on top
    // of that base — `zoom: 1.0` means "fill the frame exactly"; `zoom: 2.0`
    // means "2× zoomed in relative to auto-fit".  This eliminates the previous
    // float-equality sentinel (`cam.zoom == 1.0`) which would misfire on
    // interpolated values near 1.0 and gave no intuitive way to zoom out.
    let radius = centered.bounding_radius().max(1.0);
    let base_zoom = 0.9 * (cfg.width as f64).min(cfg.height as f64) / (2.0 * radius);
    let mut cam = camera.clone();
    cam.zoom = base_zoom * camera.zoom;

    let fb = crate::render::draw_protein(
        &centered,
        &cam,
        cfg.width as usize,
        cfg.height as usize,
        viz_mode,
        &color_scheme,
    )?;

    Ok(fb.to_rgb_image())
}

/// A per-residue color override table for hotspot highlighting.
///
/// Key: `(chain_id, residue_seq_num)` → RGB override color.
/// Reserved for future use when `ColorScheme` gains per-residue override hooks.
pub type HighlightMap = HashMap<(String, i32), [u8; 3]>;

/// Build a [`HighlightMap`] from a list of [`HotspotSpec`]s.
///
/// Note: not yet consumed by `draw_protein`; wired in Task 6.
pub fn build_highlight_map(hotspots: &[HotspotSpec]) -> HighlightMap {
    let mut map = HighlightMap::new();
    for hs in hotspots {
        for &seq_num in &hs.residues {
            map.insert((hs.chain.clone(), seq_num), hs.color);
        }
    }
    map
}

/// Top-level entrypoint: load the protein, then for each frame in [1, frames],
/// compute t = (i - 1) / max(1, frames - 1), interpolate camera, render,
/// and write `frame_NNNN.png` to `output_dir`.
pub fn run(cfg: &BatchConfig) -> anyhow::Result<()> {
    std::fs::create_dir_all(&cfg.output_dir)?;

    let protein = if cfg.input.to_lowercase().ends_with(".xyz") {
        crate::parser::xyz::load_xyz(&cfg.input)?
    } else {
        crate::parser::pdb::load_structure(&cfg.input)?
    };

    let denom = (cfg.frames.max(1) - 1).max(1) as f64;
    for i in 1..=cfg.frames {
        let t = (i.saturating_sub(1)) as f64 / denom;
        let cam = camera_at(&cfg.waypoints, t);
        let img = render_frame(&protein, &cam, cfg)?;
        let path = format!("{}/frame_{:04}.png", cfg.output_dir, i);
        img.save(&path).map_err(|e| anyhow::anyhow!("save {}: {}", path, e))?;
    }
    Ok(())
}

/// Default render mode when not specified in JSON.
fn default_render_mode() -> String {
    "fullhd".to_string()
}
/// Default color scheme when not specified in JSON.
fn default_color() -> String {
    "structure".to_string()
}
/// Default visualization mode when not specified in JSON.
fn default_viz() -> String {
    "cartoon".to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deserializes_minimal_config() {
        let json = r#"{
            "input": "examples/1AOI.pdb",
            "output_dir": "/tmp/render",
            "frames": 60,
            "width": 1920,
            "height": 1080,
            "waypoints": [
                {"t": 0.0, "rot_x": 0.0, "rot_y": 0.0, "rot_z": 0.0, "zoom": 1.0},
                {"t": 1.0, "rot_x": 0.0, "rot_y": 6.283, "rot_z": 0.0, "zoom": 1.0}
            ]
        }"#;
        let cfg: BatchConfig = serde_json::from_str(json).expect("must parse");
        assert_eq!(cfg.frames, 60);
        assert_eq!(cfg.render_mode, "fullhd");
        assert_eq!(cfg.color, "structure");
        assert_eq!(cfg.viz, "cartoon");
        assert!(cfg.hotspots.is_empty());
    }

    #[test]
    fn deserializes_with_hotspots() {
        let json = r#"{
            "input": "examples/1AOI.pdb",
            "output_dir": "/tmp/render",
            "frames": 1,
            "width": 1920,
            "height": 1080,
            "waypoints": [{"t": 0.0, "rot_x": 0.0, "rot_y": 0.0, "rot_z": 0.0, "zoom": 1.0}],
            "hotspots": [{"chain": "A", "residues": [21, 46, 48], "color": [255, 64, 0]}]
        }"#;
        let cfg: BatchConfig = serde_json::from_str(json).expect("must parse");
        assert_eq!(cfg.hotspots.len(), 1);
        assert_eq!(cfg.hotspots[0].chain, "A");
        assert_eq!(cfg.hotspots[0].residues, vec![21, 46, 48]);
    }

    #[test]
    fn interpolates_camera_at_t_zero() {
        let waypoints = vec![
            Waypoint { t: 0.0, rot_x: 0.0, rot_y: 0.0, rot_z: 0.0, zoom: 1.0, pan_x: 0.0, pan_y: 0.0 },
            Waypoint { t: 1.0, rot_x: 0.0, rot_y: std::f64::consts::TAU, rot_z: 0.0, zoom: 1.0, pan_x: 0.0, pan_y: 0.0 },
        ];
        let cam = camera_at(&waypoints, 0.0);
        assert!((cam.rot_y - 0.0).abs() < 1e-9);
        assert!((cam.zoom - 1.0).abs() < 1e-9);
    }

    #[test]
    fn interpolates_camera_midway() {
        let waypoints = vec![
            Waypoint { t: 0.0, rot_x: 0.0, rot_y: 0.0, rot_z: 0.0, zoom: 1.0, pan_x: 0.0, pan_y: 0.0 },
            Waypoint { t: 1.0, rot_x: 0.0, rot_y: 2.0, rot_z: 0.0, zoom: 2.0, pan_x: 0.0, pan_y: 0.0 },
        ];
        let cam = camera_at(&waypoints, 0.5);
        assert!((cam.rot_y - 1.0).abs() < 1e-9);
        assert!((cam.zoom - 1.5).abs() < 1e-9);
    }

    #[test]
    fn clamps_outside_range() {
        let waypoints = vec![
            Waypoint { t: 0.0, rot_x: 0.0, rot_y: 0.0, rot_z: 0.0, zoom: 1.0, pan_x: 0.0, pan_y: 0.0 },
            Waypoint { t: 1.0, rot_x: 0.0, rot_y: 1.0, rot_z: 0.0, zoom: 1.0, pan_x: 0.0, pan_y: 0.0 },
        ];
        let cam_before = camera_at(&waypoints, -0.5);
        let cam_after = camera_at(&waypoints, 1.5);
        assert!((cam_before.rot_y - 0.0).abs() < 1e-9);
        assert!((cam_after.rot_y - 1.0).abs() < 1e-9);
    }

    #[test]
    fn writes_png_sequence_to_output_dir() {
        let tmp = tempfile::tempdir().expect("tempdir");
        let out = tmp.path().to_str().unwrap().to_string();

        let cfg = BatchConfig {
            input: "examples/1UBQ.pdb".to_string(),
            output_dir: out.clone(),
            frames: 5,
            width: 160,
            height: 120,
            render_mode: "fullhd".to_string(),
            color: "structure".to_string(),
            viz: "cartoon".to_string(),
            waypoints: vec![
                Waypoint { t: 0.0, rot_x: 0.0, rot_y: 0.0, rot_z: 0.0, zoom: 1.0, pan_x: 0.0, pan_y: 0.0 },
                Waypoint { t: 1.0, rot_x: 0.0, rot_y: 1.0, rot_z: 0.0, zoom: 1.0, pan_x: 0.0, pan_y: 0.0 },
            ],
            hotspots: vec![],
        };

        run(&cfg).expect("batch run");

        // Expect 5 files: frame_0001.png .. frame_0005.png
        for i in 1..=5 {
            let path = format!("{}/frame_{:04}.png", out, i);
            assert!(std::path::Path::new(&path).exists(), "missing {path}");
        }
    }

    #[test]
    fn hotspot_residues_are_visibly_recolored() {
        let tmp = tempfile::tempdir().expect("tempdir");
        let out = tmp.path().to_str().unwrap().to_string();

        // Two configs: identical except one has a bright-red hotspot at residue 1
        let base_cfg = BatchConfig {
            input: "examples/1UBQ.pdb".to_string(),
            output_dir: out.clone(),
            frames: 1,
            width: 320,
            height: 240,
            render_mode: "fullhd".to_string(),
            color: "structure".to_string(),
            viz: "cartoon".to_string(),
            waypoints: vec![Waypoint {
                t: 0.0, rot_x: 0.0, rot_y: 0.0, rot_z: 0.0, zoom: 1.0, pan_x: 0.0, pan_y: 0.0,
            }],
            hotspots: vec![],
        };

        let protein = crate::parser::pdb::load_structure(&base_cfg.input).expect("load");
        let cam = camera_at(&base_cfg.waypoints, 0.0);
        let no_highlight = render_frame(&protein, &cam, &base_cfg).expect("render");

        let mut hot_cfg = base_cfg.clone();
        hot_cfg.hotspots = vec![HotspotSpec {
            chain: "A".to_string(),
            residues: (1..=20).collect(), // first 20 residues of ubiquitin
            color: [255, 0, 0],
        }];
        let with_highlight = render_frame(&protein, &cam, &hot_cfg).expect("render");

        // Count strongly-red pixels in each (R > 200, G < 80, B < 80)
        let count_red = |img: &image::RgbImage| -> usize {
            img.pixels().filter(|p| p.0[0] > 200 && p.0[1] < 80 && p.0[2] < 80).count()
        };
        let red_off = count_red(&no_highlight);
        let red_on = count_red(&with_highlight);

        assert!(
            red_on > red_off + 50,
            "expected hotspot to add red pixels: off={red_off} on={red_on}"
        );
    }

    #[test]
    fn renders_one_frame_to_image() {
        let cfg = BatchConfig {
            input: "examples/1UBQ.pdb".to_string(),
            output_dir: "/tmp/proteinview-batch-test".to_string(),
            frames: 1,
            width: 320,
            height: 240,
            render_mode: "fullhd".to_string(),
            color: "structure".to_string(),
            viz: "cartoon".to_string(),
            waypoints: vec![Waypoint {
                t: 0.0, rot_x: 0.0, rot_y: 0.0, rot_z: 0.0, zoom: 1.0, pan_x: 0.0, pan_y: 0.0,
            }],
            hotspots: vec![],
        };

        let protein = crate::parser::pdb::load_structure(&cfg.input).expect("load");
        let cam = camera_at(&cfg.waypoints, 0.0);
        let image = render_frame(&protein, &cam, &cfg).expect("render");

        assert_eq!(image.width(), 320);
        assert_eq!(image.height(), 240);

        // Heuristic non-empty check: at least 1% of pixels are non-black
        let total = (image.width() * image.height()) as usize;
        let non_black = image.pixels().filter(|p| p.0[0] != 0 || p.0[1] != 0 || p.0[2] != 0).count();
        assert!(
            non_black * 100 >= total,
            "rendered image too sparse: {non_black}/{total} non-black pixels"
        );
    }
}

//! Batch (headless) rendering mode for ProteinView.
//!
//! Reads a JSON config describing input PDB, camera waypoints, color scheme,
//! hotspot residues, and output dir; renders a PNG sequence without entering
//! the TUI.

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
}

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

fn default_render_mode() -> String {
    "fullhd".to_string()
}
fn default_color() -> String {
    "structure".to_string()
}
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
        assert_eq!(cfg.waypoints.len(), 2);
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
}

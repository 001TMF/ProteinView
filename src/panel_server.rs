//! Persistent headless renderer for Codex- or agent-owned terminal panels.
//!
//! The server owns parsed molecular data and camera/presentation state, but it
//! never owns the terminal. Requests arrive as size-capped NDJSON on stdin.
//! Each successful state mutation renders through ProteinView's FullHD software
//! framebuffer, atomically replaces one caller-owned PNG path, and only then
//! acknowledges the new revision on stdout.

use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use image::DynamicImage;
use serde::Deserialize;
use serde_json::{Value, json};

use crate::app::{LARGE_STRUCTURE_THRESHOLD, VizMode};
use crate::model::interface::{InterfaceAnalysis, analyze_interface_for_chain};
use crate::model::protein::Protein;
use crate::model::selection::{
    MAX_RESIDUE_COLORS, ResidueColorOverrides, ResidueColorSpec, format_rgb, parse_rgb,
    resolve_residue_colors,
};
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, ColorSchemeType};
use crate::render::hd;
use crate::render::ribbon::{RibbonTriangle, generate_ribbon_mesh};
use crate::render::snapshot::{
    MAX_SNAPSHOT_DIMENSION, MAX_SNAPSHOT_PIXELS, MIN_SNAPSHOT_DIMENSION, validate_dimensions,
    write_png_atomically,
};

const PROTOCOL_NAME: &str = "proteinview-panel";
const PROTOCOL_VERSION: u8 = 1;
pub const MAX_REQUEST_BYTES: usize = 64 * 1024;
pub const MAX_RESPONSE_BYTES: usize = 64 * 1024;
const MAX_STATE_NAME_BYTES: usize = 256;
const MAX_STATE_RESIDUE_NAME_BYTES: usize = 32;
const MAX_REQUEST_ID_BYTES: usize = 128;
const MAX_ERROR_MESSAGE_BYTES: usize = 512;

#[derive(Debug, Clone)]
pub struct PanelServerOptions {
    pub width: u32,
    pub height: u32,
    pub color_override: Option<ColorSchemeType>,
    pub residue_colors: ResidueColorOverrides,
    pub viz_mode: VizMode,
    pub user_explicit_mode: bool,
}

#[derive(Debug, Clone)]
struct PanelSettings {
    camera: Camera,
    width: u32,
    height: u32,
    base_color: ColorSchemeType,
    residue_colors: ResidueColorOverrides,
    viz_mode: VizMode,
    current_chain: usize,
    show_interface: bool,
    show_interactions: bool,
    show_ligands: bool,
}

struct PanelSession {
    protein: Protein,
    settings: PanelSettings,
    interface_cache: Vec<Option<InterfaceAnalysis>>,
    has_plddt: bool,
    color_scheme: ColorScheme,
    color_dirty: bool,
    mesh_cache: Vec<RibbonTriangle>,
    mesh_dirty: bool,
    output_path: PathBuf,
    revision: u64,
}

impl PanelSession {
    fn new(mut protein: Protein, output_path: &Path, options: PanelServerOptions) -> Result<Self> {
        validate_dimensions(options.width, options.height)?;
        let output_path = absolute_path(output_path)?;
        let output_parent = output_path.parent().unwrap_or_else(|| Path::new("."));
        if !output_parent.is_dir() {
            anyhow::bail!(
                "panel output directory does not exist: {}",
                output_parent.display()
            );
        }

        protein.center();
        let total_residues = protein.residue_count();
        let has_plddt =
            protein.has_plddt() || options.color_override == Some(ColorSchemeType::Plddt);
        let viz_mode = if total_residues > LARGE_STRUCTURE_THRESHOLD
            && !options.user_explicit_mode
            && options.viz_mode == VizMode::Cartoon
        {
            VizMode::Backbone
        } else {
            options.viz_mode
        };
        let base_color = options.color_override.unwrap_or(ColorSchemeType::Structure);

        // Focused interface analyses are expensive for large complexes, so
        // compute them lazily once per selected chain and retain the results.
        let interface_cache = std::iter::repeat_with(|| None)
            .take(protein.chains.len().max(1))
            .collect();

        let mut camera = Camera::default();
        camera.zoom = fit_zoom(&protein, options.width, options.height);
        let color_scheme = ColorScheme::new(base_color, total_residues)
            .with_residue_colors(options.residue_colors.clone());
        let (mesh_cache, mesh_dirty) = if viz_mode == VizMode::Cartoon {
            (generate_ribbon_mesh(&protein, &color_scheme), false)
        } else {
            (Vec::new(), true)
        };

        Ok(Self {
            protein,
            settings: PanelSettings {
                camera,
                width: options.width,
                height: options.height,
                base_color,
                residue_colors: options.residue_colors,
                viz_mode,
                current_chain: 0,
                show_interface: false,
                show_interactions: false,
                show_ligands: true,
            },
            interface_cache,
            has_plddt,
            color_scheme,
            color_dirty: false,
            mesh_cache,
            mesh_dirty,
            output_path,
            revision: 0,
        })
    }

    fn render_next_revision(&mut self) -> Result<()> {
        self.sync_render_state();
        let interactions = if self.settings.show_interface && self.settings.show_interactions {
            self.current_interface_analysis().interactions.as_slice()
        } else {
            &[]
        };
        let framebuffer = hd::render_hd_framebuffer(
            &self.protein,
            &self.settings.camera,
            &self.color_scheme,
            self.settings.viz_mode,
            f64::from(self.settings.width),
            f64::from(self.settings.height),
            &self.mesh_cache,
            self.settings.show_ligands,
            interactions,
            None,
        );
        let image = DynamicImage::ImageRgba8(framebuffer.to_rgba_image());
        write_png_atomically(&image, &self.output_path)?;
        self.revision = self
            .revision
            .checked_add(1)
            .context("panel frame revision exhausted")?;
        Ok(())
    }

    fn sync_render_state(&mut self) {
        if self.settings.show_interface {
            self.ensure_current_interface_analysis();
        }
        if self.color_dirty {
            self.color_scheme = if self.settings.show_interface {
                let analysis = self
                    .interface_cache
                    .get(self.settings.current_chain)
                    .and_then(Option::as_ref)
                    .expect("focused interface analysis was initialized");
                ColorScheme::new_interface(
                    self.protein.residue_count(),
                    self.settings.current_chain,
                    analysis,
                    &self.protein,
                )
            } else {
                ColorScheme::new(self.settings.base_color, self.protein.residue_count())
            }
            .with_residue_colors(self.settings.residue_colors.clone());
            self.color_dirty = false;
            self.mesh_dirty = true;
        }
        if self.settings.viz_mode == VizMode::Cartoon && self.mesh_dirty {
            self.mesh_cache = generate_ribbon_mesh(&self.protein, &self.color_scheme);
            self.mesh_dirty = false;
        }
    }

    fn ensure_current_interface_analysis(&mut self) {
        let index = self
            .settings
            .current_chain
            .min(self.interface_cache.len().saturating_sub(1));
        if self.interface_cache[index].is_none() {
            self.interface_cache[index] =
                Some(analyze_interface_for_chain(&self.protein, 4.5, index));
        }
    }

    fn current_interface_analysis(&self) -> &InterfaceAnalysis {
        self.interface_cache
            .get(self.settings.current_chain)
            .and_then(Option::as_ref)
            .expect("focused interface analysis was initialized")
    }

    fn restore_after_render_failure(&mut self, settings: PanelSettings) {
        self.settings = settings;
        self.color_dirty = true;
        self.mesh_dirty = true;
    }

    fn fit(&mut self) {
        self.settings.camera.zoom =
            fit_zoom(&self.protein, self.settings.width, self.settings.height);
    }

    fn state_json(&self) -> Value {
        let current_chain_id = self
            .protein
            .chains
            .get(self.settings.current_chain)
            .map(|chain| chain.id.as_str());
        let residue_colors = self
            .settings
            .residue_colors
            .entries()
            .iter()
            .map(|entry| {
                json!({
                    "chain": entry.chain_id,
                    "residue_number": entry.residue_number,
                    "insertion_code": entry.insertion_code,
                    "residue_name": bounded_text(
                        &entry.residue_name,
                        MAX_STATE_RESIDUE_NAME_BYTES,
                        "residue",
                    ),
                    "color": format_rgb(entry.color),
                })
            })
            .collect::<Vec<_>>();
        json!({
            "structure": {
                "name": bounded_display_name(&self.protein.name),
                "chain_count": self.protein.chains.len(),
                "residue_count": self.protein.residue_count(),
                "atom_count": self.protein.atom_count(),
                "ligand_count": self.protein.ligand_count(),
                "chains": self.protein.chains.iter().map(|chain| chain.id.as_str()).collect::<Vec<_>>(),
            },
            "camera": {
                "rot_x": self.settings.camera.rot_x,
                "rot_y": self.settings.camera.rot_y,
                "rot_z": self.settings.camera.rot_z,
                "zoom": self.settings.camera.zoom,
                "pan_x": self.settings.camera.pan_x,
                "pan_y": self.settings.camera.pan_y,
                "auto_rotate": self.settings.camera.auto_rotate,
            },
            "viewport": {
                "width": self.settings.width,
                "height": self.settings.height,
            },
            "presentation": {
                "viz_mode": viz_mode_name(self.settings.viz_mode),
                "color": color_name(self.settings.base_color),
                "effective_color": if self.settings.show_interface {
                    "interface"
                } else {
                    color_name(self.settings.base_color)
                },
                "current_chain_index": self.settings.current_chain,
                "current_chain_id": current_chain_id,
                "interface": self.settings.show_interface,
                "interactions": self.settings.show_interactions,
                "ligands": self.settings.show_ligands,
                "residue_colors": residue_colors,
            },
        })
    }

    fn frame_json(&self) -> Value {
        json!({
            "path": self.output_path.to_string_lossy(),
            "mime_type": "image/png",
            "width": self.settings.width,
            "height": self.settings.height,
        })
    }
}

fn fit_zoom(protein: &Protein, width: u32, height: u32) -> f64 {
    let radius = protein.bounding_radius().max(1.0);
    0.9 * f64::from(width.min(height)) / (2.0 * radius)
}

fn absolute_path(path: &Path) -> Result<PathBuf> {
    if path.is_absolute() {
        Ok(path.to_path_buf())
    } else {
        Ok(std::env::current_dir()
            .context("failed to resolve panel output path")?
            .join(path))
    }
}

fn bounded_display_name(value: &str) -> String {
    bounded_text(value, MAX_STATE_NAME_BYTES, "structure")
}

fn bounded_text(value: &str, max_bytes: usize, fallback: &str) -> String {
    let mut output = String::new();
    for character in value.chars().take(max_bytes) {
        let character = if character.is_control() {
            ' '
        } else {
            character
        };
        if output.len().saturating_add(character.len_utf8()) > max_bytes {
            break;
        }
        output.push(character);
    }
    let trimmed = output.trim();
    if trimmed.is_empty() {
        fallback.to_string()
    } else {
        trimmed.to_string()
    }
}

#[derive(Debug, Deserialize)]
struct PanelRequest {
    #[serde(default = "protocol_version")]
    version: u8,
    id: Value,
    #[serde(flatten)]
    command: PanelCommand,
}

fn protocol_version() -> u8 {
    PROTOCOL_VERSION
}

#[derive(Debug, Deserialize)]
#[serde(tag = "command", rename_all = "snake_case")]
enum PanelCommand {
    Render,
    Resize {
        width: u32,
        height: u32,
    },
    Rotate {
        axis: Axis,
        delta: f64,
    },
    Pan {
        dx: f64,
        dy: f64,
    },
    Zoom {
        factor: Option<f64>,
        direction: Option<ZoomDirection>,
    },
    Reset,
    Fit,
    SetColor {
        color: WireColor,
    },
    SetResidueColors {
        residues: Vec<WireResidueColor>,
    },
    CycleColor,
    SetViz {
        mode: WireVizMode,
    },
    CycleViz,
    SelectChain {
        direction: ChainDirection,
    },
    SetChain {
        chain: String,
    },
    ToggleInterface,
    SetInterface {
        enabled: bool,
    },
    ToggleInteractions,
    SetInteractions {
        enabled: bool,
    },
    ToggleLigands,
    SetLigands {
        enabled: bool,
    },
    ToggleAutoRotate,
    SetAutoRotate {
        enabled: bool,
    },
    GetState,
    Shutdown,
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "lowercase")]
enum Axis {
    X,
    Y,
    Z,
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "lowercase")]
enum ZoomDirection {
    In,
    Out,
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "lowercase")]
enum ChainDirection {
    Prev,
    Next,
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "lowercase")]
enum WireColor {
    Structure,
    Element,
    Chain,
    Bfactor,
    Rainbow,
    Plddt,
}

#[derive(Debug, Clone, Deserialize)]
struct WireResidueColor {
    chain: String,
    residue_number: i32,
    #[serde(default)]
    insertion_code: Option<String>,
    color: String,
}

impl WireResidueColor {
    fn into_spec(self) -> Result<ResidueColorSpec> {
        ResidueColorSpec::new(
            self.chain,
            self.residue_number,
            self.insertion_code,
            parse_rgb(&self.color)?,
        )
    }
}

impl From<WireColor> for ColorSchemeType {
    fn from(value: WireColor) -> Self {
        match value {
            WireColor::Structure => Self::Structure,
            WireColor::Element => Self::Element,
            WireColor::Chain => Self::Chain,
            WireColor::Bfactor => Self::BFactor,
            WireColor::Rainbow => Self::Rainbow,
            WireColor::Plddt => Self::Plddt,
        }
    }
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "lowercase")]
enum WireVizMode {
    Backbone,
    Cartoon,
    Wireframe,
}

impl From<WireVizMode> for VizMode {
    fn from(value: WireVizMode) -> Self {
        match value {
            WireVizMode::Backbone => Self::Backbone,
            WireVizMode::Cartoon => Self::Cartoon,
            WireVizMode::Wireframe => Self::Wireframe,
        }
    }
}

enum CommandOutcome {
    Render,
    Respond,
    Shutdown,
}

#[derive(Debug)]
struct ProtocolError {
    code: &'static str,
    message: String,
}

impl ProtocolError {
    fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }
}

fn apply_command(
    session: &mut PanelSession,
    command: PanelCommand,
) -> std::result::Result<CommandOutcome, ProtocolError> {
    match command {
        PanelCommand::Render => {
            if session.settings.camera.auto_rotate {
                session.settings.camera.tick();
            }
            Ok(CommandOutcome::Render)
        }
        PanelCommand::Resize { width, height } => {
            validate_dimensions(width, height)
                .map_err(|error| ProtocolError::new("invalid_params", error.to_string()))?;
            session.settings.width = width;
            session.settings.height = height;
            session.fit();
            Ok(CommandOutcome::Render)
        }
        PanelCommand::Rotate { axis, delta } => {
            validate_finite("delta", delta)?;
            match axis {
                Axis::X => {
                    let result = session.settings.camera.rot_x + delta;
                    validate_finite("resulting x rotation", result)?;
                    session.settings.camera.rot_x = result;
                }
                Axis::Y => {
                    let result = session.settings.camera.rot_y + delta;
                    validate_finite("resulting y rotation", result)?;
                    session.settings.camera.rot_y = result;
                }
                Axis::Z => {
                    let result = session.settings.camera.rot_z + delta;
                    validate_finite("resulting z rotation", result)?;
                    session.settings.camera.rot_z = result;
                }
            }
            session.settings.camera.reset_tick_timer();
            Ok(CommandOutcome::Render)
        }
        PanelCommand::Pan { dx, dy } => {
            validate_finite("dx", dx)?;
            validate_finite("dy", dy)?;
            let pan_x = session.settings.camera.pan_x + dx;
            let pan_y = session.settings.camera.pan_y + dy;
            validate_finite("resulting x pan", pan_x)?;
            validate_finite("resulting y pan", pan_y)?;
            session.settings.camera.pan_x = pan_x;
            session.settings.camera.pan_y = pan_y;
            session.settings.camera.reset_tick_timer();
            Ok(CommandOutcome::Render)
        }
        PanelCommand::Zoom { factor, direction } => {
            match (factor, direction) {
                (Some(factor), None) => {
                    validate_finite("factor", factor)?;
                    if factor <= 0.0 {
                        return Err(ProtocolError::new(
                            "invalid_params",
                            "zoom factor must be greater than zero",
                        ));
                    }
                    let zoom = session.settings.camera.zoom * factor;
                    validate_positive_finite("resulting zoom", zoom)?;
                    session.settings.camera.zoom = zoom;
                }
                (None, Some(direction)) => {
                    let factor = match direction {
                        ZoomDirection::In => 1.1,
                        ZoomDirection::Out => 0.9,
                    };
                    let zoom = session.settings.camera.zoom * factor;
                    validate_positive_finite("resulting zoom", zoom)?;
                    session.settings.camera.zoom = zoom;
                }
                _ => {
                    return Err(ProtocolError::new(
                        "invalid_params",
                        "zoom requires exactly one of factor or direction",
                    ));
                }
            }
            session.settings.camera.reset_tick_timer();
            Ok(CommandOutcome::Render)
        }
        PanelCommand::Reset => {
            session.settings.camera.reset();
            session.fit();
            Ok(CommandOutcome::Render)
        }
        PanelCommand::Fit => {
            session.fit();
            Ok(CommandOutcome::Render)
        }
        PanelCommand::SetColor { color } => {
            let color = color.into();
            if session.settings.base_color != color {
                session.settings.base_color = color;
                if !session.settings.show_interface {
                    session.color_dirty = true;
                }
            }
            Ok(CommandOutcome::Render)
        }
        PanelCommand::SetResidueColors { residues } => {
            if residues.len() > MAX_RESIDUE_COLORS {
                return Err(ProtocolError::new(
                    "invalid_params",
                    format!("at most {MAX_RESIDUE_COLORS} exact residue colors are allowed"),
                ));
            }
            let specs = residues
                .into_iter()
                .map(WireResidueColor::into_spec)
                .collect::<Result<Vec<_>>>()
                .map_err(|error| ProtocolError::new("invalid_params", error.to_string()))?;
            session.settings.residue_colors = resolve_residue_colors(&session.protein, &specs)
                .map_err(|error| ProtocolError::new("invalid_params", error.to_string()))?;
            session.color_dirty = true;
            Ok(CommandOutcome::Render)
        }
        PanelCommand::CycleColor => {
            session.settings.base_color = session.settings.base_color.next(session.has_plddt);
            if !session.settings.show_interface {
                session.color_dirty = true;
            }
            Ok(CommandOutcome::Render)
        }
        PanelCommand::SetViz { mode } => {
            session.settings.viz_mode = mode.into();
            Ok(CommandOutcome::Render)
        }
        PanelCommand::CycleViz => {
            session.settings.viz_mode = session.settings.viz_mode.next();
            Ok(CommandOutcome::Render)
        }
        PanelCommand::SelectChain { direction } => {
            let chain_count = session.protein.chains.len();
            if chain_count > 0 {
                session.settings.current_chain = match direction {
                    ChainDirection::Prev if session.settings.current_chain == 0 => chain_count - 1,
                    ChainDirection::Prev => session.settings.current_chain - 1,
                    ChainDirection::Next => (session.settings.current_chain + 1) % chain_count,
                };
                if session.settings.show_interface {
                    session.color_dirty = true;
                }
            }
            Ok(CommandOutcome::Render)
        }
        PanelCommand::SetChain { chain } => {
            let chain_index = session
                .protein
                .chains
                .iter()
                .position(|candidate| candidate.id == chain)
                .ok_or_else(|| {
                    let available = session
                        .protein
                        .chains
                        .iter()
                        .map(|candidate| candidate.id.as_str())
                        .collect::<Vec<_>>()
                        .join(", ");
                    ProtocolError::new(
                        "invalid_params",
                        format!("chain {chain:?} was not found; available chains: {available}"),
                    )
                })?;
            session.settings.current_chain = chain_index;
            if session.settings.show_interface {
                session.color_dirty = true;
            }
            Ok(CommandOutcome::Render)
        }
        PanelCommand::ToggleInterface => {
            set_interface(session, !session.settings.show_interface);
            Ok(CommandOutcome::Render)
        }
        PanelCommand::SetInterface { enabled } => {
            set_interface(session, enabled);
            Ok(CommandOutcome::Render)
        }
        PanelCommand::ToggleInteractions => {
            if !session.settings.show_interface {
                return Err(ProtocolError::new(
                    "invalid_state",
                    "interactions require interface mode",
                ));
            }
            session.settings.show_interactions = !session.settings.show_interactions;
            Ok(CommandOutcome::Render)
        }
        PanelCommand::SetInteractions { enabled } => {
            if enabled && !session.settings.show_interface {
                return Err(ProtocolError::new(
                    "invalid_state",
                    "interactions require interface mode",
                ));
            }
            session.settings.show_interactions = enabled;
            Ok(CommandOutcome::Render)
        }
        PanelCommand::ToggleLigands => {
            session.settings.show_ligands = !session.settings.show_ligands;
            Ok(CommandOutcome::Render)
        }
        PanelCommand::SetLigands { enabled } => {
            session.settings.show_ligands = enabled;
            Ok(CommandOutcome::Render)
        }
        PanelCommand::ToggleAutoRotate => {
            session.settings.camera.auto_rotate = !session.settings.camera.auto_rotate;
            session.settings.camera.reset_tick_timer();
            Ok(CommandOutcome::Render)
        }
        PanelCommand::SetAutoRotate { enabled } => {
            session.settings.camera.auto_rotate = enabled;
            session.settings.camera.reset_tick_timer();
            Ok(CommandOutcome::Render)
        }
        PanelCommand::GetState => Ok(CommandOutcome::Respond),
        PanelCommand::Shutdown => Ok(CommandOutcome::Shutdown),
    }
}

fn set_interface(session: &mut PanelSession, enabled: bool) {
    if session.settings.show_interface != enabled {
        session.settings.show_interface = enabled;
        session.color_dirty = true;
    }
    if !enabled {
        session.settings.show_interactions = false;
    }
}

fn validate_finite(name: &str, value: f64) -> std::result::Result<(), ProtocolError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(ProtocolError::new(
            "invalid_params",
            format!("{name} must be finite"),
        ))
    }
}

fn validate_positive_finite(name: &str, value: f64) -> std::result::Result<(), ProtocolError> {
    validate_finite(name, value)?;
    if value > 0.0 {
        Ok(())
    } else {
        Err(ProtocolError::new(
            "invalid_params",
            format!("{name} must be greater than zero"),
        ))
    }
}

fn viz_mode_name(mode: VizMode) -> &'static str {
    match mode {
        VizMode::Backbone => "backbone",
        VizMode::Cartoon => "cartoon",
        VizMode::Wireframe => "wireframe",
    }
}

fn color_name(color: ColorSchemeType) -> &'static str {
    match color {
        ColorSchemeType::Structure => "structure",
        ColorSchemeType::Chain => "chain",
        ColorSchemeType::Element => "element",
        ColorSchemeType::BFactor => "bfactor",
        ColorSchemeType::Rainbow => "rainbow",
        ColorSchemeType::Interface => "interface",
        ColorSchemeType::Plddt => "plddt",
    }
}

fn is_supported_command(command: &str) -> bool {
    matches!(
        command,
        "render"
            | "resize"
            | "rotate"
            | "pan"
            | "zoom"
            | "reset"
            | "fit"
            | "set_color"
            | "set_residue_colors"
            | "cycle_color"
            | "set_viz"
            | "cycle_viz"
            | "select_chain"
            | "set_chain"
            | "toggle_interface"
            | "set_interface"
            | "toggle_interactions"
            | "set_interactions"
            | "toggle_ligands"
            | "set_ligands"
            | "toggle_auto_rotate"
            | "set_auto_rotate"
            | "get_state"
            | "shutdown"
    )
}

enum CappedLine {
    Eof,
    Line(Vec<u8>),
    TooLarge,
}

fn read_capped_line(reader: &mut impl BufRead) -> io::Result<CappedLine> {
    let mut line = Vec::new();
    let mut saw_data = false;
    let mut too_large = false;

    loop {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            return if !saw_data {
                Ok(CappedLine::Eof)
            } else if too_large {
                Ok(CappedLine::TooLarge)
            } else {
                Ok(CappedLine::Line(line))
            };
        }

        let newline = available.iter().position(|byte| *byte == b'\n');
        let take = newline.unwrap_or(available.len());
        saw_data |= take > 0 || newline.is_some();
        if !too_large {
            if line.len().saturating_add(take) > MAX_REQUEST_BYTES {
                too_large = true;
                line.clear();
            } else {
                line.extend_from_slice(&available[..take]);
            }
        }
        let consumed = take + usize::from(newline.is_some());
        reader.consume(consumed);

        if newline.is_some() {
            if too_large {
                return Ok(CappedLine::TooLarge);
            }
            if line.last() == Some(&b'\r') {
                line.pop();
            }
            return Ok(CappedLine::Line(line));
        }
    }
}

fn write_record(writer: &mut impl Write, record: &Value) -> Result<()> {
    let encoded = encode_record(record)?;
    writer
        .write_all(&encoded)
        .context("failed to write panel response")?;
    writer.flush().context("failed to flush panel response")
}

fn encode_record(record: &Value) -> Result<Vec<u8>> {
    let mut encoded = serde_json::to_vec(record).context("failed to encode panel response")?;
    encoded.push(b'\n');
    if encoded.len() > MAX_RESPONSE_BYTES {
        anyhow::bail!("panel response exceeds the {MAX_RESPONSE_BYTES}-byte limit");
    }
    Ok(encoded)
}

fn ready_record(session: &PanelSession) -> Value {
    json!({
        "type": "ready",
        "protocol": PROTOCOL_NAME,
        "version": PROTOCOL_VERSION,
        "revision": session.revision,
        "state": session.state_json(),
        "frame": session.frame_json(),
        "limits": {
            "min_dimension": MIN_SNAPSHOT_DIMENSION,
            "max_dimension": MAX_SNAPSHOT_DIMENSION,
            "max_pixels": MAX_SNAPSHOT_PIXELS,
            "max_request_bytes": MAX_REQUEST_BYTES,
        },
    })
}

fn success_record(session: &PanelSession, id: Value, revision: u64) -> Value {
    json!({
        "type": "response",
        "protocol": PROTOCOL_NAME,
        "version": PROTOCOL_VERSION,
        "id": id,
        "ok": true,
        "revision": revision,
        "state": session.state_json(),
        "frame": session.frame_json(),
    })
}

fn error_record(revision: u64, id: Value, code: &'static str, message: impl Into<String>) -> Value {
    let message = bounded_text(
        &message.into(),
        MAX_ERROR_MESSAGE_BYTES,
        "panel request failed",
    );
    json!({
        "type": "response",
        "protocol": PROTOCOL_NAME,
        "version": PROTOCOL_VERSION,
        "id": id,
        "ok": false,
        "revision": revision,
        "error": {
            "code": code,
            "message": message,
        },
    })
}

fn valid_request_id(id: &Value) -> bool {
    match id {
        Value::String(value) => value.len() <= MAX_REQUEST_ID_BYTES,
        Value::Number(number) => number.is_i64() || number.is_u64(),
        Value::Null | Value::Bool(_) | Value::Array(_) | Value::Object(_) => false,
    }
}

/// Run a persistent headless panel session until shutdown or stdin EOF.
pub fn serve(
    protein: Protein,
    output_path: &Path,
    options: PanelServerOptions,
    mut input: impl BufRead,
    mut output: impl Write,
) -> Result<()> {
    let mut session = PanelSession::new(protein, output_path, options)?;
    let mut preflight = Vec::new();
    write_record(&mut preflight, &ready_record(&session))?;
    session.render_next_revision()?;
    write_record(&mut output, &ready_record(&session))?;

    loop {
        let line = match read_capped_line(&mut input)? {
            CappedLine::Eof => break,
            CappedLine::TooLarge => {
                write_record(
                    &mut output,
                    &error_record(
                        session.revision,
                        Value::Null,
                        "request_too_large",
                        format!("request exceeds the {MAX_REQUEST_BYTES}-byte limit"),
                    ),
                )?;
                continue;
            }
            CappedLine::Line(line) if line.iter().all(u8::is_ascii_whitespace) => continue,
            CappedLine::Line(line) => line,
        };

        let value: Value = match serde_json::from_slice(&line) {
            Ok(value) => value,
            Err(error) => {
                write_record(
                    &mut output,
                    &error_record(
                        session.revision,
                        Value::Null,
                        "invalid_json",
                        error.to_string(),
                    ),
                )?;
                continue;
            }
        };
        let id = value.get("id").cloned().unwrap_or(Value::Null);
        if !valid_request_id(&id) {
            write_record(
                &mut output,
                &error_record(
                    session.revision,
                    Value::Null,
                    "invalid_request",
                    format!(
                        "request id must be an integer or a string no longer than \
                         {MAX_REQUEST_ID_BYTES} bytes"
                    ),
                ),
            )?;
            continue;
        }
        let Some(command_name) = value.get("command").and_then(Value::as_str) else {
            write_record(
                &mut output,
                &error_record(
                    session.revision,
                    id,
                    "invalid_request",
                    "request command must be a string",
                ),
            )?;
            continue;
        };
        if !is_supported_command(command_name) {
            write_record(
                &mut output,
                &error_record(
                    session.revision,
                    id,
                    "unknown_command",
                    format!("unknown panel command {command_name:?}"),
                ),
            )?;
            continue;
        }

        let request: PanelRequest = match serde_json::from_value(value) {
            Ok(request) => request,
            Err(error) => {
                write_record(
                    &mut output,
                    &error_record(session.revision, id, "invalid_params", error.to_string()),
                )?;
                continue;
            }
        };
        if request.version != PROTOCOL_VERSION {
            write_record(
                &mut output,
                &error_record(
                    session.revision,
                    request.id,
                    "invalid_request",
                    format!(
                        "unsupported panel protocol version {}; expected {PROTOCOL_VERSION}",
                        request.version
                    ),
                ),
            )?;
            continue;
        }

        let previous_settings = session.settings.clone();
        let outcome = match apply_command(&mut session, request.command) {
            Ok(outcome) => outcome,
            Err(error) => {
                write_record(
                    &mut output,
                    &error_record(session.revision, request.id, error.code, error.message),
                )?;
                continue;
            }
        };

        let response_revision = if matches!(outcome, CommandOutcome::Render) {
            session
                .revision
                .checked_add(1)
                .context("panel frame revision exhausted")?
        } else {
            session.revision
        };
        if encode_record(&success_record(
            &session,
            request.id.clone(),
            response_revision,
        ))
        .is_err()
        {
            session.restore_after_render_failure(previous_settings);
            write_record(
                &mut output,
                &error_record(
                    session.revision,
                    request.id,
                    "response_too_large",
                    format!("response exceeds the {MAX_RESPONSE_BYTES}-byte limit"),
                ),
            )?;
            continue;
        }

        if matches!(outcome, CommandOutcome::Render) {
            if let Err(error) = session.render_next_revision() {
                session.restore_after_render_failure(previous_settings);
                write_record(
                    &mut output,
                    &error_record(
                        session.revision,
                        request.id,
                        "render_failed",
                        error.to_string(),
                    ),
                )?;
                continue;
            }
        }

        let should_shutdown = matches!(outcome, CommandOutcome::Shutdown);
        write_record(
            &mut output,
            &success_record(&session, request.id, session.revision),
        )?;
        if should_shutdown {
            break;
        }
    }

    Ok(())
}

#[cfg(test)]
#[path = "panel_server_tests.rs"]
mod tests;

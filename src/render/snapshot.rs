use std::fs::{self, OpenOptions};
use std::io::BufWriter;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{Context, Result, bail};
use image::{DynamicImage, ImageFormat};

use crate::app::{LARGE_STRUCTURE_THRESHOLD, VizMode};
use crate::model::interface::{Interaction, analyze_interface_for_chain};
use crate::model::protein::Protein;
use crate::model::selection::ResidueColorOverrides;
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, ColorSchemeType};
use crate::render::hd;
use crate::render::ribbon::generate_ribbon_mesh;

pub const DEFAULT_SNAPSHOT_WIDTH: u32 = 1920;
pub const DEFAULT_SNAPSHOT_HEIGHT: u32 = 1080;

pub const MIN_SNAPSHOT_DIMENSION: u32 = 64;
pub const MAX_SNAPSHOT_DIMENSION: u32 = 4096;
pub const MAX_SNAPSHOT_PIXELS: u64 = 8_388_608;
static TEMP_FILE_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Options for a single non-interactive FullHD pixel render.
#[derive(Debug, Clone)]
pub struct SnapshotOptions {
    pub width: u32,
    pub height: u32,
    pub color_override: Option<ColorSchemeType>,
    pub residue_colors: ResidueColorOverrides,
    pub viz_mode: VizMode,
    pub user_explicit_mode: bool,
    pub show_ligands: bool,
    pub interface_chain: Option<String>,
    pub show_interactions: bool,
}

/// Rasterize one ProteinView FullHD frame and write it as a PNG.
///
/// This path uses the same software framebuffer as the interactive FullHD
/// renderer, but deliberately skips terminal detection and alternate-screen
/// setup so agent tools can invoke it without nesting one TUI inside another.
pub fn save_png(mut protein: Protein, output_path: &Path, options: SnapshotOptions) -> Result<()> {
    validate_dimensions(options.width, options.height)?;
    if options.show_interactions && options.interface_chain.is_none() {
        bail!("snapshot interactions require an interface focus chain");
    }

    protein.center();
    let radius = protein.bounding_radius().max(1.0);
    let total_residues = protein.residue_count();
    let viz_mode = if total_residues > LARGE_STRUCTURE_THRESHOLD
        && !options.user_explicit_mode
        && options.viz_mode == VizMode::Cartoon
    {
        VizMode::Backbone
    } else {
        options.viz_mode
    };

    let focus_chain = options
        .interface_chain
        .as_deref()
        .map(|requested| {
            protein
                .chains
                .iter()
                .position(|chain| chain.id == requested)
                .with_context(|| {
                    let available = protein
                        .chains
                        .iter()
                        .map(|chain| chain.id.as_str())
                        .collect::<Vec<_>>()
                        .join(", ");
                    format!(
                        "snapshot interface chain {requested:?} was not found; available chains: {available}"
                    )
                })
        })
        .transpose()?;
    let interface_analysis =
        focus_chain.map(|chain| analyze_interface_for_chain(&protein, 4.5, chain));
    let color_scheme =
        if let (Some(focus_chain), Some(analysis)) = (focus_chain, interface_analysis.as_ref()) {
            ColorScheme::new_interface(total_residues, focus_chain, analysis, &protein)
        } else {
            ColorScheme::new(
                options.color_override.unwrap_or(ColorSchemeType::Structure),
                total_residues,
            )
        }
        .with_residue_colors(options.residue_colors);
    let mesh = if viz_mode == VizMode::Cartoon {
        generate_ribbon_mesh(&protein, &color_scheme)
    } else {
        Vec::new()
    };

    let mut camera = Camera::default();
    camera.zoom = 0.9 * f64::from(options.width.min(options.height)) / (2.0 * radius);

    let interactions: &[Interaction] = if options.show_interactions {
        interface_analysis
            .as_ref()
            .map(|analysis| analysis.interactions.as_slice())
            .unwrap_or_default()
    } else {
        &[]
    };
    let framebuffer = hd::render_hd_framebuffer(
        &protein,
        &camera,
        &color_scheme,
        viz_mode,
        f64::from(options.width),
        f64::from(options.height),
        &mesh,
        options.show_ligands,
        interactions,
        None,
    );
    let image = DynamicImage::ImageRgba8(framebuffer.to_rgba_image());
    write_png_atomically(&image, output_path)
}

pub(crate) fn validate_dimensions(width: u32, height: u32) -> Result<()> {
    if !(MIN_SNAPSHOT_DIMENSION..=MAX_SNAPSHOT_DIMENSION).contains(&width)
        || !(MIN_SNAPSHOT_DIMENSION..=MAX_SNAPSHOT_DIMENSION).contains(&height)
    {
        bail!(
            "snapshot dimensions must each be between {MIN_SNAPSHOT_DIMENSION} and {MAX_SNAPSHOT_DIMENSION} pixels"
        );
    }
    if u64::from(width) * u64::from(height) > MAX_SNAPSHOT_PIXELS {
        bail!("snapshot exceeds the {MAX_SNAPSHOT_PIXELS}-pixel safety limit");
    }
    Ok(())
}

pub(crate) fn write_png_atomically(image: &DynamicImage, output_path: &Path) -> Result<()> {
    let parent = output_path.parent().unwrap_or_else(|| Path::new("."));
    let output_name = output_path
        .file_name()
        .with_context(|| format!("snapshot path has no file name: {}", output_path.display()))?
        .to_string_lossy();

    let mut reserved = None;
    for _ in 0..100 {
        let sequence = TEMP_FILE_COUNTER.fetch_add(1, Ordering::Relaxed);
        let temp_path = parent.join(format!(
            ".{output_name}.{}.{}.tmp",
            std::process::id(),
            sequence
        ));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temp_path)
        {
            Ok(file) => {
                reserved = Some((temp_path, file));
                break;
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => {
                return Err(error).with_context(|| {
                    format!(
                        "failed to create temporary snapshot beside {}",
                        output_path.display()
                    )
                });
            }
        }
    }

    let (temp_path, file) = reserved.context("failed to reserve a temporary snapshot file")?;
    let write_result = (|| -> Result<()> {
        let mut writer = BufWriter::new(file);
        image
            .write_to(&mut writer, ImageFormat::Png)
            .context("failed to encode FullHD PNG")?;
        let file = writer.into_inner().context("failed to flush FullHD PNG")?;
        file.sync_all().context("failed to sync FullHD PNG")?;
        fs::rename(&temp_path, output_path).with_context(|| {
            format!(
                "failed to atomically move FullHD PNG into {}",
                output_path.display()
            )
        })
    })();

    if write_result.is_err() {
        let _ = fs::remove_file(&temp_path);
    }
    write_result
}

#[cfg(test)]
#[path = "snapshot_tests.rs"]
mod tests;

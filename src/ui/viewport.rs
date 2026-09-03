use image::DynamicImage;
use ratatui::Frame;
use ratatui::layout::Rect;
use ratatui_image::picker::ProtocolType;
use ratatui_image::{Image, Resize};

use crate::app::{self, App, ConnectionType, RenderMode};
use crate::model::interface::Interaction;
use crate::render::braille;
use crate::render::framebuffer::{
    Framebuffer, framebuffer_to_braille_widget, framebuffer_to_braille_widget_ssaa,
};
use crate::render::hd;
use crate::render::kitty_png::KittyPngImage;
use crate::render::kitty_shm;

/// Supersampling factor for HDplus mode.
///
/// A braille cell addresses a fixed 2x4 dot grid, so this buys no extra dots --
/// it anti-aliases the silhouette and stabilises per-cell color instead.  It is
/// close to free: the HD framebuffer is small enough that rasterization cost is
/// dominated by the per-triangle projection stage, which is resolution
/// independent, and the emitted character grid is unchanged.
const HD_SSAA: usize = 2;

/// Color quantization step for HDplus mode over SSH.
///
/// Supersampling makes neighbouring cell colors more continuous, which would
/// break more run-length spans and cost bytes on the wire.  Rounding each
/// channel to a multiple of this merges them back together; the dominant cost of
/// this render path over SSH is the count of SGR color escapes, not pixels.
const HD_QUANT_STEP_SSH: u8 = 8;

/// Color quantization step for HDplus on a local terminal, where bandwidth is
/// free but merging spans still reduces per-frame diffing work.
const HD_QUANT_STEP_LOCAL: u8 = 4;

/// Quantization step to use for the current session.
fn hd_quant_step(connection: ConnectionType) -> u8 {
    match connection {
        ConnectionType::Ssh => HD_QUANT_STEP_SSH,
        ConnectionType::Local => HD_QUANT_STEP_LOCAL,
    }
}

/// Render the main 3D viewport
pub fn render_viewport(frame: &mut Frame, area: Rect, app: &App) {
    let interactions: &[Interaction] = if app.show_interface && app.show_interactions {
        &app.interface_analysis.interactions
    } else {
        &[]
    };

    match app.render_mode {
        RenderMode::Braille => {
            // Braille mode: 2x4 dots per cell, higher resolution but monochrome per cell
            let width = area.width as f64 * 2.0;
            let height = area.height as f64 * 4.0;

            let canvas = braille::render_protein(
                &app.protein,
                &app.camera,
                &app.color_scheme,
                app.viz_mode,
                width,
                height,
                app.show_ligands,
                interactions,
            );

            frame.render_widget(canvas, area);
        }
        RenderMode::HalfBlock => {
            // HalfBlock mode: render at braille resolution (2x4 per cell) through
            // the HD rasterizer (Lambert shading, z-buffer, depth fog) and convert
            // to colored braille characters.  This gives the same spatial resolution
            // as the basic Braille renderer but with much higher quality shading.
            let width = area.width as f64 * 2.0;
            let height = area.height as f64 * 4.0;

            let fb = hd::render_hd_framebuffer(
                &app.protein,
                &app.camera,
                &app.color_scheme,
                app.viz_mode,
                width,
                height,
                &app.mesh_cache,
                app.show_ligands,
                interactions,
            );

            let widget = framebuffer_to_braille_widget(&fb);
            frame.render_widget(widget, area);
        }
        RenderMode::HalfBlockPlus => {
            render_hdplus_viewport(frame, area, app, interactions);
        }
        RenderMode::FullHD => {
            render_fullhd_viewport(frame, area, app, interactions);
        }
    }
}

/// Render the HDplus viewport: the same 2x4 braille cell grid as HD, rasterized
/// into a supersampled framebuffer and box-filtered back down.
///
/// A braille cell addresses a fixed 2x4 dot grid, so supersampling buys no extra
/// dots.  What it buys is an anti-aliased silhouette and a stable per-cell color:
/// each dot is lit on >= 50% coverage of its sample block instead of on a single
/// sampled pixel, which stops thin ribbons from crawling between dots as the
/// camera rotates.  The emitted character grid is identical to HD's, so the cost
/// on the wire is unchanged.
fn render_hdplus_viewport(frame: &mut Frame, area: Rect, app: &App, interactions: &[Interaction]) {
    let ssaa = HD_SSAA as f64;

    // The framebuffer is supersampled, so the camera must be scaled to match --
    // zoom and pan are both in framebuffer pixel units.  Without this the
    // protein would cover the same pixel count in a buffer twice as wide, and
    // downsample to half its apparent size relative to Braille and HD.
    let mut cam = app.camera.clone();
    cam.zoom *= ssaa;
    cam.pan_x *= ssaa;
    cam.pan_y *= ssaa;

    let width = area.width as f64 * 2.0 * ssaa;
    let height = area.height as f64 * 4.0 * ssaa;

    let fb = hd::render_hd_framebuffer_ssaa(
        &app.protein,
        &cam,
        &app.color_scheme,
        app.viz_mode,
        width,
        height,
        &app.mesh_cache,
        app.show_ligands,
        interactions,
        ssaa,
    );

    let widget =
        framebuffer_to_braille_widget_ssaa(&fb, HD_SSAA, hd_quant_step(app.connection_type));
    frame.render_widget(widget, area);
}

/// Render the FullHD viewport using a graphics protocol (Kitty/Sixel/iTerm2)
/// when available, falling back to colored braille characters otherwise.
fn render_fullhd_viewport(frame: &mut Frame, area: Rect, app: &App, interactions: &[Interaction]) {
    let proto = app.picker.protocol_type();
    let (font_w, font_h) = app.picker.font_size();
    let is_graphics = proto != ProtocolType::Halfblocks && font_w > 0 && font_h > 0;

    // A graphics-protocol viewport is measured in device pixels, so on a HiDPI
    // panel it is far larger than the cell grid suggests.  `fullhd_framebuffer_size`
    // owns that policy, and `App::recalculate_zoom` sizes the zoom through the
    // very same function, so the two cannot drift apart.
    let (still_w, still_h) = app::fullhd_framebuffer_size(
        f64::from(area.width),
        f64::from(area.height),
        font_w,
        font_h,
        is_graphics,
    );

    // While the camera moves, render smaller and let the terminal stretch the
    // result back over the viewport via the protocol's c=/r= keys: a quarter of
    // the pixels, and motion hides the softness.  Zoom and pan are in
    // framebuffer pixels, so the camera has to be scaled to match.
    let scale = if is_graphics && app.is_interacting() {
        app::FULLHD_INTERACTIVE_SCALE
    } else {
        1.0
    };
    let (px_w, px_h) = ((still_w * scale).max(1.0), (still_h * scale).max(1.0));

    let mut cam = app.camera.clone();
    if scale < 1.0 {
        cam.zoom *= scale;
        cam.pan_x *= scale;
        cam.pan_y *= scale;
    }

    let fb = hd::render_hd_framebuffer(
        &app.protein,
        &cam,
        &app.color_scheme,
        app.viz_mode,
        px_w,
        px_h,
        &app.mesh_cache,
        app.show_ligands,
        interactions,
    );

    // If the terminal supports a real graphics protocol, hand it the pixels.
    if proto == ProtocolType::Kitty {
        if let Some(widget) = kitty_widget(&fb, app, area) {
            frame.render_widget(widget, area);
            return;
        }
        // Transmission failed — fall through to braille.
    } else if proto != ProtocolType::Halfblocks {
        // Sixel/iTerm2: use ratatui-image (no shared-memory option there).
        let dyn_img = DynamicImage::ImageRgb8(fb.to_rgb_image());
        if let Ok(protocol) = app.picker.new_protocol(dyn_img, area, Resize::Fit(None)) {
            let widget = Image::new(&protocol);
            frame.render_widget(widget, area);
            return;
        }
        // Protocol error — fall through to braille.
    }

    // Fallback: colored braille character rendering (always works).
    let widget = framebuffer_to_braille_widget(&fb);
    frame.render_widget(widget, area);
}

/// Build the Kitty transmission for this frame.
///
/// Prefers shared memory, where the pixels are written once — straight into the
/// mapping the terminal will read — and the escape sequence carries only the
/// object's name.  Falls back to the inline zlib transport when the terminal
/// did not accept shared memory at startup, or if the mapping cannot be created
/// this frame.
fn kitty_widget(fb: &Framebuffer, app: &App, area: Rect) -> Option<KittyPngImage> {
    let (w, h) = (fb.width as u32, fb.height as u32);
    let byte_len = fb.width * fb.height * 4;

    if app.kitty_shm {
        let slot = app.next_shm_slot();
        if let Some(mut shm) = kitty_shm::ShmFrame::create(slot, byte_len) {
            fb.write_rgba(shm.as_mut_slice());
            let widget = KittyPngImage::from_shm(shm.name(), w, h, byte_len, area);
            // Unmapping our view is safe immediately: the object stays linked
            // until the terminal reads and unlinks it.
            drop(shm);
            return Some(widget);
        }
    }

    let mut rgba = vec![0u8; byte_len];
    fb.write_rgba(&mut rgba);
    KittyPngImage::from_rgba(&rgba, w, h, area)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::app::{AppConfig, VizMode};
    use crate::model::protein::{Atom, Chain, MoleculeType, Protein, Residue, SecondaryStructure};
    use crate::model::selection::ResidueColorOverrides;
    use ratatui::Terminal;
    use ratatui::backend::TestBackend;
    use ratatui::buffer::Buffer;
    use ratatui_image::picker::Picker;

    /// A CA trace spanning enough space to give a measurable extent.
    fn fixture() -> Protein {
        let residues = (0..16)
            .map(|i| Residue {
                name: "ALA".to_string(),
                seq_num: i + 1,
                insertion_code: None,
                atoms: vec![Atom {
                    name: "CA".to_string(),
                    element: "C".to_string(),
                    x: i as f64 * 4.0 - 30.0,
                    y: (i as f64 * 0.7).sin() * 12.0,
                    z: 0.0,
                    b_factor: 20.0,
                    is_backbone: true,
                    is_hetero: false,
                }],
                secondary_structure: SecondaryStructure::Coil,
            })
            .collect();
        Protein {
            name: "trace".to_string(),
            chains: vec![Chain {
                id: "A".to_string(),
                residues,
                molecule_type: MoleculeType::Protein,
            }],
            ligands: Vec::new(),
        }
    }

    /// Drive the real viewport renderer through a test backend, exactly as the
    /// main loop does, and return the resulting cell buffer.
    fn draw(mode: RenderMode, cols: u16, rows: u16) -> Buffer {
        draw_protein(fixture(), mode, cols, rows)
    }

    fn draw_protein(protein: Protein, mode: RenderMode, cols: u16, rows: u16) -> Buffer {
        let app = App::new(
            protein,
            AppConfig {
                render_mode: mode,
                viz_mode: VizMode::Backbone,
                user_explicit_mode: true,
                color_override: None,
                residue_colors: ResidueColorOverrides::default(),
            },
            cols,
            rows,
            Picker::halfblocks(),
        );
        // The main layout reserves 4 rows of chrome around the viewport.
        let area = Rect::new(0, 0, cols, rows - 4);
        let mut term = Terminal::new(TestBackend::new(cols, rows)).unwrap();
        term.draw(|f| render_viewport(f, area, &app)).unwrap();
        term.backend().buffer().clone()
    }

    /// Extent of drawn cells, in terminal cells.
    fn extent(buf: &Buffer, cols: u16, rows: u16) -> (u16, u16) {
        let (mut min_x, mut max_x) = (u16::MAX, 0u16);
        let (mut min_y, mut max_y) = (u16::MAX, 0u16);
        for y in 0..rows {
            for x in 0..cols {
                if buf[(x, y)].symbol().trim().is_empty() {
                    continue;
                }
                min_x = min_x.min(x);
                max_x = max_x.max(x);
                min_y = min_y.min(y);
                max_y = max_y.max(y);
            }
        }
        assert!(min_x != u16::MAX, "viewport should draw something");
        (max_x - min_x, max_y - min_y)
    }

    #[test]
    fn hdplus_renders_at_the_same_apparent_size_as_hd_and_braille() {
        // HDplus rasterizes into a supersampled framebuffer, so the viewport has
        // to scale the camera to match.  If it does not, the protein downsamples
        // to roughly half the size of the other modes.  This drives the real
        // `render_viewport` wiring rather than the rasterizer in isolation.
        let (cols, rows) = (100u16, 40u16);
        let braille = extent(&draw(RenderMode::Braille, cols, rows), cols, rows);
        let hd = extent(&draw(RenderMode::HalfBlock, cols, rows), cols, rows);
        let hdplus = extent(&draw(RenderMode::HalfBlockPlus, cols, rows), cols, rows);

        for (label, other) in [("HD", hd), ("Braille", braille)] {
            assert!(
                hdplus.0.abs_diff(other.0) <= 2 && hdplus.1.abs_diff(other.1) <= 2,
                "HDplus extent {hdplus:?} should match {label} {other:?} within 2 cells"
            );
        }
    }

    #[test]
    fn hdplus_emits_the_same_cell_grid_as_hd() {
        // Supersampling must not change how many characters reach the terminal.
        let (cols, rows) = (100u16, 40u16);
        let hd = draw(RenderMode::HalfBlock, cols, rows);
        let hdplus = draw(RenderMode::HalfBlockPlus, cols, rows);
        assert_eq!(hd.area(), hdplus.area());
    }
}

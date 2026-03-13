use image::{RgbImage, RgbaImage};
use ratatui::style::{Color, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;

/// RGB pixel framebuffer with z-buffer for software rasterization.
///
/// Pixel coordinates: (0,0) is top-left, x increases right, y increases down.
/// The framebuffer dimensions are in *pixels*, not terminal cells.
/// For half-block rendering each terminal row maps to 2 pixel rows,
/// so `height` should typically be `terminal_rows * 2`.
pub struct Framebuffer {
    pub width: usize,
    pub height: usize,
    /// RGB color per pixel, row-major: index = y * width + x
    pub color: Vec<[u8; 3]>,
    /// Depth (z) per pixel for z-buffer tests. Smaller z = closer to viewer.
    pub depth: Vec<f64>,
}

/// A triangle in screen space, ready for rasterization.
pub struct Triangle {
    /// Three vertices in screen-space [x, y, z].
    /// x,y are pixel coordinates; z is depth for z-buffering.
    pub verts: [[f64; 3]; 3],
    /// Base RGB color before shading is applied.
    pub color: [u8; 3],
    /// Unit face normal in world/view space for Lambert shading.
    pub normal: [f64; 3],
}

impl Framebuffer {
    /// Create a new framebuffer initialized to black with infinite depth.
    ///
    /// Width and height are clamped to at least 1 to avoid creating an empty
    /// framebuffer where any `set_pixel` call would panic with index-out-of-bounds.
    pub fn new(width: usize, height: usize) -> Self {
        let width = width.max(1);
        let height = height.max(1);
        let n = width * height;
        Self {
            width,
            height,
            color: vec![[0, 0, 0]; n],
            depth: vec![f64::INFINITY; n],
        }
    }

    /// Reset the framebuffer to black pixels and infinite depth.
    #[cfg(test)]
    pub fn clear(&mut self) {
        for c in self.color.iter_mut() {
            *c = [0, 0, 0];
        }
        for d in self.depth.iter_mut() {
            *d = f64::INFINITY;
        }
    }

    /// Set a single pixel if it passes the z-buffer test.
    #[inline]
    fn set_pixel(&mut self, x: usize, y: usize, z: f64, color: [u8; 3]) {
        let idx = y * self.width + x;
        if z < self.depth[idx] {
            self.depth[idx] = z;
            self.color[idx] = color;
        }
    }

    /// Rasterize a single triangle with Lambert shading and z-buffering.
    ///
    /// `light_dir` should be a *unit* vector pointing toward the light source.
    /// The triangle's `normal` is expected to be a unit vector as well.
    ///
    /// Shading: `intensity = max(ambient, dot(normal, light_dir))` where
    /// ambient = 0.15. Each color channel is scaled by the intensity.
    ///
    /// Convenience wrapper around `rasterize_triangle_depth` for tests.
    #[cfg(test)]
    pub fn rasterize_triangle(&mut self, tri: &Triangle, light_dir: [f64; 3]) {
        self.rasterize_triangle_depth(tri, light_dir);
    }

    /// Rasterize a triangle with Lambert shading and z-buffering.
    ///
    /// `light_dir` should be a *unit* vector pointing toward the light source.
    /// The triangle's `normal` is expected to be a unit vector as well.
    ///
    /// Shading uses two-sided half-Lambert wrap lighting with an ambient term.
    /// Depth fog is handled separately via [`apply_depth_tint`] as a post-pass.
    pub fn rasterize_triangle_depth(
        &mut self,
        tri: &Triangle,
        light_dir: [f64; 3],
    ) {
        const AMBIENT: f64 = 0.55;

        // --- Two-sided Lambert shading with wrap lighting ---
        // Use abs(dot) so back-facing triangles also get proper lighting,
        // then apply a half-Lambert wrap to soften the falloff
        let dot = tri.normal[0] * light_dir[0]
            + tri.normal[1] * light_dir[1]
            + tri.normal[2] * light_dir[2];
        let half_lambert = dot.abs() * 0.4 + 0.6;
        let intensity = AMBIENT + (1.0 - AMBIENT) * half_lambert;

        // Precompute flat shaded color.
        let shaded: [u8; 3] = [
            (tri.color[0] as f64 * intensity).min(255.0) as u8,
            (tri.color[1] as f64 * intensity).min(255.0) as u8,
            (tri.color[2] as f64 * intensity).min(255.0) as u8,
        ];

        // --- Extract vertices ---
        let [v0, v1, v2] = tri.verts;

        // --- Bounding box (clamped to framebuffer) ---
        let min_x = v0[0].min(v1[0]).min(v2[0]).floor() as isize;
        let max_x = v0[0].max(v1[0]).max(v2[0]).ceil() as isize;
        let min_y = v0[1].min(v1[1]).min(v2[1]).floor() as isize;
        let max_y = v0[1].max(v1[1]).max(v2[1]).ceil() as isize;

        let min_x = min_x.max(0) as usize;
        let max_x = max_x.max(0).min(self.width as isize - 1) as usize;
        let min_y = min_y.max(0) as usize;
        let max_y = max_y.max(0).min(self.height as isize - 1) as usize;

        // Triangle entirely off-screen
        if min_x > max_x || min_y > max_y {
            return;
        }

        // --- Precompute barycentric denominator ---
        // For vertices A(v0), B(v1), C(v2), the signed area * 2:
        //   denom = (B.y - C.y)*(A.x - C.x) + (C.x - B.x)*(A.y - C.y)
        let denom = (v1[1] - v2[1]) * (v0[0] - v2[0]) + (v2[0] - v1[0]) * (v0[1] - v2[1]);
        if denom.abs() < 1e-12 {
            return; // degenerate triangle
        }
        let inv_denom = 1.0 / denom;

        // --- Rasterize pixels in bounding box ---
        for py in min_y..=max_y {
            let pf_y = py as f64 + 0.5; // pixel center
            for px in min_x..=max_x {
                let pf_x = px as f64 + 0.5; // pixel center

                // Barycentric coordinates
                let u =
                    ((v1[1] - v2[1]) * (pf_x - v2[0]) + (v2[0] - v1[0]) * (pf_y - v2[1]))
                        * inv_denom;
                let v =
                    ((v2[1] - v0[1]) * (pf_x - v2[0]) + (v0[0] - v2[0]) * (pf_y - v2[1]))
                        * inv_denom;
                let w = 1.0 - u - v;

                // Inside test (with a tiny epsilon for edge cases)
                if u >= -1e-6 && v >= -1e-6 && w >= -1e-6 {
                    // Interpolate z
                    let z = u * v0[2] + v * v1[2] + w * v2[2];
                    self.set_pixel(px, py, z, shaded);
                }
            }
        }
    }

    /// Apply a depth-based color tint to all rasterized pixels in the framebuffer.
    ///
    /// This is a post-pass that runs after all geometry has been rasterized.
    /// For each pixel with a valid depth (not `f64::INFINITY`), its color is
    /// lerped toward `fog_color` based on how far it is from the camera:
    ///
    /// - Nearest pixels (z == z_min) keep their original color
    /// - Farthest pixels (z == z_max) are blended most toward `fog_color`
    /// - The `fog_strength` parameter (0.0..=1.0) controls the maximum blend
    ///
    /// Background pixels (depth == INFINITY) remain unchanged (black).
    pub fn apply_depth_tint(&mut self, fog_color: [u8; 3], fog_strength: f64) {
        // Find z_min and z_max across all valid (non-background) pixels.
        let mut z_min = f64::INFINITY;
        let mut z_max = f64::NEG_INFINITY;
        for &d in &self.depth {
            if d < f64::INFINITY {
                if d < z_min { z_min = d; }
                if d > z_max { z_max = d; }
            }
        }

        // No valid pixels, or all at the same depth — nothing to tint.
        let z_range = z_max - z_min;
        if z_range.abs() < 1e-12 {
            return;
        }

        let inv_range = 1.0 / z_range;

        for i in 0..self.depth.len() {
            let d = self.depth[i];
            if d >= f64::INFINITY {
                continue; // background pixel — leave black
            }
            let t = ((d - z_min) * inv_range).clamp(0.0, 1.0);
            let blend = t * fog_strength;
            let c = &mut self.color[i];
            c[0] = (c[0] as f64 + (fog_color[0] as f64 - c[0] as f64) * blend).clamp(0.0, 255.0) as u8;
            c[1] = (c[1] as f64 + (fog_color[1] as f64 - c[1] as f64) * blend).clamp(0.0, 255.0) as u8;
            c[2] = (c[2] as f64 + (fog_color[2] as f64 - c[2] as f64) * blend).clamp(0.0, 255.0) as u8;
        }
    }

    /// Draw a 3D line with Bresenham's algorithm and z-interpolation.
    ///
    /// Useful for wireframe and ball-and-stick rendering modes.
    /// `p1` and `p2` are `[x, y, z]` in screen/pixel space.
    pub fn draw_line_3d(&mut self, p1: [f64; 3], p2: [f64; 3], color: [u8; 3]) {
        let mut x0 = p1[0].round() as isize;
        let mut y0 = p1[1].round() as isize;
        let x1 = p2[0].round() as isize;
        let y1 = p2[1].round() as isize;

        // Skip lines where both endpoints are entirely off the same side of the screen.
        let w = self.width as isize;
        let h = self.height as isize;
        if (x0 < 0 && x1 < 0) || (x0 >= w && x1 >= w) ||
           (y0 < 0 && y1 < 0) || (y0 >= h && y1 >= h) {
            return;
        }

        let dx = (x1 - x0).abs();
        let dy = -(y1 - y0).abs();
        let sx: isize = if x0 < x1 { 1 } else { -1 };
        let sy: isize = if y0 < y1 { 1 } else { -1 };
        let mut err = dx + dy;

        // Total Manhattan-ish distance for z interpolation
        let total_steps = dx.max(-dy) as f64;

        loop {
            // Compute interpolation parameter t
            let t = if total_steps > 0.0 {
                let from_start_x = (x0 - p1[0].round() as isize).unsigned_abs() as f64;
                let from_start_y = (y0 - p1[1].round() as isize).unsigned_abs() as f64;
                from_start_x.max(from_start_y) / total_steps
            } else {
                0.0
            };
            let z = p1[2] * (1.0 - t) + p2[2] * t;

            if x0 >= 0 && y0 >= 0 && (x0 as usize) < self.width && (y0 as usize) < self.height {
                self.set_pixel(x0 as usize, y0 as usize, z, color);
            }

            if x0 == x1 && y0 == y1 {
                break;
            }

            let e2 = 2 * err;
            if e2 >= dy {
                err += dy;
                x0 += sx;
            }
            if e2 <= dx {
                err += dx;
                y0 += sy;
            }
        }
    }

    /// Draw a 3D line with the given pixel thickness.
    ///
    /// For each pixel along the main Bresenham line, a filled circle of the
    /// given radius (`thickness / 2`) is drawn perpendicular to the line so
    /// that the resulting stroke has the requested width.  For `thickness`
    /// values <= 1.0 this falls back to the regular single-pixel `draw_line_3d`.
    pub fn draw_thick_line_3d(
        &mut self,
        p1: [f64; 3],
        p2: [f64; 3],
        color: [u8; 3],
        thickness: f64,
    ) {
        if thickness <= 1.0 {
            self.draw_line_3d(p1, p2, color);
            return;
        }

        let half = thickness / 2.0;

        // Direction vector in screen-space (xy only).
        let dx = p2[0] - p1[0];
        let dy = p2[1] - p1[1];
        let len = (dx * dx + dy * dy).sqrt();

        if len < 1e-6 {
            // Degenerate (zero-length) line – just draw a dot.
            self.draw_circle_z(p1[0], p1[1], p1[2], half, color);
            return;
        }

        // Perpendicular unit vector in screen-space.
        let px = -dy / len;
        let py = dx / len;

        // Draw the line at several perpendicular offsets to fill out the
        // thickness.  Step size of 0.5 gives smooth coverage.
        let steps = (half / 0.5).ceil() as isize;
        for i in -steps..=steps {
            let offset = i as f64 * 0.5;
            if offset * offset > half * half {
                continue;
            }
            let off_p1 = [p1[0] + px * offset, p1[1] + py * offset, p1[2]];
            let off_p2 = [p2[0] + px * offset, p2[1] + py * offset, p2[2]];
            self.draw_line_3d(off_p1, off_p2, color);
        }
    }

    /// Draw a filled circle at pixel coordinates `(cx, cy)` with the given radius and color.
    #[cfg(test)]
    pub fn draw_circle(&mut self, cx: f64, cy: f64, radius: f64, color: [u8; 3]) {
        self.draw_circle_z(cx, cy, 0.0, radius, color);
    }

    /// Convert this framebuffer's pixel data into an `image::RgbImage`.
    ///
    /// The resulting image has the same width and height as the framebuffer,
    /// with each pixel's RGB channels copied directly.  This is used by the
    /// ratatui-image integration to send the framebuffer to the terminal via
    /// Sixel, Kitty, or other graphics protocols.
    pub fn to_rgb_image(&self) -> RgbImage {
        let mut img = RgbImage::new(self.width as u32, self.height as u32);
        for y in 0..self.height {
            for x in 0..self.width {
                let c = self.color[y * self.width + x];
                img.put_pixel(x as u32, y as u32, image::Rgb(c));
            }
        }
        img
    }

    /// Convert this framebuffer into an `image::RgbaImage` with transparency.
    ///
    /// Background pixels (depth == INFINITY, color == black) get alpha = 0 so
    /// the terminal background shows through.  Drawn pixels get alpha = 255.
    pub fn to_rgba_image(&self) -> RgbaImage {
        let mut img = RgbaImage::new(self.width as u32, self.height as u32);
        for y in 0..self.height {
            for x in 0..self.width {
                let idx = y * self.width + x;
                let c = self.color[idx];
                let alpha = if self.depth[idx] >= f64::INFINITY { 0 } else { 255 };
                img.put_pixel(x as u32, y as u32, image::Rgba([c[0], c[1], c[2], alpha]));
            }
        }
        img
    }

    /// Draw a filled circle with a specific z-depth for z-buffer testing.
    ///
    /// Pixels at the center of the circle are written at `z`; pixels at the
    /// edge are also written at `z` (flat disk). For a sphere-like appearance,
    /// callers can submit multiple concentric circles with varying z.
    pub fn draw_circle_z(&mut self, cx: f64, cy: f64, z: f64, radius: f64, color: [u8; 3]) {
        let ix_min = ((cx - radius).floor() as isize).max(0) as usize;
        let ix_max = ((cx + radius).ceil() as isize).max(0).min(self.width as isize - 1) as usize;
        let iy_min = ((cy - radius).floor() as isize).max(0) as usize;
        let iy_max = ((cy + radius).ceil() as isize).max(0).min(self.height as isize - 1) as usize;

        // Circle entirely off-screen
        if ix_min > ix_max || iy_min > iy_max {
            return;
        }

        let r_sq = radius * radius;

        for py in iy_min..=iy_max {
            let dy = py as f64 + 0.5 - cy;
            let dy_sq = dy * dy;
            if dy_sq > r_sq {
                continue;
            }
            for px in ix_min..=ix_max {
                let dx = px as f64 + 0.5 - cx;
                if dx * dx + dy_sq <= r_sq {
                    self.set_pixel(px, py, z, color);
                }
            }
        }
    }
}

/// Quantize a single color channel by rounding to the nearest multiple of
/// `step`.  This reduces the number of distinct RGB triples in the output so
/// that more adjacent cells share the same color and get merged into longer
/// runs, dramatically cutting the number of ANSI escape sequences emitted per
/// frame.
///
/// A `step` of 1 is a no-op (full precision).  A `step` of 8 reduces 256
/// levels to 32 distinct values -- visually almost imperceptible but can cut
/// the span count (and therefore terminal output size) by 3-5x.
#[inline]
fn quantize_channel(v: u8, step: u8) -> u8 {
    if step <= 1 {
        return v;
    }
    let half = step / 2;
    // Round to nearest multiple of `step`, clamped to 255.
    let q = ((v as u16 + half as u16) / step as u16) * step as u16;
    q.min(255) as u8
}

/// Quantize an RGB triple.  Black `[0,0,0]` is kept exactly black so that the
/// blank-cell optimisation still fires.
#[inline]
fn quantize_color(c: [u8; 3], step: u8) -> [u8; 3] {
    if c == [0, 0, 0] {
        return c;
    }
    let q = [
        quantize_channel(c[0], step),
        quantize_channel(c[1], step),
        quantize_channel(c[2], step),
    ];
    // Avoid rounding a near-black color *to* black, which would make it
    // invisible.  Clamp to at least `step` in the brightest channel.
    if q == [0, 0, 0] {
        [step, step, step]
    } else {
        q
    }
}

/// Convert a [`Framebuffer`] into a ratatui [`Paragraph`] widget using half-block characters.
///
/// Each terminal row maps to two pixel rows:
/// - Top pixel  -> foreground color
/// - Bottom pixel -> background color
/// - Character: `'▀'` (upper half block, U+2580)
///
/// Consecutive cells with identical (fg, bg) pairs are merged into a single
/// [`Span`] to reduce the number of styled segments ratatui needs to process.
///
/// Color quantization is available via `quant_step` but currently disabled
/// (`quant_step = 1`, a no-op).  Set it to e.g. 4 or 8 to reduce distinct
/// colors and increase run-length merging at the expense of color precision.
#[cfg(test)]
pub fn framebuffer_to_widget(fb: &Framebuffer) -> Paragraph<'static> {
    // Quantization step: 4 gives 64 levels per channel -- preserves smooth
    // shading gradients while still merging runs.  Use 1 (no quantization)
    // for tiny framebuffers where output is already small.
    let quant_step: u8 = 1;

    let term_rows = (fb.height + 1) / 2;
    let mut lines: Vec<Line<'static>> = Vec::with_capacity(term_rows);

    for tr in 0..term_rows {
        let top_row = tr * 2;
        let bot_row = top_row + 1;

        let mut spans: Vec<Span<'static>> = Vec::new();

        // We classify each cell into one of 4 cases and track runs of same type
        // Case 0: blank (both black) → space, no styling (terminal bg shows through)
        // Case 1: both have color → '▀' with fg=top, bg=bottom
        // Case 2: top only → '▀' with fg=top, bg=Reset
        // Case 3: bottom only → '▄' with fg=bottom, bg=Reset
        #[derive(PartialEq, Clone, Copy)]
        enum CellKind { Blank, Both([u8;3],[u8;3]), TopOnly([u8;3]), BotOnly([u8;3]) }

        let mut run_text = String::new();
        let mut run_kind = CellKind::Blank;
        let mut run_started = false;

        let flush = |spans: &mut Vec<Span<'static>>, text: &str, kind: &CellKind| {
            if text.is_empty() { return; }
            let style = match kind {
                CellKind::Blank => Style::default(),
                CellKind::Both(top, bot) => Style::default()
                    .fg(Color::Rgb(top[0], top[1], top[2]))
                    .bg(Color::Rgb(bot[0], bot[1], bot[2])),
                CellKind::TopOnly(top) => Style::default()
                    .fg(Color::Rgb(top[0], top[1], top[2])),
                CellKind::BotOnly(bot) => Style::default()
                    .fg(Color::Rgb(bot[0], bot[1], bot[2])),
            };
            spans.push(Span::styled(text.to_string(), style));
        };

        for col in 0..fb.width {
            let top = quantize_color(fb.color[top_row * fb.width + col], quant_step);
            let bot = if bot_row < fb.height {
                quantize_color(fb.color[bot_row * fb.width + col], quant_step)
            } else {
                [0, 0, 0]
            };

            let top_black = top == [0, 0, 0];
            let bot_black = bot == [0, 0, 0];

            let kind = if top_black && bot_black {
                CellKind::Blank
            } else if !top_black && !bot_black {
                CellKind::Both(top, bot)
            } else if !top_black {
                CellKind::TopOnly(top)
            } else {
                CellKind::BotOnly(bot)
            };

            if run_started && kind == run_kind {
                match kind {
                    CellKind::Blank => run_text.push(' '),
                    CellKind::Both(..) | CellKind::TopOnly(_) => run_text.push('\u{2580}'),
                    CellKind::BotOnly(_) => run_text.push('\u{2584}'),
                }
            } else {
                if run_started {
                    flush(&mut spans, &run_text, &run_kind);
                }
                run_text.clear();
                match kind {
                    CellKind::Blank => run_text.push(' '),
                    CellKind::Both(..) | CellKind::TopOnly(_) => run_text.push('\u{2580}'),
                    CellKind::BotOnly(_) => run_text.push('\u{2584}'),
                }
                run_kind = kind;
                run_started = true;
            }
        }

        if run_started {
            flush(&mut spans, &run_text, &run_kind);
        }

        lines.push(Line::from(spans));
    }

    Paragraph::new(lines)
}

/// Convert a [`Framebuffer`] rendered at braille resolution into a ratatui
/// [`Paragraph`] widget using colored Unicode braille characters.
///
/// The framebuffer is expected to have dimensions `(cols * 2, rows * 4)` where
/// `cols` and `rows` are the target terminal cell dimensions.  Each terminal
/// cell maps to a 2x4 block of pixels.  Non-black pixels become "on" braille
/// dots; their average RGB color is used as the cell's foreground color.
///
/// This gives 4x the spatial resolution of half-block rendering at the cost of
/// per-cell (rather than per-pixel) coloring.
///
/// Consecutive cells with the same foreground color are merged into a single
/// [`Span`] for performance (run-length encoding).  Color quantization is
/// available via `quant_step` but currently disabled (`quant_step = 1`, a
/// no-op).  Set it to e.g. 4 or 8 to reduce distinct colors and increase
/// run-length merging at the expense of color precision.
pub fn framebuffer_to_braille_widget(fb: &Framebuffer) -> Paragraph<'static> {
    // Terminal cell grid dimensions derived from the framebuffer.
    let term_cols = (fb.width + 1) / 2;
    let term_rows = (fb.height + 3) / 4;

    if term_cols == 0 || term_rows == 0 {
        return Paragraph::new("");
    }

    // Braille dot bit values indexed by (dx, dy) within the 2x4 cell block.
    // Layout:
    //   Col 0  Col 1
    //   bit 0  bit 3   (row 0)
    //   bit 1  bit 4   (row 1)
    //   bit 2  bit 5   (row 2)
    //   bit 6  bit 7   (row 3)
    const BRAILLE_BITS: [[u8; 4]; 2] = [
        [0x01, 0x02, 0x04, 0x40], // column 0: rows 0-3
        [0x08, 0x10, 0x20, 0x80], // column 1: rows 0-3
    ];

    let quant_step: u8 = 1;

    let mut lines: Vec<Line<'static>> = Vec::with_capacity(term_rows);

    for tr in 0..term_rows {
        let py_base = tr * 4;

        let mut spans: Vec<Span<'static>> = Vec::new();
        let mut run_text = String::new();
        // Track the current run's color; None means blank (space) run.
        let mut run_color: Option<[u8; 3]> = None;
        let mut run_started = false;

        let flush =
            |spans: &mut Vec<Span<'static>>, text: &str, color: &Option<[u8; 3]>| {
                if text.is_empty() {
                    return;
                }
                let style = match color {
                    Some(c) => Style::default().fg(Color::Rgb(c[0], c[1], c[2])),
                    None => Style::default(),
                };
                spans.push(Span::styled(text.to_string(), style));
            };

        for tc in 0..term_cols {
            let px_base = tc * 2;

            // Build braille bit pattern and accumulate color of "on" dots.
            let mut bits: u8 = 0;
            let mut r_sum: u32 = 0;
            let mut g_sum: u32 = 0;
            let mut b_sum: u32 = 0;
            let mut on_count: u32 = 0;

            for dx in 0..2usize {
                let px = px_base + dx;
                if px >= fb.width {
                    continue;
                }
                for dy in 0..4usize {
                    let py = py_base + dy;
                    if py >= fb.height {
                        continue;
                    }
                    let c = fb.color[py * fb.width + px];
                    if c != [0, 0, 0] {
                        bits |= BRAILLE_BITS[dx][dy];
                        r_sum += c[0] as u32;
                        g_sum += c[1] as u32;
                        b_sum += c[2] as u32;
                        on_count += 1;
                    }
                }
            }

            if bits == 0 {
                // All dots off — emit a space.
                let cell_color: Option<[u8; 3]> = None;
                if run_started && run_color == cell_color {
                    run_text.push(' ');
                } else {
                    if run_started {
                        flush(&mut spans, &run_text, &run_color);
                    }
                    run_text.clear();
                    run_text.push(' ');
                    run_color = cell_color;
                    run_started = true;
                }
            } else {
                // Compute average color of "on" pixels.
                let avg = quantize_color(
                    [
                        (r_sum / on_count) as u8,
                        (g_sum / on_count) as u8,
                        (b_sum / on_count) as u8,
                    ],
                    quant_step,
                );
                let cell_color = Some(avg);

                let braille_char =
                    char::from_u32(0x2800u32 + bits as u32).unwrap_or(' ');

                if run_started && run_color == cell_color {
                    run_text.push(braille_char);
                } else {
                    if run_started {
                        flush(&mut spans, &run_text, &run_color);
                    }
                    run_text.clear();
                    run_text.push(braille_char);
                    run_color = cell_color;
                    run_started = true;
                }
            }
        }

        if run_started {
            flush(&mut spans, &run_text, &run_color);
        }

        lines.push(Line::from(spans));
    }

    Paragraph::new(lines)
}

/// Normalize a 3-component vector in place and return the result.
/// If the vector has zero length, returns `[0.0, 0.0, 0.0]`.
pub fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

/// Default light direction (normalized) pointing from the upper-right-front.
///
/// The X component is negative to compensate for the camera's X-negation
/// (which corrects chirality so L-amino acids render correctly).  In screen
/// space this still appears as upper-right lighting.
pub fn default_light_dir() -> [f64; 3] {
    normalize([-0.3, 0.8, 0.5])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_and_clear() {
        let mut fb = Framebuffer::new(4, 4);
        assert_eq!(fb.color.len(), 16);
        assert_eq!(fb.depth.len(), 16);
        assert!(fb.depth[0].is_infinite());
        assert_eq!(fb.color[0], [0, 0, 0]);

        // Write something, then clear
        fb.color[0] = [255, 0, 0];
        fb.depth[0] = 1.0;
        fb.clear();
        assert_eq!(fb.color[0], [0, 0, 0]);
        assert!(fb.depth[0].is_infinite());
    }

    #[test]
    fn test_zbuffer() {
        let mut fb = Framebuffer::new(4, 4);
        fb.set_pixel(1, 1, 5.0, [255, 0, 0]);
        assert_eq!(fb.color[1 * 4 + 1], [255, 0, 0]);

        // Closer fragment wins
        fb.set_pixel(1, 1, 3.0, [0, 255, 0]);
        assert_eq!(fb.color[1 * 4 + 1], [0, 255, 0]);

        // Farther fragment is rejected
        fb.set_pixel(1, 1, 4.0, [0, 0, 255]);
        assert_eq!(fb.color[1 * 4 + 1], [0, 255, 0]);
    }

    #[test]
    fn test_normalize() {
        let v = normalize([3.0, 0.0, 4.0]);
        let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-9);
        assert!((v[0] - 0.6).abs() < 1e-9);
        assert!((v[2] - 0.8).abs() < 1e-9);
    }

    #[test]
    fn test_default_light_dir_is_unit() {
        let d = default_light_dir();
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_rasterize_covers_pixels() {
        let mut fb = Framebuffer::new(10, 10);
        let tri = Triangle {
            verts: [[2.0, 2.0, 1.0], [8.0, 2.0, 1.0], [5.0, 8.0, 1.0]],
            color: [200, 100, 50],
            normal: [0.0, 0.0, 1.0],
        };
        fb.rasterize_triangle(&tri, normalize([0.0, 0.0, 1.0]));

        // The centroid (5,4) should definitely be filled
        let idx = 4 * 10 + 5;
        assert_ne!(fb.color[idx], [0, 0, 0]);
        // A corner outside the triangle should remain black
        assert_eq!(fb.color[0], [0, 0, 0]);
    }

    #[test]
    fn test_draw_line_3d() {
        let mut fb = Framebuffer::new(10, 10);
        fb.draw_line_3d([0.0, 0.0, 1.0], [9.0, 9.0, 2.0], [255, 255, 255]);
        // Diagonal line should touch (0,0) and (9,9)
        assert_eq!(fb.color[0], [255, 255, 255]);
        assert_eq!(fb.color[9 * 10 + 9], [255, 255, 255]);
    }

    #[test]
    fn test_draw_circle() {
        let mut fb = Framebuffer::new(20, 20);
        fb.draw_circle(10.0, 10.0, 3.0, [128, 64, 32]);
        // Center pixel should be filled
        assert_eq!(fb.color[10 * 20 + 10], [128, 64, 32]);
        // Far corner should not be
        assert_eq!(fb.color[0], [0, 0, 0]);
    }

    #[test]
    fn test_framebuffer_to_widget_basic() {
        let mut fb = Framebuffer::new(2, 4);
        fb.color[0] = [255, 0, 0]; // row 0, col 0 (top pixel of term row 0)
        fb.color[2] = [0, 255, 0]; // row 1, col 0 (bottom pixel of term row 0)
        // The widget should produce 2 terminal rows for 4 pixel rows
        let _widget = framebuffer_to_widget(&fb);
        // Just ensure it doesn't panic; visual inspection is manual
    }

    #[test]
    fn test_quantize_channel_no_op() {
        // step=1 should be a no-op
        assert_eq!(quantize_channel(0, 1), 0);
        assert_eq!(quantize_channel(127, 1), 127);
        assert_eq!(quantize_channel(255, 1), 255);
    }

    #[test]
    fn test_quantize_channel_step_8() {
        // 0 stays 0
        assert_eq!(quantize_channel(0, 8), 0);
        // 4 rounds to 8 (nearest multiple of 8)
        assert_eq!(quantize_channel(4, 8), 8);
        // 3 rounds to 0
        assert_eq!(quantize_channel(3, 8), 0);
        // 128 stays 128
        assert_eq!(quantize_channel(128, 8), 128);
        // 255 should quantize to a valid u8 (256 clamped to 255)
        let q255 = quantize_channel(255, 8);
        assert!(q255 == 248 || q255 == 255, "unexpected quantize(255, 8): {}", q255);
        // 252 should round up to 248 or be clamped
        let q252 = quantize_channel(252, 8);
        assert!(q252 == 248 || q252 == 255, "unexpected quantize(252, 8): {}", q252);
    }

    #[test]
    fn test_quantize_color_black_unchanged() {
        assert_eq!(quantize_color([0, 0, 0], 8), [0, 0, 0]);
    }

    #[test]
    fn test_quantize_color_near_black_stays_visible() {
        // A dim color like [3, 2, 1] would quantize to [0,0,0] without the
        // guard.  The function should keep it visible.
        let q = quantize_color([3, 2, 1], 8);
        assert_ne!(q, [0, 0, 0], "near-black color should not vanish");
    }

    #[test]
    fn test_quantize_color_normal() {
        let q = quantize_color([200, 100, 50], 8);
        // Each channel should be a multiple of 8 (or 255 if clamped)
        for &c in &q {
            assert!(c % 8 == 0 || c == 255, "channel {} not quantized", c);
        }
    }

    #[test]
    fn test_apply_depth_tint_blends_colors() {
        let mut fb = Framebuffer::new(4, 1);
        // Place pixels at different depths: near (z=1) and far (z=10)
        let near_color = [200, 100, 50];
        let far_color = [200, 100, 50];
        fb.color[0] = near_color;
        fb.depth[0] = 1.0;
        fb.color[1] = far_color;
        fb.depth[1] = 10.0;
        // Pixels 2 and 3 remain at INFINITY (background)

        let fog = [40, 50, 70];
        fb.apply_depth_tint(fog, 0.5);

        // Near pixel (z=1, t=0.0) should stay unchanged
        assert_eq!(fb.color[0], near_color, "nearest pixel should keep original color");

        // Far pixel (z=10, t=1.0) should be blended halfway toward fog
        // new = base + (fog - base) * 1.0 * 0.5
        // R: 200 + (40 - 200) * 0.5 = 200 - 80 = 120
        // G: 100 + (50 - 100) * 0.5 = 100 - 25 = 75
        // B:  50 + (70 -  50) * 0.5 =  50 + 10 = 60
        assert_eq!(fb.color[1], [120, 75, 60], "farthest pixel should blend toward fog");
    }

    #[test]
    fn test_apply_depth_tint_skips_background() {
        let mut fb = Framebuffer::new(4, 1);
        // Set one valid pixel and leave others at INFINITY
        fb.color[0] = [200, 100, 50];
        fb.depth[0] = 5.0;
        fb.color[1] = [180, 90, 40];
        fb.depth[1] = 10.0;
        // Pixels 2 and 3 are background (depth = INFINITY, color = [0,0,0])

        fb.apply_depth_tint([40, 50, 70], 0.5);

        // Background pixels must remain [0,0,0]
        assert_eq!(fb.color[2], [0, 0, 0], "background pixel at index 2 should stay black");
        assert_eq!(fb.color[3], [0, 0, 0], "background pixel at index 3 should stay black");
    }
}

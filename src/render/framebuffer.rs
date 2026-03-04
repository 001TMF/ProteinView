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
    pub fn new(width: usize, height: usize) -> Self {
        let n = width * height;
        Self {
            width,
            height,
            color: vec![[0, 0, 0]; n],
            depth: vec![f64::INFINITY; n],
        }
    }

    /// Reset the framebuffer to black pixels and infinite depth.
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
    pub fn rasterize_triangle(&mut self, tri: &Triangle, light_dir: [f64; 3]) {
        const AMBIENT: f64 = 0.15;

        // --- Lambert shading intensity ---
        let dot = tri.normal[0] * light_dir[0]
            + tri.normal[1] * light_dir[1]
            + tri.normal[2] * light_dir[2];
        let intensity = if dot > AMBIENT { dot } else { AMBIENT };

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
        let max_x = (max_x as usize).min(self.width.saturating_sub(1));
        let min_y = min_y.max(0) as usize;
        let max_y = (max_y as usize).min(self.height.saturating_sub(1));

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

    /// Draw a 3D line with Bresenham's algorithm and z-interpolation.
    ///
    /// Useful for wireframe and ball-and-stick rendering modes.
    /// `p1` and `p2` are `[x, y, z]` in screen/pixel space.
    pub fn draw_line_3d(&mut self, p1: [f64; 3], p2: [f64; 3], color: [u8; 3]) {
        let mut x0 = p1[0].round() as isize;
        let mut y0 = p1[1].round() as isize;
        let x1 = p2[0].round() as isize;
        let y1 = p2[1].round() as isize;

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

    /// Draw a filled circle at pixel coordinates `(cx, cy)` with the given radius and color.
    ///
    /// Used for rendering atoms as dots in ball-and-stick mode.
    /// No z-buffering on the circle itself (it writes at z = 0); callers should
    /// provide screen-depth by using `draw_circle_z` if depth is needed.
    pub fn draw_circle(&mut self, cx: f64, cy: f64, radius: f64, color: [u8; 3]) {
        self.draw_circle_z(cx, cy, 0.0, radius, color);
    }

    /// Draw a filled circle with a specific z-depth for z-buffer testing.
    ///
    /// Pixels at the center of the circle are written at `z`; pixels at the
    /// edge are also written at `z` (flat disk). For a sphere-like appearance,
    /// callers can submit multiple concentric circles with varying z.
    pub fn draw_circle_z(&mut self, cx: f64, cy: f64, z: f64, radius: f64, color: [u8; 3]) {
        let ix_min = ((cx - radius).floor() as isize).max(0) as usize;
        let ix_max = ((cx + radius).ceil() as isize).min(self.width as isize - 1) as usize;
        let iy_min = ((cy - radius).floor() as isize).max(0) as usize;
        let iy_max = ((cy + radius).ceil() as isize).min(self.height as isize - 1) as usize;

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

/// Convert a [`Framebuffer`] into a ratatui [`Paragraph`] widget using half-block characters.
///
/// Each terminal row maps to two pixel rows:
/// - Top pixel  -> foreground color
/// - Bottom pixel -> background color
/// - Character: `'▀'` (upper half block, U+2580)
///
/// Consecutive cells with identical (fg, bg) pairs are merged into a single
/// [`Span`] to reduce the number of styled segments ratatui needs to process.
pub fn framebuffer_to_widget(fb: &Framebuffer) -> Paragraph<'static> {
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
            let top = fb.color[top_row * fb.width + col];
            let bot = if bot_row < fb.height {
                fb.color[bot_row * fb.width + col]
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
pub fn default_light_dir() -> [f64; 3] {
    normalize([0.3, 0.8, 0.5])
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
}

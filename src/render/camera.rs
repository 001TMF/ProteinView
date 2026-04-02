use std::time::Instant;

/// 3D camera for viewing protein structures
#[derive(Debug, Clone)]
pub struct Camera {
    pub rot_x: f64,
    pub rot_y: f64,
    pub rot_z: f64,
    pub zoom: f64,
    pub pan_x: f64,
    pub pan_y: f64,
    pub auto_rotate: bool,
    last_tick: Instant,
}

impl Default for Camera {
    fn default() -> Self {
        Self {
            rot_x: 0.0,
            rot_y: 0.0,
            rot_z: 0.0,
            zoom: 1.0,
            pan_x: 0.0,
            pan_y: 0.0,
            auto_rotate: false,
            last_tick: Instant::now(),
        }
    }
}

/// A projected 2D point with depth
#[derive(Debug, Clone, Copy)]
pub struct Projected {
    pub x: f64,
    pub y: f64,
    pub z: f64, // depth for z-buffering
}

/// Pre-computed sin/cos values for a single frame.
///
/// Computing `sin_cos()` for each of the three rotation angles is expensive
/// when done per-vertex (millions of calls per frame).  `ProjectionCache`
/// evaluates the six trig values once and then reuses them for every
/// `project()` and `rotate_normal()` call within the frame.
#[derive(Debug, Clone, Copy)]
pub struct ProjectionCache {
    sin_x: f64,
    cos_x: f64,
    sin_y: f64,
    cos_y: f64,
    sin_z: f64,
    cos_z: f64,
    zoom: f64,
    pan_x: f64,
    pan_y: f64,
}

impl ProjectionCache {
    /// Project a 3D point to 2D using the cached trig values.
    ///
    /// The math is identical to [`Camera::project()`] but avoids redundant
    /// `sin_cos()` evaluations.
    #[inline]
    pub fn project(&self, x: f64, y: f64, z: f64) -> Projected {
        // Rotation around X axis
        let y1 = y * self.cos_x - z * self.sin_x;
        let z1 = y * self.sin_x + z * self.cos_x;

        // Rotation around Y axis
        let x2 = x * self.cos_y + z1 * self.sin_y;
        let z2 = -x * self.sin_y + z1 * self.cos_y;

        // Rotation around Z axis
        let x3 = x2 * self.cos_z - y1 * self.sin_z;
        let y3 = x2 * self.sin_z + y1 * self.cos_z;

        Projected {
            x: -x3 * self.zoom + self.pan_x,
            y: y3 * self.zoom + self.pan_y,
            z: z2,
        }
    }

    /// Apply the camera rotation to a direction vector (no zoom/pan).
    #[inline]
    pub fn rotate_normal(&self, nx: f64, ny: f64, nz: f64) -> [f64; 3] {
        let y1 = ny * self.cos_x - nz * self.sin_x;
        let z1 = ny * self.sin_x + nz * self.cos_x;

        let x2 = nx * self.cos_y + z1 * self.sin_y;
        let z2 = -nx * self.sin_y + z1 * self.cos_y;

        let x3 = x2 * self.cos_z - y1 * self.sin_z;
        let y3 = x2 * self.sin_z + y1 * self.cos_z;

        [x3, y3, z2]
    }
}

impl Camera {
    const ROT_STEP: f64 = 0.1;
    const ZOOM_STEP: f64 = 0.1;
    const PAN_STEP: f64 = 2.0;

    pub fn rotate_x(&mut self, dir: f64) {
        self.rot_x += dir * Self::ROT_STEP;
    }
    pub fn rotate_y(&mut self, dir: f64) {
        self.rot_y -= dir * Self::ROT_STEP;
    }
    pub fn rotate_z(&mut self, dir: f64) {
        self.rot_z -= dir * Self::ROT_STEP;
    }
    pub fn zoom_in(&mut self) {
        self.zoom *= 1.0 + Self::ZOOM_STEP;
    }
    pub fn zoom_out(&mut self) {
        self.zoom *= 1.0 - Self::ZOOM_STEP;
    }
    pub fn pan(&mut self, dx: f64, dy: f64) {
        self.pan_x += dx * Self::PAN_STEP;
        self.pan_y += dy * Self::PAN_STEP;
    }
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    /// Auto-rotate speed in radians per second (~0.6 rad/s = one full turn in ~10s).
    const AUTO_ROTATE_SPEED: f64 = 0.6;

    /// Maximum dt (in seconds) that a single tick can apply.  This prevents
    /// the protein from "jumping" when a frame takes longer than expected or
    /// when frames are skipped.  At 30 FPS the nominal interval is ~0.033s;
    /// we allow up to 2x that to accommodate occasional slow frames while
    /// still clamping large gaps.
    const MAX_DT: f64 = 0.066;

    pub fn tick(&mut self) {
        let now = Instant::now();
        let raw_dt = now.duration_since(self.last_tick).as_secs_f64();
        self.last_tick = now;
        // Clamp dt so that a long gap (frame skip, slow draw, debugger pause)
        // never produces a visible jump in auto-rotation.
        let dt = raw_dt.min(Self::MAX_DT);
        if self.auto_rotate {
            self.rot_y -= Self::AUTO_ROTATE_SPEED * dt;
        }
    }

    /// Reset the internal tick timer without applying any rotation.
    /// Call this when skipping frames so the next real tick starts from a
    /// fresh baseline rather than accumulating all the skipped time.
    pub fn reset_tick_timer(&mut self) {
        self.last_tick = Instant::now();
    }

    /// Project a 3D point to 2D using rotation matrices + orthographic projection
    pub fn project(&self, x: f64, y: f64, z: f64) -> Projected {
        // Rotation around X axis
        let (sin_x, cos_x) = self.rot_x.sin_cos();
        let y1 = y * cos_x - z * sin_x;
        let z1 = y * sin_x + z * cos_x;

        // Rotation around Y axis
        let (sin_y, cos_y) = self.rot_y.sin_cos();
        let x2 = x * cos_y + z1 * sin_y;
        let z2 = -x * sin_y + z1 * cos_y;

        // Rotation around Z axis
        let (sin_z, cos_z) = self.rot_z.sin_cos();
        let x3 = x2 * cos_z - y1 * sin_z;
        let y3 = x2 * sin_z + y1 * cos_z;

        // Apply zoom and pan (orthographic projection)
        Projected {
            x: -x3 * self.zoom + self.pan_x,
            y: y3 * self.zoom + self.pan_y,
            z: z2,
        }
    }

    /// Pre-compute sin/cos for the current rotation angles.
    ///
    /// Call this once per frame and then use [`ProjectionCache::project()`]
    /// and [`ProjectionCache::rotate_normal()`] for all per-vertex work.
    pub fn projection_cache(&self) -> ProjectionCache {
        let (sin_x, cos_x) = self.rot_x.sin_cos();
        let (sin_y, cos_y) = self.rot_y.sin_cos();
        let (sin_z, cos_z) = self.rot_z.sin_cos();
        ProjectionCache {
            sin_x,
            cos_x,
            sin_y,
            cos_y,
            sin_z,
            cos_z,
            zoom: self.zoom,
            pan_x: self.pan_x,
            pan_y: self.pan_y,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn project_identity_negates_x() {
        // With identity rotation (all angles zero) and unit zoom,
        // a point at positive world-X should project to negative screen-X.
        // This ensures a right-handed coordinate system so L-amino acids
        // are not rendered as their mirror-image D-amino acids.
        let cam = Camera::default();
        let p = cam.project(1.0, 0.0, 0.0);
        assert!(
            p.x < 0.0,
            "positive world-X should project to negative screen-X, got {}",
            p.x
        );
        assert!(
            (p.y).abs() < 1e-12,
            "Y should be zero for a point on the X axis"
        );
    }

    #[test]
    fn project_identity_preserves_y() {
        // Y axis should pass through without negation.
        let cam = Camera::default();
        let p = cam.project(0.0, 1.0, 0.0);
        assert!(
            p.y > 0.0,
            "positive world-Y should project to positive screen-Y, got {}",
            p.y
        );
        assert!(
            (p.x).abs() < 1e-12,
            "X should be zero for a point on the Y axis"
        );
    }

    #[test]
    fn project_respects_zoom_and_pan() {
        let mut cam = Camera::default();
        cam.zoom = 2.0;
        cam.pan_x = 5.0;
        cam.pan_y = 3.0;
        let p = cam.project(1.0, 1.0, 0.0);
        // x = -1.0 * 2.0 + 5.0 = 3.0
        assert!((p.x - 3.0).abs() < 1e-12, "expected x=3.0, got {}", p.x);
        // y = 1.0 * 2.0 + 3.0 = 5.0
        assert!((p.y - 5.0).abs() < 1e-12, "expected y=5.0, got {}", p.y);
    }

    #[test]
    fn rotate_y_produces_negative_delta() {
        let mut cam = Camera::default();
        cam.rotate_y(1.0);
        assert!(
            cam.rot_y < 0.0,
            "rotate_y(+1) should decrease rot_y, got {}",
            cam.rot_y
        );
    }

    #[test]
    fn rotate_z_produces_negative_delta() {
        let mut cam = Camera::default();
        cam.rotate_z(1.0);
        assert!(
            cam.rot_z < 0.0,
            "rotate_z(+1) should decrease rot_z, got {}",
            cam.rot_z
        );
    }

    #[test]
    fn tick_clamps_large_dt() {
        // Simulate a long gap (e.g. frame skip) by creating a camera whose
        // last_tick is far in the past, then calling tick().  The rotation
        // should be clamped to MAX_DT worth of movement.
        let mut cam = Camera::default();
        cam.auto_rotate = true;
        // Manually set last_tick 500ms in the past (way more than MAX_DT)
        cam.last_tick = Instant::now() - std::time::Duration::from_millis(500);
        cam.tick();
        // Expected rotation ≈ -AUTO_ROTATE_SPEED * MAX_DT = -0.6 * 0.066 = -0.0396
        // Without clamping it would be -0.6 * 0.5 = -0.3
        let expected_max = Camera::AUTO_ROTATE_SPEED * Camera::MAX_DT;
        assert!(
            cam.rot_y.abs() <= expected_max + 0.001,
            "rotation should be clamped to at most {:.4} rad, got {:.4}",
            expected_max,
            cam.rot_y.abs()
        );
    }

    #[test]
    fn reset_tick_timer_prevents_jump() {
        // After reset_tick_timer(), the next tick() should see near-zero dt
        // and apply negligible rotation.
        let mut cam = Camera::default();
        cam.auto_rotate = true;
        // Set last_tick far in the past
        cam.last_tick = Instant::now() - std::time::Duration::from_secs(2);
        // Reset timer (as the main loop does during frame skips)
        cam.reset_tick_timer();
        // Immediately tick — dt should be ~0
        cam.tick();
        assert!(
            cam.rot_y.abs() < 0.001,
            "rotation after reset_tick_timer + immediate tick should be ~0, got {}",
            cam.rot_y
        );
    }

    #[test]
    fn projection_cache_matches_camera_project() {
        // ProjectionCache::project() must produce identical results to
        // Camera::project() for any rotation/zoom/pan combination.
        let mut cam = Camera::default();
        cam.rot_x = 0.7;
        cam.rot_y = -1.2;
        cam.rot_z = 0.3;
        cam.zoom = 2.5;
        cam.pan_x = 10.0;
        cam.pan_y = -5.0;

        let cache = cam.projection_cache();
        let points = [
            (1.0, 2.0, 3.0),
            (-4.0, 0.5, -1.0),
            (0.0, 0.0, 0.0),
            (100.0, -200.0, 50.0),
        ];
        for (x, y, z) in points {
            let a = cam.project(x, y, z);
            let b = cache.project(x, y, z);
            assert!((a.x - b.x).abs() < 1e-12, "x mismatch for ({x},{y},{z})");
            assert!((a.y - b.y).abs() < 1e-12, "y mismatch for ({x},{y},{z})");
            assert!((a.z - b.z).abs() < 1e-12, "z mismatch for ({x},{y},{z})");
        }
    }

    #[test]
    fn projection_cache_rotate_normal_matches() {
        // Verify rotate_normal matches the standalone rotate_normal logic.
        let mut cam = Camera::default();
        cam.rot_x = 0.5;
        cam.rot_y = -0.8;
        cam.rot_z = 1.1;

        let cache = cam.projection_cache();

        // Manually compute the expected rotation (same math as Camera::project
        // but without zoom/pan, and without X negation).
        let normals = [
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.577, 0.577, 0.577),
        ];
        for (nx, ny, nz) in normals {
            let result = cache.rotate_normal(nx, ny, nz);

            // Reproduce the rotation manually
            let (sin_x, cos_x) = cam.rot_x.sin_cos();
            let y1 = ny * cos_x - nz * sin_x;
            let z1 = ny * sin_x + nz * cos_x;
            let (sin_y, cos_y) = cam.rot_y.sin_cos();
            let x2 = nx * cos_y + z1 * sin_y;
            let z2 = -nx * sin_y + z1 * cos_y;
            let (sin_z, cos_z) = cam.rot_z.sin_cos();
            let x3 = x2 * cos_z - y1 * sin_z;
            let y3 = x2 * sin_z + y1 * cos_z;

            assert!((result[0] - x3).abs() < 1e-12, "nx mismatch");
            assert!((result[1] - y3).abs() < 1e-12, "ny mismatch");
            assert!((result[2] - z2).abs() < 1e-12, "nz mismatch");
        }
    }
}

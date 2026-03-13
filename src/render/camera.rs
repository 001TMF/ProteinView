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

impl Camera {
    const ROT_STEP: f64 = 0.1;
    const ZOOM_STEP: f64 = 0.1;
    const PAN_STEP: f64 = 2.0;

    pub fn rotate_x(&mut self, dir: f64) { self.rot_x += dir * Self::ROT_STEP; }
    pub fn rotate_y(&mut self, dir: f64) { self.rot_y -= dir * Self::ROT_STEP; }
    pub fn rotate_z(&mut self, dir: f64) { self.rot_z -= dir * Self::ROT_STEP; }
    pub fn zoom_in(&mut self) { self.zoom *= 1.0 + Self::ZOOM_STEP; }
    pub fn zoom_out(&mut self) { self.zoom *= 1.0 - Self::ZOOM_STEP; }
    pub fn pan(&mut self, dx: f64, dy: f64) { self.pan_x += dx * Self::PAN_STEP; self.pan_y += dy * Self::PAN_STEP; }
    pub fn reset(&mut self) { *self = Self::default(); }

    /// Auto-rotate speed in radians per second (~0.6 rad/s = one full turn in ~10s).
    const AUTO_ROTATE_SPEED: f64 = 0.6;

    pub fn tick(&mut self) {
        let now = Instant::now();
        let dt = now.duration_since(self.last_tick).as_secs_f64();
        self.last_tick = now;
        if self.auto_rotate {
            self.rot_y -= Self::AUTO_ROTATE_SPEED * dt;
        }
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
        assert!((p.y).abs() < 1e-12, "Y should be zero for a point on the X axis");
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
        assert!((p.x).abs() < 1e-12, "X should be zero for a point on the Y axis");
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
        assert!(cam.rot_y < 0.0, "rotate_y(+1) should decrease rot_y, got {}", cam.rot_y);
    }

    #[test]
    fn rotate_z_produces_negative_delta() {
        let mut cam = Camera::default();
        cam.rotate_z(1.0);
        assert!(cam.rot_z < 0.0, "rotate_z(+1) should decrease rot_z, got {}", cam.rot_z);
    }
}

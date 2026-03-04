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
    pub fn rotate_y(&mut self, dir: f64) { self.rot_y += dir * Self::ROT_STEP; }
    pub fn rotate_z(&mut self, dir: f64) { self.rot_z += dir * Self::ROT_STEP; }
    pub fn zoom_in(&mut self) { self.zoom *= 1.0 + Self::ZOOM_STEP; }
    pub fn zoom_out(&mut self) { self.zoom *= 1.0 - Self::ZOOM_STEP; }
    pub fn pan(&mut self, dx: f64, dy: f64) { self.pan_x += dx * Self::PAN_STEP; self.pan_y += dy * Self::PAN_STEP; }
    pub fn reset(&mut self) { *self = Self::default(); }

    pub fn tick(&mut self) {
        if self.auto_rotate {
            self.rot_y += 0.02;
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
            x: x3 * self.zoom + self.pan_x,
            y: y3 * self.zoom + self.pan_y,
            z: z2,
        }
    }
}

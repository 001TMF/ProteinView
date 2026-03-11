use std::time::Instant;

/// 3D camera for viewing protein structures
#[derive(Debug, Clone)]
pub struct Camera {
    orientation: [[f64; 3]; 3],
    pub zoom: f64,
    pub pan_x: f64,
    pub pan_y: f64,
    pivot: [f64; 3],
    pub auto_rotate: bool,
    last_tick: Instant,
}

impl Default for Camera {
    fn default() -> Self {
        Self {
            orientation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            zoom: 1.0,
            pan_x: 0.0,
            pan_y: 0.0,
            pivot: [0.0, 0.0, 0.0],
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

    pub fn rotate_x(&mut self, dir: f64) {
        self.apply_local_rotation([1.0, 0.0, 0.0], dir * Self::ROT_STEP);
    }
    pub fn rotate_y(&mut self, dir: f64) {
        self.apply_local_rotation([0.0, 1.0, 0.0], dir * Self::ROT_STEP);
    }
    pub fn rotate_z(&mut self, dir: f64) {
        self.apply_local_rotation([0.0, 0.0, 1.0], dir * Self::ROT_STEP);
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
    pub fn set_pivot(&mut self, pivot: [f64; 3]) {
        self.pivot = pivot;
    }
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    /// Auto-rotate speed in radians per second (~0.6 rad/s = one full turn in ~10s).
    const AUTO_ROTATE_SPEED: f64 = 0.6;

    pub fn tick(&mut self) {
        let now = Instant::now();
        let dt = now.duration_since(self.last_tick).as_secs_f64();
        self.last_tick = now;
        if self.auto_rotate {
            self.apply_local_rotation([0.0, 1.0, 0.0], -Self::AUTO_ROTATE_SPEED * dt);
        }
    }

    fn apply_local_rotation(&mut self, axis: [f64; 3], angle: f64) {
        let rotation = rotation_matrix(axis, angle);
        self.orientation = mat_mul(rotation, self.orientation);
    }

    pub fn rotate_vector(&self, x: f64, y: f64, z: f64) -> [f64; 3] {
        mat_vec_mul(self.orientation, [x, y, z])
    }

    /// Project a 3D point to 2D using the current camera orientation + orthographic projection.
    pub fn project(&self, x: f64, y: f64, z: f64) -> Projected {
        let [x_view, y_view, z_view] =
            self.rotate_vector(x - self.pivot[0], y - self.pivot[1], z - self.pivot[2]);

        // Apply zoom and pan (orthographic projection)
        Projected {
            // Flip screen-space X so the rendered scene is not mirrored.
            x: -x_view * self.zoom + self.pan_x,
            y: y_view * self.zoom + self.pan_y,
            z: z_view,
        }
    }
}

fn rotation_matrix(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let [x, y, z] = axis;
    let (sin_a, cos_a) = angle.sin_cos();
    let one_minus_cos = 1.0 - cos_a;

    [
        [
            cos_a + x * x * one_minus_cos,
            x * y * one_minus_cos - z * sin_a,
            x * z * one_minus_cos + y * sin_a,
        ],
        [
            y * x * one_minus_cos + z * sin_a,
            cos_a + y * y * one_minus_cos,
            y * z * one_minus_cos - x * sin_a,
        ],
        [
            z * x * one_minus_cos - y * sin_a,
            z * y * one_minus_cos + x * sin_a,
            cos_a + z * z * one_minus_cos,
        ],
    ]
}

fn mat_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for row in 0..3 {
        for col in 0..3 {
            out[row][col] = a[row][0] * b[0][col] + a[row][1] * b[1][col] + a[row][2] * b[2][col];
        }
    }
    out
}

fn mat_vec_mul(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(lhs: [f64; 3], rhs: [f64; 3]) {
        for i in 0..3 {
            assert!(
                (lhs[i] - rhs[i]).abs() < 1e-9,
                "index {i}: left={} right={}",
                lhs[i],
                rhs[i]
            );
        }
    }

    #[test]
    fn rotate_x_uses_current_view_basis() {
        let mut camera = Camera::default();
        camera.apply_local_rotation([0.0, 1.0, 0.0], std::f64::consts::FRAC_PI_2);
        camera.apply_local_rotation([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_2);

        let rotated = camera.rotate_vector(1.0, 1.0, 0.0);
        let expected = mat_vec_mul(
            mat_mul(
                rotation_matrix([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_2),
                rotation_matrix([0.0, 1.0, 0.0], std::f64::consts::FRAC_PI_2),
            ),
            [1.0, 1.0, 0.0],
        );

        assert_close(rotated, expected);
    }

    #[test]
    fn project_uses_pivot_as_rotation_center() {
        let mut camera = Camera::default();
        camera.set_pivot([5.0, 0.0, 0.0]);

        let projected = camera.project(6.0, 0.0, 0.0);
        assert!((projected.x + 1.0).abs() < 1e-9);
        assert!(projected.y.abs() < 1e-9);
    }
}

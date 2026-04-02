# Research: Mirroring Bug Fix

## 1. Root Cause Analysis

### Current projection in `src/render/camera.rs` (lines 65-87)

The `project()` method applies three sequential Euler angle rotations (Rx, Ry, Rz) then outputs:

```rust
Projected {
    x: x3 * self.zoom + self.pan_x,   // x3 is the rotated X component
    y: y3 * self.zoom + self.pan_y,    // y3 is the rotated Y component
    z: z2,                              // depth
}
```

**The problem:** The standard Euler rotation matrices used here produce a right-handed coordinate system where +X points right when viewed from behind the Z-axis. But screen-space convention (and viewer expectation) is to look *down* the -Z axis, which means the X-axis appears mirrored -- left and right are swapped. This is equivalent to viewing the scene through a mirror.

In concrete terms: at identity rotation (rot_x = rot_y = rot_z = 0), a point at world position (+1, 0, 0) projects to screen position (+1, 0) -- it appears on the right. But the viewer is implicitly looking along +Z (since there is no explicit view transform that flips the Z axis). In a standard right-handed "look-at" camera looking along -Z, that same point should project to screen position (-1, 0) -- it should appear on the *left* side, or equivalently the X axis of the output needs to be negated.

**The mirrored axis is X.** The projection is essentially doing a right-handed rotation but displaying the result with the wrong handedness convention, causing a horizontal mirror.

### Mathematical verification

At identity rotation (all angles = 0): sin=0, cos=1 for all axes.
- Rx: y1=y, z1=z
- Ry: x2=x, z2=z
- Rz: x3=x, y3=y

So `project(1, 0, 0)` yields `Projected { x: 1.0, y: 0.0, z: 0.0 }`. The `to_pixel` function in hd.rs (line 97) then maps this to `[1.0 + half_w, half_h - 0.0, 0.0]`, placing it to the right of center. A point at world +X appearing screen-right is correct for a camera looking along +Z, but protein viewers conventionally look along -Z (or equivalently, the +Z axis points toward the viewer). This flips the handedness and makes L-amino acids appear as D-amino acids.

## 2. The Fork's Fix (branch `desperadus-master`)

### What the fork did

The fork (`desperadus-master`) did a complete camera rewrite:
- Replaced Euler angles (`rot_x`, `rot_y`, `rot_z`) with an orientation matrix (`[[f64; 3]; 3]`)
- Added Rodrigues' rotation formula via `rotation_matrix()` for axis-angle rotation
- Added matrix multiplication infrastructure (`mat_mul`, `mat_vec_mul`)
- Added a pivot point for rotation center

### The actual mirroring fix

Despite the large rewrite, the mirroring fix itself is a single line in `project()`:

```rust
Projected {
    x: -x_view * self.zoom + self.pan_x,   // NOTE: negated x_view
    y: y_view * self.zoom + self.pan_y,
    z: z_view,
}
```

The key change is **negating the X component** of the projection output: `-x_view` instead of `x_view`.

### Minimal fix for our Euler camera

In our current `project()` method (line 82-86), the fix is simply:

```rust
Projected {
    x: -x3 * self.zoom + self.pan_x,   // was: x3 * self.zoom + self.pan_x
    y: y3 * self.zoom + self.pan_y,
    z: z2,
}
```

That is one character change: add a `-` before `x3`.

### Compensating changes the fork also made

The fork made several other changes that compensate for side effects of the X negation:

1. **Auto-rotate direction** (camera.rs `tick()`): Changed from `self.rot_y += ...` to `self.apply_local_rotation([0.0, 1.0, 0.0], -Self::AUTO_ROTATE_SPEED * dt)` -- note the negated angle. In our Euler camera this would mean negating the auto-rotate increment.

2. **j/k rotation direction** (main.rs): Swapped the direction arguments:
   - `j`/Down: changed from `rotate_x(1.0)` to `rotate_x(-1.0)`
   - `k`/Up: changed from `rotate_x(-1.0)` to `rotate_x(1.0)`

3. **rotate_normal in hd.rs**: Changed from inline Euler rotation to `camera.rotate_vector(nx, ny, nz)` which uses the orientation matrix. The normal rotation itself does NOT negate X, because the normal just needs to be rotated, not projected to screen space.

4. **Pan direction**: The fork did NOT change pan direction -- `pan(-1.0, 0.0)` for 'a' and `pan(1.0, 0.0)` for 'd' remain the same.

## 3. Downstream Effects Analysis

### 3.1 `src/render/hd.rs` -- Normal rotation and lighting

**`rotate_normal()` (lines 78-92):** This function applies the same Euler rotation to surface normals for Lambert shading. It does NOT include zoom/pan, and crucially it does NOT include the X negation from `project()`. This is correct behavior -- normals should be rotated by the same rotation as the geometry, but should not be mirrored.

**After applying the fix:** `rotate_normal()` should remain UNCHANGED. The normal rotation is independent of the screen-space X flip. The light direction (`default_light_dir() = normalize([0.3, 0.8, 0.5])`) is defined in view space (after rotation), so the dot product with rotated normals will be the same regardless of the X flip in the final projection.

**`to_pixel()` (line 96-98):** Converts from centered projection coordinates to pixel coordinates:
```rust
[proj_x + half_w, half_h - proj_y, proj_z]
```
This adds half_w to X and flips Y (since screen Y goes down, but world Y goes up). This function does NOT need to change -- the X negation happens upstream in `project()`.

**Risk: LOW.** No changes needed in hd.rs.

### 3.2 `src/render/braille.rs` -- X-axis assumptions

The braille renderer uses `camera.project()` to get 2D coordinates and passes them directly to the ratatui Canvas widget. The Canvas has its own coordinate system with `x_bounds` set to `[-width/2, width/2]`, so projections are already in centered coordinates.

**After applying the fix:** The braille renderer will automatically render the correct (un-mirrored) image because it relies entirely on `project()` for its coordinates. No X-axis assumptions in the code.

**Risk: NONE.** No changes needed.

### 3.3 `src/render/ribbon.rs` -- Triangle winding order

**THIS IS THE PRIMARY RISK AREA.**

The ribbon mesh is generated in world space. Triangle winding order determines which side of a triangle is the "front" face. The winding order is set by `emit_strip()` (lines 215-252) and `emit_quad()` (lines 868-888).

In `emit_strip()`, for each quad formed by two consecutive cross-section rings:
```rust
// Triangle 1: a0, a1, b0
// Triangle 2: a1, b1, b0
```

The `triangle_normal()` function (lines 203-207) computes the normal via cross product:
```rust
fn triangle_normal(v0: V3, v1: V3, v2: V3) -> V3 {
    let e1 = v3_sub(v1, v0);
    let e2 = v3_sub(v2, v0);
    v3_normalize(v3_cross(e1, e2))
}
```

**Critical analysis:** When we negate X in the projection, the screen-space winding of triangles reverses (what was CW becomes CCW and vice versa). However, our rasterizer uses **two-sided lighting** (`dot.abs()` in `rasterize_triangle_depth`, framebuffer.rs line 98):

```rust
let half_lambert = dot.abs() * 0.4 + 0.6;
```

The `.abs()` means back-facing triangles get the same lighting as front-facing triangles. Additionally, the rasterizer uses a **signed-area barycentric test** (framebuffer.rs lines 131-135):

```rust
let denom = (v1[1] - v2[1]) * (v0[0] - v2[0]) + (v2[0] - v1[0]) * (v0[1] - v2[1]);
if denom.abs() < 1e-12 { return; }  // degenerate check uses abs
```

The `denom` is the signed area of the triangle. For reversed winding, `denom` becomes negative, but the rasterizer still works because the barycentric coordinate computation uses `inv_denom = 1.0 / denom`, and the inside test is `u >= -1e-6 && v >= -1e-6 && w >= -1e-6`. When winding reverses, `u` and `v` become negative for points that were previously inside -- **BUT WAIT**: with negative `denom`, the sign of `u` and `v` also flips, so the test still works correctly. Let me verify:

With negative `denom`, `inv_denom` is negative. The barycentric weights are multiplied by `inv_denom`, and with the opposite winding, the numerator also flips sign. So `u`, `v`, `w` remain positive for interior points. The test passes correctly regardless of winding.

**Conclusion on winding:** The triangle rasterizer is winding-order-agnostic due to:
1. Two-sided lighting (`dot.abs()`)
2. Barycentric rasterization that handles both CW and CCW triangles

**Risk: LOW.** No changes needed to ribbon winding order. The fork also did NOT change ribbon winding order (confirmed by examining the diff -- the winding in `emit_strip` and `emit_quad` is unchanged).

### 3.4 `src/render/framebuffer.rs` -- Depth buffer and rasterization

The depth buffer uses Z values directly from `project()`. In our current code, `z2 = -x * sin_y + z1 * cos_y` (the rotated Z after Y rotation). This is not affected by negating the output X.

The rasterizer itself (as analyzed in 3.3) is winding-order-agnostic.

**Risk: NONE.** No changes needed.

### 3.5 `src/app.rs` -- Pan direction

Pan is implemented as:
```rust
pub fn pan(&mut self, dx: f64, dy: f64) {
    self.pan_x += dx * Self::PAN_STEP;
    self.pan_y += dy * Self::PAN_STEP;
}
```

In the projection: `x: -x3 * self.zoom + self.pan_x` (after fix).

When the user presses 'a' (mapped to `pan(-1.0, 0.0)`), `pan_x` decreases by `PAN_STEP`. In the projection, `pan_x` is added directly to the output X. So pressing 'a' moves the image left on screen. This is the correct/intuitive behavior: 'a' = pan view left = image moves left.

Before the fix: `x: x3 * self.zoom + self.pan_x` -- pressing 'a' also moved the image left. So pan behavior is the same before and after the fix because pan operates in screen space (added after the X negation).

**However**, consider the conceptual interpretation: pressing 'a' should move the *camera* left, which should make the *object* appear to shift right. Under the current mapping (`pan(-1.0, 0.0)` for 'a'), pan_x decreases, which shifts the projection left. This means pressing 'a' shifts the object left, as if moving the viewport right (or the object left). This is "drag-to-pan" semantics (the object follows the cursor direction), which is the standard for most viewers.

**After applying the fix, pan continues to work correctly and intuitively.** No negation of `pan_x` is needed.

**Risk: NONE.** No changes needed.

### 3.6 Key rotation direction (h/l for Y, j/k for X)

**Y-axis rotation (h/l):**
- `h` maps to `rotate_y(-1.0)` which decreases `rot_y`
- `l` maps to `rotate_y(1.0)` which increases `rot_y`

Before the fix, increasing `rot_y` rotates the object such that points on the +X side move toward the viewer (toward -Z). On screen (where +X is right), this appears as the right side moving toward the viewer = clockwise rotation when viewed from above. Pressing `l` (right) causes clockwise rotation from above.

After the fix, with negated screen X, the same `rot_y` increase still causes the same 3D rotation. But now +X world maps to -X screen (left). So the left side of the screen moves toward the viewer. Pressing `l` (right) now causes what appears to be counter-clockwise rotation from above, or equivalently, the object appears to rotate to the LEFT. This is reversed from the expected behavior (pressing right should rotate right).

The fork compensates for this by... actually, looking at the fork's main.rs diff, h/l mappings are unchanged: `h` = `rotate_y(-1.0)`, `l` = `rotate_y(1.0)`. But the fork uses a completely different rotation system (orientation matrices with `apply_local_rotation`), and critically, the auto-rotate direction is negated. Let me re-examine.

**Wait -- re-analysis:** Actually, let me trace through more carefully with our Euler system.

At identity, `project(1, 0, 0)` gives `x = 1` (before fix) or `x = -1` (after fix).

After pressing `l` once: `rot_y += 0.1`. Now for point (1, 0, 0):
- Ry: x2 = 1*cos(0.1) + 0*sin(0.1) = cos(0.1) ~= 0.995; z2 = -sin(0.1) ~= -0.0998
- After fix: projected x = -0.995 * zoom

So a point originally at x=-1 (screen) moves to x=-0.995 (screen) -- it moved slightly to the right. This means pressing `l` makes the object rotate so the left side comes toward the viewer and the right side goes away. From the user's perspective looking at the screen, this looks like the object rotates to the right (the right side recedes, the left side comes forward). That IS the correct/intuitive behavior for pressing `l` (right).

Wait, let me reconsider. If the right side goes away and the left side comes forward, in a horizontal rotation, this means the object is rotating counter-clockwise when viewed from above, which for a viewer looking at the front, means the object's left side comes toward them. This could be described as "rotating right" in the sense that the object turns to show its left side, which is what you'd see if you rotated the camera to the right. So actually this might be fine.

Let me try a simpler analysis. Before the fix, auto-rotate does `rot_y += speed * dt` (positive direction). This makes the molecule appear to spin. After the fix with negated X, the same auto-rotation would spin in the opposite visual direction. The fork negates the auto-rotate: `self.apply_local_rotation([0.0, 1.0, 0.0], -Self::AUTO_ROTATE_SPEED * dt)`.

**Conclusion on Y rotation:** After negating X, the Y-axis rotation direction appears reversed. We need to **negate `rot_y` increments** to compensate, or equivalently:
- Negate the auto-rotate direction: `self.rot_y -= Self::AUTO_ROTATE_SPEED * dt`
- Swap h/l rotation directions, OR negate the rotation step in `rotate_y`

The simplest approach: negate the sign in `rotate_y`:
```rust
pub fn rotate_y(&mut self, dir: f64) { self.rot_y -= dir * Self::ROT_STEP; }
//                                             ^^ was +=
```
This reverses Y rotation to compensate for the X flip.

**X-axis rotation (j/k):**
The fork swapped j/k directions:
- j/Down: changed from `rotate_x(1.0)` to `rotate_x(-1.0)`
- k/Up: changed from `rotate_x(-1.0)` to `rotate_x(1.0)`

However, X-axis rotation (tilting up/down) should NOT be affected by a screen-space X flip. Let me trace through:

At identity + small rot_x=0.1, for point (0, 1, 0):
- Rx: y1 = cos(0.1), z1 = sin(0.1)
- Ry (identity): x2=0, z2=sin(0.1)
- Rz (identity): x3=0, y3=cos(0.1)
- After fix: projected x=0, y=cos(0.1)*zoom

So pressing j (which does rot_x += 0.1) makes a point at +Y move slightly downward (cos decreases from 1). This looks like tilting the object forward (top tilts away). Pressing j (down arrow) tilting the top away is correct/intuitive.

After the X flip: x3 is negated to -x3, but x3=0 here, so no effect on Y-only points. For a general point, X rotation affects Y and Z, not X directly. The rotated X component (x3) IS affected by X rotation through the Z rotation coupling, but for small angles at identity Z rotation, x3 ~= x2 which is the Y-rotated x. So X rotation primarily moves things vertically on screen.

**I believe the fork swapped j/k for a different reason** -- possibly because their orientation matrix system composes rotations differently (pre-multiply vs post-multiply). In our Euler angle system, negating screen X should not require swapping j/k. However, this needs testing.

**Risk: MEDIUM.** Y rotation direction (h/l) and auto-rotate will likely need direction reversal.

### 3.7 Z-axis rotation (u/i)

Similar to Y rotation, Z rotation (`u`/`i`) may also need direction reversal after the X flip. When screen X is negated, a clockwise rotation in 3D maps to counter-clockwise on screen. This means `rotate_z` should also be negated.

**Risk: MEDIUM.** May need direction reversal -- test empirically.

## 4. Pitfall Analysis

### Pitfall 1: Triangle winding order (RISK: LOW)

**Assessment:** Not a problem. As analyzed in section 3.3, the rasterizer is winding-order-agnostic due to:
- Two-sided lighting with `dot.abs()`
- Barycentric rasterization that handles both CW and CCW winding

The fork did NOT change any winding order in ribbon.rs, confirming this is not an issue.

### Pitfall 2: Normal vectors for Lambert shading (RISK: LOW)

**Assessment:** Not a problem. `rotate_normal()` in hd.rs applies the same rotation to normals as `project()` applies to positions (minus the X negation and zoom/pan). Since the lighting uses `dot.abs()` (two-sided), even if the normal direction flips relative to the light, the shading result is identical.

The fork replaced `rotate_normal()` with `camera.rotate_vector()` because they changed the rotation system to orientation matrices. In our case, `rotate_normal()` stays as-is.

### Pitfall 3: Pan controls (RISK: NONE)

**Assessment:** Not a problem. Pan operates in screen space (added after the rotation and X negation), so `pan_x` changes have the same screen-space effect regardless of the X flip. Pressing 'a' still moves the image left.

### Pitfall 4: Rotation controls (RISK: MEDIUM)

**Assessment:** This IS a real issue that needs compensation.

After negating X in the projection:
- **Y rotation (h/l):** Visual direction reverses. Fix: negate `rot_y` increment direction.
- **Auto-rotate:** Visual direction reverses. Fix: negate auto-rotate increment.
- **X rotation (j/k):** Probably unaffected for our Euler system, but needs empirical testing.
- **Z rotation (u/i):** Likely reverses. Fix: negate `rot_z` increment direction.

**Safest approach:** Negate the sign in `rotate_y` and `rotate_z` methods, and negate the auto-rotate increment. Leave `rotate_x` unchanged initially, test, and adjust if needed.

Alternatively, a cleaner approach: instead of negating the projection output, negate the rotation angle for Y rotation only. But this is less principled -- the X flip IS the correct fix.

## 5. Verification Approach

### Visual verification methods

1. **Known chirality test:** Load a well-known protein with obvious chirality. PDB 1CRN (crambin) is a small protein (46 residues) where the backbone trace has a clear handedness. In the corrected view, alpha-helices should coil in the right-handed direction when viewed from N-terminus to C-terminus.

2. **Alpha-helix handedness:** All natural alpha-helices are right-handed. View a helix along its axis (from N to C terminus). It should wind clockwise going away from the viewer. This is the quickest visual test.

3. **Asymmetric structure test:** Use any PDB file that contains a structure with a clearly asymmetric feature visible in a specific orientation. Compare with the same structure viewed in PyMOL, ChimeraX, or the RCSB 3D viewer.

4. **Simple coordinate test:** Create a test PDB with atoms at known coordinates (e.g., three atoms forming an L-shape at +X, origin, +Y). Verify that the +X atom appears on the LEFT side of the screen (since in standard biochemistry viewing convention, +X is to the viewer's left when looking along -Z).

5. **Comparison with fork:** Run both our fixed version and the fork's version side by side on the same PDB file and confirm they produce visually identical (non-mirrored) output.

### Automated regression test

After the fix, add to camera tests:

```rust
#[test]
fn project_identity_x_is_negative() {
    // At identity rotation, a point at +X world should project to
    // -X screen (left side) when looking along -Z.
    let camera = Camera::default();
    let p = camera.project(1.0, 0.0, 0.0);
    assert!(p.x < 0.0, "positive world X should project to negative screen X");
}
```

## 6. Test Strategy

### Unit tests for projection correctness

```rust
#[test]
fn project_not_mirrored() {
    // At identity rotation, world +X should map to screen -X.
    // This ensures the projection matches the standard right-handed
    // viewing convention (camera looks along -Z, +X is to the left).
    let cam = Camera::default();
    let p = cam.project(1.0, 0.0, 0.0);
    assert!(p.x < 0.0, "world +X should project to screen left (negative X)");
    assert!((p.y).abs() < 1e-9, "world +X should have zero Y projection");
}

#[test]
fn project_y_unchanged() {
    // World +Y should still project to screen +Y.
    let cam = Camera::default();
    let p = cam.project(0.0, 1.0, 0.0);
    assert!(p.y > 0.0, "world +Y should project to screen top (positive Y)");
    assert!((p.x).abs() < 1e-9, "world +Y should have zero X projection");
}

#[test]
fn project_z_depth() {
    // World +Z should project to positive depth (farther from viewer).
    let cam = Camera::default();
    let p = cam.project(0.0, 0.0, 1.0);
    assert!(p.z > 0.0, "world +Z should be positive depth");
}
```

### Integration test for rotation compensation

```rust
#[test]
fn rotate_y_visual_direction() {
    // After a small positive Y rotation (pressing 'l'), the visible
    // rotation should be such that the right side of the object moves
    // away from the viewer (positive depth).
    let mut cam = Camera::default();
    cam.rotate_y(1.0);  // simulate pressing 'l'

    // A point on the right side of the screen (world -X after fix)
    // should gain positive depth after right rotation.
    let p = cam.project(-1.0, 0.0, 0.0);
    // After the fix + rotation compensation, this should behave intuitively.
    // Specific assertions depend on the final rotation sign convention.
}
```

### Visual regression test (manual)

1. Load `examples/1AOI.pdb` (or any available example)
2. Take a screenshot at identity rotation
3. Compare with reference image from PyMOL/ChimeraX at the same orientation
4. Verify helices, sheets, and overall fold chirality match

## 7. Summary of Required Changes

### Minimal fix (in priority order)

1. **camera.rs line 83**: Negate X in projection output
   ```
   x: -x3 * self.zoom + self.pan_x
   ```

2. **camera.rs line 45**: Negate Y rotation direction to compensate
   ```
   pub fn rotate_y(&mut self, dir: f64) { self.rot_y -= dir * Self::ROT_STEP; }
   ```

3. **camera.rs line 60**: Negate auto-rotate direction
   ```
   self.rot_y -= Self::AUTO_ROTATE_SPEED * dt;
   ```

4. **camera.rs line 46**: Likely negate Z rotation direction (test to confirm)
   ```
   pub fn rotate_z(&mut self, dir: f64) { self.rot_z -= dir * Self::ROT_STEP; }
   ```

5. **No changes needed:**
   - `rotate_normal()` in hd.rs (stays the same)
   - `to_pixel()` in hd.rs (stays the same)
   - Triangle winding in ribbon.rs (stays the same)
   - Rasterizer in framebuffer.rs (stays the same)
   - Pan in camera.rs / main.rs (stays the same)
   - Braille renderer (stays the same)
   - X rotation direction (probably stays the same -- test to confirm)

### Files to modify

| File | Change | Risk |
|------|--------|------|
| `src/render/camera.rs` | Negate X in `project()` output (line 83) | Core fix |
| `src/render/camera.rs` | Negate Y rotation direction (line 45) | Compensating |
| `src/render/camera.rs` | Negate auto-rotate direction (line 60) | Compensating |
| `src/render/camera.rs` | Possibly negate Z rotation direction (line 46) | Test first |

### Files confirmed safe (no changes needed)

| File | Reason |
|------|--------|
| `src/render/hd.rs` | `rotate_normal()` and `to_pixel()` unaffected |
| `src/render/braille.rs` | Uses `project()` output directly, no X assumptions |
| `src/render/ribbon.rs` | World-space mesh, winding-order-agnostic rasterizer |
| `src/render/framebuffer.rs` | Two-sided lighting, winding-agnostic barycentric rasterization |
| `src/app.rs` | Pan is in screen space, unaffected by X flip |
| `src/main.rs` | Key bindings unchanged (rotation sign changes are in camera.rs) |

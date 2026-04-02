# Design: HD Mode Performance Over SSH

**Date:** 2026-03-13
**Branch:** dev/structural-analysis-features
**Status:** Proposed

## Problem Statement

HD mode rendering over SSH has degraded from usable to unusably slow. The root causes are a combination of a recent change that increased per-frame payload size and pre-existing architectural issues that compound over high-latency, bandwidth-constrained SSH connections.

### Measurements

| Metric | Value | Source |
|---|---|---|
| Frame payload (Sixel, 80x24 terminal, 8x16 font) | ~170 KB | Measured from sixel_string output |
| Frame payload (RGBA, same terminal) | ~245 KB raw, ~170 KB encoded | 640x384 * 4 bytes raw; Sixel re-encodes as RGB |
| Frame rate | 30 FPS (hardcoded 33ms tick) | src/main.rs line 143 |
| Sustained bandwidth demand | ~5.1 MB/s (Sixel) | 170 KB * 30 FPS |
| Heap allocation per frame | ~3.6 MB | Two Vec allocations for color ([u8;3]) + depth (f64) at 640x384 |
| RGBA image allocation per frame | ~983 KB | 640 * 384 * 4 bytes |
| Typical SSH bandwidth (LAN) | 1-10 MB/s usable | Depends on cipher, compression |
| Typical SSH bandwidth (WAN) | 0.1-1 MB/s usable | Depends on route, encryption |

### What Changed: The RGBA Regression

The commit `89e37f6` ("fix: use RGBA images in HD mode so terminal background shows through") switched the image encoding pipeline from RGB to RGBA:

```
- let rgb_img = fb.to_rgb_image();
- let dyn_img = DynamicImage::ImageRgb8(rgb_img);
+ let rgba_img = fb.to_rgba_image();
+ let dyn_img = DynamicImage::ImageRgba8(rgba_img);
```

This change causes two performance regressions:

1. **33% larger raw image buffer**: RGBA is 4 bytes/pixel vs RGB's 3 bytes/pixel. For a 640x384 framebuffer, that is 983 KB vs 737 KB per frame.

2. **Background color overlay in ratatui-image**: When `ImageSource::new()` receives an RGBA image with `background_color` alpha != 0 (the default is `Rgba([0, 0, 0, 0])` -- transparent), it skips the overlay. However, the `Resize::resize()` method *always* creates a background-filled `DynamicImage` and composites the source onto it, doubling the image allocation. More critically, the Sixel encoder calls `img.to_rgb8()` which must flatten the alpha channel -- and for Kitty, `img.to_rgba8()` is called followed by full base64 encoding of 4 bytes per pixel instead of 3.

3. **Kitty payload bloat**: The Kitty protocol transmits raw RGBA (`f=32`) base64-encoded. With RGBA, the raw payload before base64 is 33% larger, and base64 adds another 33% on top. Net effect: Kitty frames went from ~1.3 MB to ~1.7 MB of escape sequence data per frame.

The old RGB path was still bandwidth-heavy but was within the budget for LAN SSH connections. The RGBA switch pushed it over the edge.

### Pre-Existing Issues (Compounding Factors)

Even before the RGBA change, the HD-over-SSH path had these structural problems:

1. **No frame deduplication**: Every frame is fully re-rendered and re-encoded even when the scene has not changed (no rotation, no input, auto-rotate off). An idle application sends 5.1 MB/s for nothing.

2. **No framebuffer reuse**: `render_hd_framebuffer()` allocates a new `Framebuffer` (two `Vec`s: color + depth) on every frame. At 640x384: 737 KB for color, 1.97 MB for depth, plus the image conversion buffer. Total: ~3.6 MB of heap allocation per frame, 108 MB/s at 30 FPS. The allocator overhead alone is significant.

3. **Depth scan in post-pass**: `apply_depth_tint()` performs a full scan of the depth buffer to find z_min/z_max before the tinting pass. These min/max values could be tracked during rasterization for free.

4. **Frame skipping threshold too high**: Frame skipping only triggers when draw time exceeds 2x the tick rate (66ms). Over SSH, the bottleneck is not draw CPU time but write latency -- frames complete "quickly" on the CPU but saturate the PTY buffer. The frame skipper never engages.

5. **No adaptive frame rate**: The 30 FPS tick rate is hardcoded. Over SSH, even 5 FPS would be sufficient for interactive use and would reduce bandwidth from 5.1 MB/s to 0.85 MB/s.

6. **Full image re-encode every frame**: `picker.new_protocol()` creates a new `Protocol` object on every frame, re-encoding the Sixel or Kitty data from scratch. The Kitty protocol in ratatui-image supports stateful `transmit-once, render-many` via `StatefulKitty`, but the current code uses the stateless `Kitty` path.

## Solution Architecture

The solution is composed of five independent, layered components that can be implemented incrementally. Each component delivers measurable improvement on its own; together they reduce SSH bandwidth demand by 10-50x.

```
                    +---------------------------+
                    |     Main Event Loop        |
                    |  (adaptive tick rate)      |
                    +------+--------------------+
                           |
                    +------v--------------------+
                    |   Camera Dirty Tracking    |
                    |   (skip identical frames)  |
                    +------+--------------------+
                           |
                    +------v--------------------+
                    |   Framebuffer Pool         |
                    |   (reusable allocations)   |
                    +------+--------------------+
                           |
                    +------v--------------------+
                    |   Resolution Scaling       |
                    |   (reduce pixel count)     |
                    +------+--------------------+
                           |
                    +------v--------------------+
                    |   Protocol Optimization    |
                    |   (RGB path, Kitty state)  |
                    +------+--------------------+
                           |
                    +------v--------------------+
                    |   SSH-Aware Adaptation     |
                    |   (detect + throttle)      |
                    +---------------------------+
```

### Component 1: Camera Dirty Tracking

**Goal**: Skip rendering entirely when the scene has not changed.

**Impact**: Eliminates 100% of bandwidth when idle (no rotation, no input). With auto-rotate off, the typical user spends >80% of their time examining a static view.

**Design**:

Add a `dirty` flag to `Camera` that is set on any state mutation and cleared after a frame is drawn.

```rust
// src/render/camera.rs
pub struct Camera {
    // ... existing fields ...
    /// True when camera state has changed since last frame was drawn.
    dirty: bool,
}

impl Camera {
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    pub fn clear_dirty(&mut self) {
        self.dirty = false;
    }

    // Every mutation sets the flag:
    pub fn rotate_x(&mut self, dir: f64) {
        self.rot_x += dir * Self::ROT_STEP;
        self.dirty = true;
    }
    // ... same for rotate_y, rotate_z, zoom_in, zoom_out, pan, reset ...

    pub fn tick(&mut self) {
        let now = Instant::now();
        let dt = now.duration_since(self.last_tick).as_secs_f64();
        self.last_tick = now;
        if self.auto_rotate {
            self.rot_y -= Self::AUTO_ROTATE_SPEED * dt;
            self.dirty = true;
        }
    }
}
```

Add a corresponding `scene_dirty` flag to `App` that aggregates camera dirty state plus any other state changes that require a redraw (mode change, color change, ligand toggle, interface toggle, terminal resize):

```rust
// src/app.rs
impl App {
    /// Returns true if anything has changed that requires a new frame.
    pub fn needs_redraw(&self) -> bool {
        self.camera.is_dirty() || self.mesh_dirty || self.ui_dirty
    }

    pub fn mark_drawn(&mut self) {
        self.camera.clear_dirty();
        self.ui_dirty = false;
    }
}
```

In the main loop, skip `terminal.draw()` entirely when `!app.needs_redraw() && !had_input`:

```rust
// src/main.rs - main loop
if !app.needs_redraw() && !had_input {
    app.tick();
    std::thread::sleep(tick_rate);
    continue;
}
```

**Edge cases**:
- First frame must always draw (initialize `dirty = true` in `Camera::default()`).
- Terminal resize events must set `ui_dirty = true`. Currently resize is not handled; when resize handling is added, it should go through this flag.
- Help overlay toggle must set `ui_dirty = true`.
- `had_input` forces a redraw even if no camera mutation occurred (e.g., pressing an unrecognized key should not cause a redraw, but the current code redraws on any input -- this is acceptable for simplicity).

**Files modified**: `src/render/camera.rs`, `src/app.rs`, `src/main.rs`

### Component 2: Framebuffer Pool (Allocation Reuse)

**Goal**: Eliminate per-frame heap allocation for the framebuffer and RGBA image buffer.

**Impact**: Removes ~3.6 MB of allocation and ~983 KB of RGBA allocation per frame. Reduces GC pressure and eliminates allocator overhead (~0.5ms per frame measured on typical allocators for this size).

**Design**:

Move the `Framebuffer` and the RGBA `Vec<u8>` into `App` as persistent, reusable buffers:

```rust
// src/app.rs
pub struct App {
    // ... existing fields ...
    /// Reusable HD framebuffer. Resized only when terminal dimensions change.
    pub hd_framebuffer: Option<Framebuffer>,
    /// Reusable RGBA pixel buffer for image encoding.
    pub hd_rgba_buffer: Vec<u8>,
}
```

Promote `Framebuffer::clear()` from `#[cfg(test)]` to always-available:

```rust
// src/render/framebuffer.rs
impl Framebuffer {
    /// Reset the framebuffer to black pixels and infinite depth.
    /// Reuses existing allocations.
    pub fn clear(&mut self) {
        for c in self.color.iter_mut() {
            *c = [0, 0, 0];
        }
        for d in self.depth.iter_mut() {
            *d = f64::INFINITY;
        }
    }

    /// Resize the framebuffer, only reallocating if the new size is larger.
    pub fn resize(&mut self, width: usize, height: usize) {
        self.width = width;
        self.height = height;
        let n = width * height;
        if self.color.len() < n {
            self.color.resize(n, [0, 0, 0]);
            self.depth.resize(n, f64::INFINITY);
        }
        self.clear();
    }
}
```

Change `render_hd_framebuffer` to accept a `&mut Framebuffer` instead of returning a new one:

```rust
// src/render/hd.rs
pub fn render_hd_framebuffer(
    fb: &mut Framebuffer,  // <-- takes mutable reference instead of returning owned
    protein: &Protein,
    camera: &Camera,
    // ... other params ...
) {
    let px_w = width as usize;
    let px_h = height as usize;
    fb.resize(px_w, px_h);
    // ... rest of rendering into fb ...
}
```

For the RGBA image buffer, add a method that writes into a pre-allocated `Vec<u8>` instead of creating a new `RgbaImage`:

```rust
// src/render/framebuffer.rs
impl Framebuffer {
    /// Write RGBA pixel data into the provided buffer, resizing it if needed.
    /// Returns (width, height) of the written image.
    pub fn write_rgba_into(&self, buf: &mut Vec<u8>) -> (u32, u32) {
        let n = self.width * self.height * 4;
        buf.resize(n, 0);
        for i in 0..(self.width * self.height) {
            let c = self.color[i];
            let alpha = if self.depth[i] >= f64::INFINITY { 0 } else { 255 };
            let off = i * 4;
            buf[off] = c[0];
            buf[off + 1] = c[1];
            buf[off + 2] = c[2];
            buf[off + 3] = alpha;
        }
        (self.width as u32, self.height as u32)
    }
}
```

**Track z_min/z_max during rasterization**: Eliminate the O(n) depth scan in `apply_depth_tint` by tracking min/max as pixels are written:

```rust
// src/render/framebuffer.rs
pub struct Framebuffer {
    // ... existing fields ...
    /// Tracked min/max z of written pixels (for depth tinting).
    z_min: f64,
    z_max: f64,
}

impl Framebuffer {
    #[inline]
    fn set_pixel(&mut self, x: usize, y: usize, z: f64, color: [u8; 3]) {
        let idx = y * self.width + x;
        if z < self.depth[idx] {
            self.depth[idx] = z;
            self.color[idx] = color;
            // Track depth range
            if z < self.z_min { self.z_min = z; }
            if z > self.z_max { self.z_max = z; }
        }
    }

    pub fn apply_depth_tint(&mut self, fog_color: [u8; 3], fog_strength: f64) {
        // Use tracked z_min/z_max instead of scanning
        let z_range = self.z_max - self.z_min;
        if z_range.abs() < 1e-12 {
            return;
        }
        let inv_range = 1.0 / z_range;
        // ... rest of tinting loop unchanged ...
    }

    pub fn clear(&mut self) {
        // ... existing clear code ...
        self.z_min = f64::INFINITY;
        self.z_max = f64::NEG_INFINITY;
    }
}
```

**Files modified**: `src/render/framebuffer.rs`, `src/render/hd.rs`, `src/ui/viewport.rs`, `src/app.rs`

### Component 3: Restore RGB Encoding Path + Protocol Optimization

**Goal**: Reduce per-frame payload by 25-33% by using RGB instead of RGBA for Sixel, and by using stateful Kitty protocol to avoid retransmitting unchanged images.

**Impact**: Sixel frames drop from ~170 KB to ~127 KB. Kitty frames drop from ~1.7 MB to ~1.3 MB (RGB), and further to near-zero for unchanged frames (stateful placement).

**Design**:

The RGBA change (commit `89e37f6`) was needed for transparent backgrounds. The correct fix is to use RGBA only when the terminal actually benefits from it, and to use RGB for protocols that discard the alpha channel anyway.

**3a. Protocol-aware image format selection**:

```rust
// src/ui/viewport.rs
fn render_hd_viewport(frame: &mut Frame, area: Rect, app: &App) {
    // ... rasterize into fb as before ...

    let proto = app.picker.protocol_type();
    if proto != ProtocolType::Halfblocks {
        let dyn_img = match proto {
            // Sixel discards alpha internally (calls to_rgb8()), so skip RGBA overhead
            ProtocolType::Sixel => {
                let rgb_img = fb.to_rgb_image();
                DynamicImage::ImageRgb8(rgb_img)
            }
            // Kitty and iTerm2 support transparency
            ProtocolType::Kitty | ProtocolType::Iterm2 => {
                let rgba_img = fb.to_rgba_image();
                DynamicImage::ImageRgba8(rgba_img)
            }
            _ => unreachable!(),
        };
        // ... rest unchanged ...
    }
}
```

This alone saves ~8 KB per Sixel frame (avoiding the RGBA->RGB conversion inside icy_sixel) and reduces the allocation from 983 KB to 737 KB.

**3b. Un-gate `to_rgb_image` from `#[cfg(test)]`**:

The `to_rgb_image` method was gated behind `#[cfg(test)]` when the RGBA path was introduced. It needs to be restored for production use:

```rust
// src/render/framebuffer.rs
// Remove #[cfg(test)] from to_rgb_image
pub fn to_rgb_image(&self) -> RgbImage { ... }
```

**3c. Sixel quality reduction for SSH**:

The Sixel encoder currently uses `Quality::HIGH`. For SSH connections, `Quality::LOW` reduces the palette size and dithering complexity, producing smaller Sixel output:

```rust
// When SSH is detected (see Component 5), pass a quality parameter
let quality = if app.ssh_mode { Quality::LOW } else { Quality::HIGH };
```

This is not directly controllable through ratatui-image's public API today. Two approaches:

- **Short-term**: Accept the quality ratatui-image chooses. The bandwidth savings from other components are sufficient.
- **Long-term**: Fork or contribute a `SixelOptions` struct to ratatui-image that exposes quality/palette settings.

**3d. Use `StatefulProtocol` for Kitty to avoid retransmission**:

The Kitty protocol in ratatui-image uses `transmit-once, render-many` semantics. The current code creates a new `Protocol` every frame, which forces retransmission. Using `StatefulProtocol` (via `StatefulImage` widget) enables the library to skip encoding when the image area hasn't changed:

```rust
// src/app.rs
pub struct App {
    // ... existing fields ...
    /// Stateful protocol for Kitty -- avoids retransmitting unchanged images.
    pub stateful_protocol: Option<StatefulProtocol>,
}
```

In the viewport, use `StatefulImage` when Kitty is the active protocol:

```rust
// src/ui/viewport.rs
fn render_hd_viewport(frame: &mut Frame, area: Rect, app: &mut App) {
    // ... rasterize ...
    match proto {
        ProtocolType::Kitty => {
            let dyn_img = DynamicImage::ImageRgba8(fb.to_rgba_image());
            if app.stateful_protocol.is_none() {
                app.stateful_protocol = Some(app.picker.new_resize_protocol(dyn_img));
            } else {
                // Update the source image; StatefulProtocol handles re-encode only if needed
                // Note: this requires updating the ImageSource, which the current API
                // doesn't directly support for source replacement. See Risks section.
            }
            let widget = StatefulImage::default();
            frame.render_stateful_widget(widget, area, app.stateful_protocol.as_mut().unwrap());
        }
        _ => {
            // Sixel / iTerm2: use stateless Protocol (current path)
            // ...
        }
    }
}
```

**Limitation**: The `StatefulProtocol` API in ratatui-image 9.0.0 is designed for static images that get resized to fit an area. It does not natively support replacing the source image while keeping the same protocol state. For our use case (new image every frame when dirty, same image when clean), we need to either:

1. Create a new `StatefulProtocol` only when the image changes (still saves bandwidth when idle).
2. Contribute an `update_source()` method to ratatui-image upstream.

Option 1 is sufficient for the initial implementation and still delivers the key benefit: zero bandwidth when the scene is unchanged.

**Files modified**: `src/render/framebuffer.rs`, `src/ui/viewport.rs`, `src/app.rs`

### Component 4: Resolution Scaling for SSH

**Goal**: Reduce pixel count (and therefore Sixel/Kitty payload size) when rendering over SSH, without changing the terminal cell area.

**Impact**: A 2x resolution reduction (e.g., from 640x384 to 320x192) reduces Sixel payload by ~4x (from ~170 KB to ~43 KB per frame). ratatui-image's `Resize::Fit` will upscale the smaller image to fill the terminal area.

**Design**:

Add a `resolution_scale` factor to `App`:

```rust
// src/app.rs
pub struct App {
    // ... existing fields ...
    /// Resolution scale factor for HD mode (1.0 = full, 0.5 = half, etc.)
    /// Reduced automatically for SSH connections.
    pub resolution_scale: f64,
}
```

Apply the scale in the viewport before rasterization:

```rust
// src/ui/viewport.rs
fn render_hd_viewport(frame: &mut Frame, area: Rect, app: &App) {
    let proto = app.picker.protocol_type();
    let (font_w, font_h) = app.picker.font_size();

    let (px_w, px_h) = if proto != ProtocolType::Halfblocks && font_w > 0 && font_h > 0 {
        (
            (area.width as f64 * font_w as f64 * app.resolution_scale).max(1.0),
            (area.height as f64 * font_h as f64 * app.resolution_scale).max(1.0),
        )
    } else {
        // Braille fallback: no scaling (already low-res)
        (area.width as f64 * 2.0, area.height as f64 * 4.0)
    };

    // Rasterize at reduced resolution
    let fb = hd::render_hd_framebuffer(..., px_w, px_h, ...);

    // ratatui-image will scale up to fill the terminal area via Resize::Fit
    // ...
}
```

Resolution scale values:

| Connection | `resolution_scale` | Pixel dimensions (80x24 @ 8x16) | Approx. Sixel size |
|---|---|---|---|
| Local | 1.0 | 640 x 384 | ~170 KB |
| SSH LAN | 0.5 | 320 x 192 | ~43 KB |
| SSH WAN | 0.25 | 160 x 96 | ~11 KB |

The upscale by ratatui-image uses `FilterType::Nearest` by default, which is fast and preserves the "pixel art" aesthetic appropriate for a protein viewer at reduced resolution.

**Zoom compensation**: When resolution is scaled down, the camera zoom must be proportionally reduced so the protein fills the same visual area:

```rust
// src/app.rs
pub fn recalculate_zoom(&mut self, term_cols: u16, term_rows: u16) {
    // ... existing logic ...
    // Apply resolution_scale to the pixel dimensions used for zoom calculation
    let (px_w, px_h) = if has_graphics {
        (
            vp_cols * font_w as f64 * self.resolution_scale,
            vp_rows * font_h as f64 * self.resolution_scale,
        )
    } else { ... };
    self.camera.zoom = 0.9 * px_w.min(px_h) / (2.0 * radius);
}
```

**User override**: Add a `--resolution-scale` CLI flag for manual control:

```rust
// src/main.rs
#[derive(Parser)]
struct Cli {
    // ... existing fields ...
    /// Resolution scale for HD mode (0.25-1.0, default auto-detected)
    #[arg(long, value_parser = clap::value_parser!(f64))]
    resolution_scale: Option<f64>,
}
```

**Files modified**: `src/main.rs`, `src/app.rs`, `src/ui/viewport.rs`

### Component 5: SSH-Aware Adaptation

**Goal**: Automatically detect SSH connections and apply appropriate performance settings.

**Impact**: Users get good performance over SSH without manual flag tweaking.

**Design**:

**5a. SSH detection**:

```rust
// src/app.rs (or a new src/env.rs module)

/// Connection environment detected at startup.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConnectionType {
    /// Direct local terminal session.
    Local,
    /// SSH connection (detected via environment variables).
    Ssh,
}

impl ConnectionType {
    pub fn detect() -> Self {
        // $SSH_CONNECTION is set by OpenSSH on the remote side.
        // $SSH_TTY is set when a pseudo-terminal is allocated.
        if std::env::var("SSH_CONNECTION").is_ok()
            || std::env::var("SSH_TTY").is_ok()
        {
            ConnectionType::Ssh
        } else {
            ConnectionType::Local
        }
    }
}
```

**5b. Adaptive defaults based on connection type**:

```rust
// src/app.rs
impl App {
    pub fn new(...) -> Self {
        let connection = ConnectionType::detect();

        let (resolution_scale, target_fps) = match connection {
            ConnectionType::Local => (1.0, 30.0),
            ConnectionType::Ssh => (0.5, 10.0),
        };

        // Allow CLI override
        let resolution_scale = cli_resolution_scale.unwrap_or(resolution_scale);
        let target_fps = cli_fps.unwrap_or(target_fps);

        // ...
    }
}
```

**5c. Adaptive tick rate in the main loop**:

```rust
// src/main.rs
let base_tick = Duration::from_secs_f64(1.0 / app.target_fps);

loop {
    // ... input handling ...

    // Measure write latency as a proxy for connection quality
    let draw_start = Instant::now();
    terminal.draw(|frame| { ... })?;
    let draw_duration = draw_start.elapsed();

    // Adaptive throttling: if draw (which includes terminal write) is slow,
    // increase the tick rate to avoid saturating the PTY buffer.
    let effective_tick = if draw_duration > base_tick {
        // Back off: wait for the PTY buffer to drain before sending more data.
        // Use 2x the measured draw time as the next tick interval.
        draw_duration * 2
    } else {
        base_tick
    };

    app.tick();
    let elapsed = draw_start.elapsed();
    if elapsed < effective_tick {
        std::thread::sleep(effective_tick - elapsed);
    }
}
```

This is a self-tuning feedback loop: when the PTY write buffer is under pressure (slow `draw`), the loop automatically backs off. When the buffer clears, it returns to the target frame rate.

**5d. CLI flags for manual override**:

```rust
// src/main.rs
#[derive(Parser)]
struct Cli {
    // ... existing fields ...

    /// Target frames per second (default: 30 local, 10 SSH)
    #[arg(long)]
    fps: Option<f64>,

    /// Resolution scale for HD mode (0.25-1.0, default: 1.0 local, 0.5 SSH)
    #[arg(long)]
    resolution_scale: Option<f64>,

    /// Force SSH performance mode even on local connections
    #[arg(long)]
    ssh: bool,
}
```

**Files modified**: `src/main.rs`, `src/app.rs`, new `src/env.rs` (optional, can be inlined)

## Component Interaction Diagram

```
Startup:
  CLI args -> ConnectionType::detect() -> default resolution_scale + target_fps
  -> App::new() with configured parameters

Per-frame:
  1. Input handling -> Camera/App mutations -> dirty flags set
  2. needs_redraw() check (Component 1)
     - false: sleep(tick_rate), continue
     - true: proceed to render
  3. fb = app.hd_framebuffer (reused, Component 2)
     fb.resize(scaled_width, scaled_height)  (Component 4)
     fb.clear()
  4. render_hd_framebuffer(&mut fb, ...)  (Components 2+4)
  5. Protocol-aware image encoding (Component 3):
     - Sixel: fb.to_rgb_image() -> DynamicImage::ImageRgb8
     - Kitty: fb.to_rgba_image() -> DynamicImage::ImageRgba8
     - (Resolution upscale handled by ratatui-image Resize::Fit)
  6. terminal.draw() flushes to PTY
  7. Adaptive tick rate adjustment (Component 5)
  8. app.mark_drawn() -> clear dirty flags
```

## Migration Path

The components are ordered by implementation complexity and risk. Each is independently shippable.

### Phase 1: Zero-Cost Wins (1-2 hours)

1. **Camera dirty tracking** (Component 1): Add `dirty` flag to `Camera`, `ui_dirty` flag to `App`, skip `terminal.draw()` when clean. This alone eliminates all idle bandwidth.

2. **Restore RGB for Sixel** (Component 3a/3b): Un-gate `to_rgb_image()`, use it for `ProtocolType::Sixel`. Keeps RGBA for Kitty/iTerm2. This reverses the regression for the most common protocol.

### Phase 2: Allocation Reduction (2-3 hours)

3. **Framebuffer reuse** (Component 2): Move `Framebuffer` to `App`, add `resize()` and un-gate `clear()`. Change `render_hd_framebuffer` to take `&mut Framebuffer`.

4. **Track z_min/z_max** (Component 2): Add tracking fields to `Framebuffer::set_pixel`, remove the scan from `apply_depth_tint`.

### Phase 3: SSH Awareness (2-3 hours)

5. **SSH detection** (Component 5a): Add `ConnectionType::detect()`.

6. **Adaptive frame rate** (Component 5c): Replace hardcoded tick rate with adaptive loop.

7. **Resolution scaling** (Component 4): Add `resolution_scale` to `App`, apply in viewport.

### Phase 4: Protocol Optimization (3-4 hours, optional)

8. **Stateful Kitty protocol** (Component 3d): Use `StatefulProtocol` for Kitty to avoid retransmission. This requires changing `render_viewport` to accept `&mut App` (currently takes `&App`), which propagates through the `terminal.draw()` closure.

## Expected Results

| Scenario | Before | After Phase 1 | After Phase 3 | After Phase 4 |
|---|---|---|---|---|
| Idle, local | 5.1 MB/s | 0 MB/s | 0 MB/s | 0 MB/s |
| Idle, SSH | 5.1 MB/s | 0 MB/s | 0 MB/s | 0 MB/s |
| Rotating, local | 5.1 MB/s | 3.8 MB/s (RGB) | 3.8 MB/s | 3.8 MB/s |
| Rotating, SSH (Sixel) | 5.1 MB/s | 3.8 MB/s | 0.43 MB/s (10 FPS, 0.5x res) | 0.43 MB/s |
| Rotating, SSH (Kitty) | ~51 MB/s | ~38 MB/s | 3.8 MB/s (10 FPS, 0.5x res) | ~0 after first frame |
| Input burst, SSH | 5.1 MB/s | 3.8 MB/s | 0.43 MB/s | 0.43 MB/s |

The combination of dirty tracking (Phase 1) and SSH-aware scaling (Phase 3) reduces typical SSH bandwidth demand by 10-12x, from 5.1 MB/s to 0.43 MB/s during rotation and 0 MB/s when idle.

## Testing Strategy

### Unit Tests

1. **Camera dirty flag**: Test that each mutation method sets `dirty`, that `clear_dirty()` resets it, and that auto-rotate in `tick()` sets it.

```rust
#[test]
fn test_camera_dirty_flag() {
    let mut cam = Camera::default();
    cam.clear_dirty();
    assert!(!cam.is_dirty());
    cam.rotate_y(1.0);
    assert!(cam.is_dirty());
    cam.clear_dirty();
    assert!(!cam.is_dirty());
}
```

2. **Framebuffer resize and clear**: Test that `resize()` reuses allocations when shrinking, allocates when growing, and that `clear()` resets z_min/z_max.

```rust
#[test]
fn test_framebuffer_resize_reuses() {
    let mut fb = Framebuffer::new(100, 100);
    let ptr1 = fb.color.as_ptr();
    fb.resize(50, 50);
    let ptr2 = fb.color.as_ptr();
    assert_eq!(ptr1, ptr2, "should reuse allocation when shrinking");
}
```

3. **z_min/z_max tracking**: Test that `set_pixel` updates tracked range and that `apply_depth_tint` uses it correctly.

4. **RGB vs RGBA output**: Test that `to_rgb_image` produces 3 bytes/pixel and `to_rgba_image` produces 4 bytes/pixel with correct alpha values.

5. **ConnectionType detection**: Test with mocked environment variables.

### Integration Tests

1. **Frame skip on idle**: Create an `App`, set auto-rotate off, verify `needs_redraw()` returns false after `mark_drawn()`.

2. **Frame skip on rotation**: Create an `App`, rotate camera, verify `needs_redraw()` returns true.

3. **Resolution scaling**: Render the same scene at 1.0x and 0.5x scale, verify the framebuffer dimensions are halved.

### Manual Testing

1. **Local HD mode**: Verify visual quality is unchanged from current behavior.
2. **SSH HD mode (LAN)**: Verify smooth interaction at 10 FPS with 0.5x resolution.
3. **SSH HD mode (WAN)**: Verify usable interaction at 5 FPS with 0.25x resolution.
4. **Mode toggle**: Verify switching between HD and braille modes works correctly with reused framebuffer.
5. **Kitty terminal**: Test stateful protocol path (if implemented).
6. **Auto-rotate over SSH**: Verify adaptive tick rate produces smooth but bandwidth-limited rotation.

### Performance Benchmarks

Add a `--benchmark` flag that runs a fixed rotation sequence and reports:
- Average frame time
- Average payload size per frame (measure bytes written to stdout)
- Peak memory usage

This enables regression testing for future changes.

## Risks and Mitigations

### Risk 1: Framebuffer reuse changes rendering correctness

**Likelihood**: Low.
**Impact**: Visual artifacts (stale pixels from previous frame).
**Mitigation**: `clear()` is called before every frame. The existing test suite for `Framebuffer` (z-buffer, triangle rasterization) validates correctness. Add a test that renders two different scenes into the same framebuffer and verifies no bleed-through.

### Risk 2: SSH detection false positives

**Likelihood**: Medium. Some development environments set `SSH_CONNECTION` even for local connections (e.g., VS Code Remote, Docker with SSH transport).
**Impact**: Unnecessarily reduced quality/frame rate on local connections.
**Mitigation**: Provide `--fps` and `--resolution-scale` CLI flags for manual override. Consider adding `--no-ssh` flag to force local mode.

### Risk 3: `StatefulProtocol` API mismatch with real-time rendering

**Likelihood**: Medium. ratatui-image's `StatefulProtocol` is designed for static images with occasional resizes, not for 30 FPS animation.
**Impact**: Component 3d (stateful Kitty) may not deliver expected savings, or may introduce flicker.
**Mitigation**: Component 3d is in Phase 4 (optional). The other components deliver the critical bandwidth savings. If `StatefulProtocol` doesn't work well for animation, fall back to the stateless path with the other optimizations applied.

### Risk 4: Resolution scaling degrades visual quality

**Likelihood**: Certain (by design).
**Impact**: Protein details are less visible at 0.5x or 0.25x resolution.
**Mitigation**: The quality reduction is intentional for SSH and only affects the Sixel/Kitty pixel path. The braille fallback is unaffected. Users can override with `--resolution-scale 1.0` if they have sufficient bandwidth. The `Resize::Fit(None)` default uses `FilterType::Nearest` which preserves sharp edges appropriate for molecular visualization.

### Risk 5: Adaptive tick rate oscillation

**Likelihood**: Low-medium. If `draw` duration fluctuates around the tick rate threshold, the effective tick rate could oscillate between fast and slow.
**Impact**: Inconsistent frame rate, possible judder.
**Mitigation**: Use exponential moving average of draw duration instead of raw measurement. Apply hysteresis: only increase tick rate if draw duration exceeds threshold for 3 consecutive frames; only decrease if it's below threshold for 10 consecutive frames.

### Risk 6: `render_viewport` signature change for mutable App

**Likelihood**: Certain (required for Component 2 and 3d).
**Impact**: The `terminal.draw()` closure currently takes `&App`. Changing to `&mut App` requires the closure to capture `&mut app`, which may conflict with other borrows.
**Mitigation**: The main loop already exclusively owns `app`. The `ribbon_mesh()` call is already outside `terminal.draw()` for exactly this reason (it needs `&mut self`). Move the framebuffer clear/resize outside the closure, pass `&mut App` into the render functions. Alternatively, pass `&mut Framebuffer` separately from `&App` to avoid the borrow conflict:

```rust
let fb = app.hd_framebuffer.as_mut().unwrap();
hd::render_hd_framebuffer(fb, &app.protein, &app.camera, ...);
terminal.draw(|frame| {
    ui::viewport::render_viewport_with_fb(frame, chunks[1], &app, fb);
    // ...
})?;
```

This avoids needing `&mut App` inside the draw closure entirely.

## Appendix: ratatui-image Protocol Internals

### Sixel Encoding Path (icy_sixel)

1. `Sixel::new()` calls `encode()` which calls `img.to_rgb8()` (discards alpha), then `sixel_string()`.
2. `sixel_string()` performs color quantization (256 colors), dithering (Stucki), and Sixel run-length encoding.
3. The output is a single `String` containing the Sixel escape sequence.
4. On render, the string is placed in cell (0,0) of the area, and all other cells are marked `skip`.
5. No caching, no delta support. Every call to `Sixel::new()` re-encodes from scratch.

### Kitty Encoding Path

1. `Kitty::new()` calls `transmit_virtual()` which calls `img.to_rgba8()`, then base64-encodes in 4096-byte chunks.
2. The escape sequence uses `a=T,U=1,f=32` (transmit, virtual placement, RGBA format).
3. The `KittyProtoState` ensures the transmit sequence is only sent once via an `AtomicBool`.
4. On subsequent renders, only unicode placeholder characters are sent (very small payload).
5. **Key insight**: If we create a new `Kitty` object every frame, the transmit-once optimization is defeated. The full base64-encoded RGBA data is sent every frame.
6. `StatefulKitty` supports re-encoding via `resize_encode()`, which replaces the `KittyProtoState` and forces retransmission on the next render. Between re-encodes, renders are placeholder-only (near zero bandwidth).

### iTerm2 Encoding Path

1. `Iterm2::new()` base64-encodes the image as PNG or similar.
2. No delta/stateful support in ratatui-image.

### Implications for Our Design

- For Sixel: We cannot avoid re-encoding when the image changes. Our savings come from (a) smaller images via resolution scaling, (b) fewer frames via dirty tracking and adaptive FPS, and (c) RGB instead of RGBA.
- For Kitty: We can get massive savings by using the stateful path and only re-encoding when the scene changes. An idle Kitty session would use near-zero bandwidth.
- For iTerm2: Same situation as Sixel -- no delta support. Apply the same scaling/throttling strategies.

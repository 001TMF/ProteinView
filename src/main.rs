mod app;
mod event;
mod model;
mod panel_server;
mod parser;
mod render;
mod ui;

use anyhow::Result;
use clap::Parser;
use crossterm::{
    event::KeyCode,
    execute,
    terminal::{EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode, enable_raw_mode},
};
use ratatui::layout::{Constraint, Direction, Layout};
use ratatui::prelude::*;
use std::io;
use std::path::PathBuf;
use std::sync::atomic::Ordering;
use std::time::{Duration, Instant};

use app::{App, AppConfig, ConnectionType, RenderMode, VizMode};
use model::selection::{ResidueColorSpec, resolve_residue_colors};

macro_rules! log {
    ($file:expr, $($arg:tt)*) => {
        if let Some(f) = $file.as_mut() {
            use std::io::Write;
            let _ = writeln!(f, $($arg)*);
            let _ = f.flush();
        }
    };
}

/// Terminal protein structure viewer
#[derive(Parser)]
#[command(name = "proteinview", version, about = "TUI protein structure viewer")]
struct Cli {
    /// Path to PDB, mmCIF, or XYZ file
    file: Option<String>,

    /// Use HD rendering (HalfBlock over SSH, FullHD locally)
    #[arg(long)]
    hd: bool,

    /// Force full pixel graphics (Sixel/Kitty/iTerm2) regardless of SSH
    #[arg(long, alias = "pixel")]
    fullhd: bool,

    /// Render mode: braille, halfblock (or hd), hdplus (or hd+), fullhd (or pixel)
    #[arg(long = "render", value_name = "MODE")]
    render_mode: Option<String>,

    /// Color scheme: structure, element, chain, plddt, bfactor, rainbow
    #[arg(long, default_value = "structure")]
    color: String,

    /// Override one exact polymer residue color: CHAIN:RES[ICODE]=RRGGBB
    #[arg(long, value_name = "SELECTOR=RRGGBB")]
    residue_color: Vec<ResidueColorSpec>,

    /// Palette file (TOML). Defaults to ~/.config/proteinview/palette.toml when present
    #[arg(long, value_name = "FILE")]
    palette: Option<PathBuf>,

    /// Visualization mode: cartoon, backbone, wireframe
    #[arg(long, default_value = "cartoon")]
    mode: String,

    /// Fetch structure from RCSB PDB by ID
    #[arg(long)]
    fetch: Option<String>,

    /// Run the persistent headless panel server over NDJSON stdin/stdout
    #[arg(long, requires = "output", conflicts_with = "snapshot")]
    panel_server: bool,

    /// PNG path atomically replaced by the panel server
    #[arg(long, value_name = "PNG", requires = "panel_server")]
    output: Option<PathBuf>,

    /// Initial panel framebuffer width in pixels
    #[arg(long, default_value_t = 960, value_name = "PX")]
    panel_width: u32,

    /// Initial panel framebuffer height in pixels
    #[arg(long, default_value_t = 540, value_name = "PX")]
    panel_height: u32,

    /// Render one FullHD pixel frame to a PNG and exit without starting the TUI
    #[arg(long, value_name = "PNG", conflicts_with = "panel_server")]
    snapshot: Option<PathBuf>,

    /// Snapshot width in pixels
    #[arg(
        long,
        default_value_t = render::snapshot::DEFAULT_SNAPSHOT_WIDTH,
        value_name = "PX"
    )]
    snapshot_width: u32,

    /// Snapshot height in pixels
    #[arg(
        long,
        default_value_t = render::snapshot::DEFAULT_SNAPSHOT_HEIGHT,
        value_name = "PX"
    )]
    snapshot_height: u32,

    /// Highlight the interface for this focus-chain ID in a snapshot
    #[arg(long, value_name = "CHAIN", requires = "snapshot")]
    snapshot_interface_chain: Option<String>,

    /// Overlay classified interface interaction lines in a snapshot
    #[arg(long, requires = "snapshot_interface_chain", requires = "snapshot")]
    snapshot_interactions: bool,

    /// Hide ligands and ions in a snapshot
    #[arg(long, requires = "snapshot")]
    snapshot_hide_ligands: bool,

    /// Write debug log to file (e.g. --log debug.log)
    #[arg(long)]
    log: Option<String>,

    /// Number of render threads (default: one per core)
    #[arg(long)]
    threads: Option<usize>,
}

/// Rebuild the sequence panel's wrapped layout for the current terminal size.
///
/// Cheap unless the width actually changed, so it is called both before and
/// after input handling.
fn sync_sequence_viewport(app: &mut App, fallback: (u16, u16)) {
    if !app.show_sequence {
        return;
    }
    let (cols, rows) = crossterm::terminal::size().unwrap_or(fallback);
    let panel_width = if app.show_interface {
        cols.saturating_sub(ui::interface_panel::SIDEBAR_WIDTH)
    } else {
        cols
    };
    let panel_height = ui::sequence_panel::height_for(app, rows);
    app.set_sequence_viewport(panel_width, panel_height);
}

/// Keys the sequence panel owns while it is open.
///
/// Returns `true` when the key was consumed, so anything the panel does not
/// claim still reaches the normal bindings and the camera stays live.
fn handle_sequence_key(app: &mut App, key: crossterm::event::KeyEvent) -> bool {
    use crossterm::event::KeyModifiers;

    // Shift plus an arrow extends the selection from where the cursor was, the
    // way a text editor does.
    let extend = key.modifiers.contains(KeyModifiers::SHIFT);
    match key.code {
        KeyCode::Left => app.seq_move_horizontal(-1, extend),
        KeyCode::Right => app.seq_move_horizontal(1, extend),
        KeyCode::Up => app.seq_move_vertical(-1, extend),
        KeyCode::Down => app.seq_move_vertical(1, extend),
        KeyCode::PageUp => app.seq_page(false, extend),
        KeyCode::PageDown => app.seq_page(true, extend),
        KeyCode::Home => app.seq_move_to_chain_edge(false, extend),
        KeyCode::End => app.seq_move_to_chain_edge(true, extend),
        KeyCode::Enter => app.seq_toggle_selection(),
        KeyCode::Char('A') => app.seq_toggle_chain_selection(),
        KeyCode::Char('x') => app.clear_selection(),
        KeyCode::Char('<') | KeyCode::Char(',') => app.resize_sequence_panel(-1),
        KeyCode::Char('>') | KeyCode::Char('.') => app.resize_sequence_panel(1),
        _ => return false,
    }
    true
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Resolve the color palette before anything renders.  A bad palette is a
    // hard error rather than a silent fallback: ignoring a file the user wrote
    // is worse than refusing to start.
    render::palette::init(cli.palette.as_deref())?;

    // Rasterization splits the framebuffer into bands, so it scales with cores
    // until it becomes memory-bound.  Default to one thread per core: with the
    // shared-memory transport the terminal no longer has a frame to decompress
    // on every tick, so there is no longer a reason to hold cores back for it.
    // Cap at 16 -- beyond that the bands get too thin to be worth a thread.
    let num_threads = cli.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map_or(4, std::num::NonZeroUsize::get)
            .min(16)
    });
    let num_threads = num_threads.max(1);
    match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
    {
        Ok(()) => {}
        Err(e) => eprintln!("Warning: failed to initialize rayon thread pool: {e}"),
    }

    // Determine the file path. Fetched structures are temporary inputs and are
    // removed immediately after parsing, including when parsing fails.
    let fetched_path = cli
        .fetch
        .as_deref()
        .map(parser::fetch::fetch_pdb)
        .transpose()?;
    let file_path = if let Some(path) = fetched_path.as_ref() {
        path.clone()
    } else if let Some(path) = &cli.file {
        path.clone()
    } else {
        eprintln!("Error: provide a file path or use --fetch <PDB_ID>");
        std::process::exit(1);
    };

    // Load protein structure (dispatch by file extension)
    let lower = file_path.to_lowercase();
    let is_xyz = lower.ends_with(".xyz");
    let protein_result = if is_xyz {
        parser::xyz::load_xyz(&file_path)
    } else {
        parser::pdb::load_structure(&file_path)
    };
    if let Some(path) = fetched_path {
        if let Err(error) = std::fs::remove_file(&path) {
            eprintln!(
                "Warning: failed to remove temporary fetched structure '{}': {error}",
                path
            );
        }
    }
    let protein = protein_result?;
    let residue_colors = resolve_residue_colors(&protein, &cli.residue_color)?;
    eprintln!(
        "Loaded: {} ({} chains, {} residues, {} atoms{})",
        protein.name,
        protein.chains.len(),
        protein.residue_count(),
        protein.atom_count(),
        if protein.ligands.is_empty() {
            String::new()
        } else {
            format!(", {} ligands", protein.ligand_count())
        },
    );

    // Open log file if requested
    let mut logfile: Option<std::fs::File> = match &cli.log {
        Some(path) => {
            let f = std::fs::File::create(path)
                .map_err(|e| anyhow::anyhow!("cannot create log file '{}': {}", path, e))?;
            Some(f)
        }
        None => None,
    };

    // Detect connection type
    let connection_type = ConnectionType::detect();
    log!(logfile, "connection type: {:?}", connection_type);

    // Determine render mode from CLI flags
    let render_mode = if let Some(mode_str) = &cli.render_mode {
        match mode_str.to_ascii_lowercase().as_str() {
            "braille" => RenderMode::Braille,
            "halfblock" | "hd" | "half-block" => RenderMode::HalfBlock,
            "hdplus" | "hd+" | "halfblockplus" | "half-block-plus" => RenderMode::HalfBlockPlus,
            "fullhd" | "pixel" | "full-hd" => RenderMode::FullHD,
            _ => {
                eprintln!("Warning: unknown render mode '{}', using default", mode_str);
                RenderMode::Braille
            }
        }
    } else if cli.fullhd {
        // --fullhd / --pixel always forces FullHD regardless of SSH
        RenderMode::FullHD
    } else if cli.hd {
        // --hd is SSH-aware: FullHD locally, HalfBlock over SSH
        match connection_type {
            ConnectionType::Local => RenderMode::FullHD,
            ConnectionType::Ssh => RenderMode::HalfBlock,
        }
    } else {
        RenderMode::Braille
    };

    let user_explicit_color =
        std::env::args().any(|argument| argument == "--color" || argument.starts_with("--color="));

    // Parse CLI color scheme override.
    let color_override = match cli.color.to_ascii_lowercase().as_str() {
        "structure" => None, // default, no override needed
        "element" => Some(render::color::ColorSchemeType::Element),
        "chain" => Some(render::color::ColorSchemeType::Chain),
        "bfactor" | "b-factor" => Some(render::color::ColorSchemeType::BFactor),
        "rainbow" => Some(render::color::ColorSchemeType::Rainbow),
        "plddt" => Some(render::color::ColorSchemeType::Plddt),
        _ => {
            eprintln!(
                "Warning: unknown color scheme '{}', using structure",
                cli.color
            );
            None
        }
    };

    // Parse CLI visualization mode override.
    let user_explicit_mode = !cli.mode.eq_ignore_ascii_case("cartoon")
        || std::env::args().any(|a| a == "--mode" || a.starts_with("--mode="));
    let viz_mode = match cli.mode.to_ascii_lowercase().as_str() {
        "cartoon" => VizMode::Cartoon,
        "backbone" => VizMode::Backbone,
        "wireframe" => VizMode::Wireframe,
        _ => {
            eprintln!(
                "Warning: unknown visualization mode '{}', using cartoon",
                cli.mode
            );
            VizMode::Cartoon
        }
    };

    // XYZ files default to Element coloring + Wireframe mode unless overridden.
    let (color_override, viz_mode) = if is_xyz {
        let color = if color_override.is_none() && !user_explicit_color {
            Some(render::color::ColorSchemeType::Element)
        } else {
            color_override
        };
        let viz = if !user_explicit_mode {
            VizMode::Wireframe
        } else {
            viz_mode
        };
        (color, viz)
    } else {
        (color_override, viz_mode)
    };

    // The panel server owns a persistent, terminal-independent render session.
    // It must start before raw mode, terminal probing, alternate-screen setup,
    // or the crossterm input thread so its stdout remains strict NDJSON.
    if cli.panel_server {
        let output_path = cli
            .output
            .as_deref()
            .expect("clap enforces --output with --panel-server");
        let stdin = io::stdin();
        let stdout = io::stdout();
        panel_server::serve(
            protein,
            output_path,
            panel_server::PanelServerOptions {
                width: cli.panel_width,
                height: cli.panel_height,
                color_override,
                residue_colors: residue_colors.clone(),
                viz_mode,
                user_explicit_mode,
            },
            stdin.lock(),
            stdout.lock(),
        )?;
        return Ok(());
    }

    // A snapshot is the non-interactive form of ProteinView's FullHD pixel
    // renderer. Exit before raw mode, terminal probing, or alternate-screen
    // setup so it is safe to call from an agent tool running inside a TUI.
    if let Some(snapshot_path) = cli.snapshot.as_deref() {
        if cli.snapshot_interface_chain.is_some() && user_explicit_color {
            anyhow::bail!(
                "--snapshot-interface-chain uses ProteinView's interface palette and cannot be combined with --color"
            );
        }
        render::snapshot::save_png(
            protein,
            snapshot_path,
            render::snapshot::SnapshotOptions {
                width: cli.snapshot_width,
                height: cli.snapshot_height,
                color_override,
                residue_colors: residue_colors.clone(),
                viz_mode,
                user_explicit_mode,
                show_ligands: !cli.snapshot_hide_ligands,
                interface_chain: cli.snapshot_interface_chain.clone(),
                show_interactions: cli.snapshot_interactions,
            },
        )?;
        eprintln!(
            "Rendered FullHD PNG: {} ({}x{})",
            snapshot_path.display(),
            cli.snapshot_width,
            cli.snapshot_height
        );
        return Ok(());
    }

    // Get terminal dimensions before entering alternate screen
    let (term_cols, term_rows) = crossterm::terminal::size().unwrap_or((80, 24));
    log!(logfile, "terminal size: {}x{}", term_cols, term_rows);

    // Install panic hook that restores the terminal
    let original_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = disable_raw_mode();
        let _ = execute!(io::stderr(), LeaveAlternateScreen);
        render::kitty_shm::unlink_all();
        original_hook(info);
    }));

    // Setup terminal — must happen before Picker::from_query_stdio()
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    // Detect terminal graphics protocol (Sixel/Kitty/iTerm2) and font size.
    // Must be called after entering alternate screen but before spawning the
    // input thread (which reads from stdin).
    let picker = ratatui_image::picker::Picker::from_query_stdio()
        .unwrap_or_else(|_| ratatui_image::picker::Picker::halfblocks());
    log!(
        logfile,
        "picker: protocol={:?} font_size={:?}",
        picker.protocol_type(),
        picker.font_size()
    );

    // Ask the terminal whether it can read pixels straight out of shared
    // memory.  Must happen here: after the picker's own query has drained its
    // responses, and before the input thread starts consuming stdin.  Over SSH
    // there is no shared filesystem to share memory through, so don't ask.
    let kitty_shm = picker.protocol_type() == ratatui_image::picker::ProtocolType::Kitty
        && connection_type == ConnectionType::Local
        && render::kitty_shm::probe(Duration::from_millis(500));
    log!(logfile, "kitty shared-memory transport: {}", kitty_shm);

    // Create app with actual terminal dimensions for dynamic zoom
    let mut app = App::new(
        protein,
        AppConfig {
            render_mode,
            viz_mode,
            user_explicit_mode,
            color_override,
            residue_colors,
        },
        term_cols,
        term_rows,
        picker,
    );
    app.kitty_shm = kitty_shm;
    log!(
        logfile,
        "app created: render_mode={:?} chains={} zoom={:.2}",
        app.render_mode,
        app.protein.chains.len(),
        app.camera.zoom
    );

    // Spawn dedicated input thread — decouples input from rendering so
    // quit always works even when HD rendering is slow
    let (input_rx, quit_flag) = event::spawn_input_thread();

    // Main loop
    let tick_rate = Duration::from_millis(33); // ~30 FPS
    let mut frame_count: u64 = 0;
    // Track how long the previous terminal.draw() took so we can skip frames
    // when rendering is too slow (prevents PTY buffer saturation & freezes).
    let mut last_draw_duration = Duration::ZERO;
    let mut frames_to_skip: u32 = 0;
    // Tracks the interaction state of the previous drawn frame, so the
    // transition back to a still view can trigger one full-resolution redraw.
    let mut was_interacting = false;

    loop {
        // Cursor movement is expressed in layout rows, so the layout has to be
        // current before keys are handled.
        sync_sequence_viewport(&mut app, (term_cols, term_rows));

        // Drain all queued input from the dedicated input thread
        let mut had_input = false;
        while let Ok(app_event) = input_rx.try_recv() {
            had_input = true;
            match app_event {
                event::AppEvent::Resize(cols, rows) => {
                    log!(logfile, "resize: {}x{}", cols, rows);
                    app.recalculate_zoom(cols, rows);
                    app.mesh_dirty_flag();
                }
                event::AppEvent::Key(key) => {
                    log!(logfile, "key: {:?}", key.code);
                    // While the panel is open it owns the arrow keys and the
                    // selection keys; h/j/k/l keep driving the camera, so the
                    // view stays steerable with the sequence in front of you.
                    if app.show_sequence && handle_sequence_key(&mut app, key) {
                        continue;
                    }
                    match key.code {
                        KeyCode::Char('q') => app.should_quit = true,
                        KeyCode::Char('c')
                            if key
                                .modifiers
                                .contains(crossterm::event::KeyModifiers::CONTROL) =>
                        {
                            app.should_quit = true
                        }
                        // Every camera-moving key notes an interaction, which
                        // drops FullHD to its reduced resolution until the view
                        // settles again.
                        KeyCode::Char('h') | KeyCode::Left => {
                            app.camera.rotate_y(-1.0);
                            app.note_camera_change();
                        }
                        KeyCode::Char('l') | KeyCode::Right => {
                            app.camera.rotate_y(1.0);
                            app.note_camera_change();
                        }
                        KeyCode::Char('j') | KeyCode::Down => {
                            app.camera.rotate_x(1.0);
                            app.note_camera_change();
                        }
                        KeyCode::Char('k') | KeyCode::Up => {
                            app.camera.rotate_x(-1.0);
                            app.note_camera_change();
                        }
                        KeyCode::Char('u') => {
                            app.camera.rotate_z(-1.0);
                            app.note_camera_change();
                        }
                        KeyCode::Char('i') => {
                            app.camera.rotate_z(1.0);
                            app.note_camera_change();
                        }
                        KeyCode::Char('+') | KeyCode::Char('=') => {
                            app.camera.zoom_in();
                            app.note_camera_change();
                        }
                        KeyCode::Char('-') => {
                            app.camera.zoom_out();
                            app.note_camera_change();
                        }
                        KeyCode::Char('w') => {
                            app.camera.pan(0.0, 1.0);
                            app.note_camera_change();
                        }
                        KeyCode::Char('s') => {
                            app.camera.pan(0.0, -1.0);
                            app.note_camera_change();
                        }
                        KeyCode::Char('a') => {
                            app.camera.pan(-1.0, 0.0);
                            app.note_camera_change();
                        }
                        KeyCode::Char('d') => {
                            app.camera.pan(1.0, 0.0);
                            app.note_camera_change();
                        }
                        KeyCode::Char('r') => {
                            let (cols, rows) =
                                crossterm::terminal::size().unwrap_or((term_cols, term_rows));
                            app.camera.reset();
                            app.recalculate_zoom(cols, rows);
                            app.note_camera_change();
                        }
                        KeyCode::Char('c') => app.cycle_color(),
                        KeyCode::Char('v') => app.cycle_viz_mode(),
                        KeyCode::Char('m') => {
                            let (cols, rows) =
                                crossterm::terminal::size().unwrap_or((term_cols, term_rows));
                            app.toggle_hd(cols, rows);
                        }
                        KeyCode::Char('M') => {
                            let (cols, rows) =
                                crossterm::terminal::size().unwrap_or((term_cols, term_rows));
                            app.toggle_fullhd(cols, rows);
                        }
                        KeyCode::Char('[') => app.prev_chain(),
                        KeyCode::Char(']') => app.next_chain(),
                        KeyCode::Char(' ') => app.camera.auto_rotate = !app.camera.auto_rotate,
                        KeyCode::Char('f') => app.toggle_interface(),
                        KeyCode::Char('I') => app.toggle_interactions(),
                        KeyCode::Char('g') => app.toggle_ligands(),
                        KeyCode::Char('?') => app.show_help = !app.show_help,
                        KeyCode::Char('S') => app.toggle_sequence_panel(),
                        KeyCode::Char('b') => app.toggle_ball_stick(),
                        KeyCode::Char('z') => {
                            app.focus_on_selection();
                        }
                        KeyCode::Esc => {
                            if app.show_help {
                                app.show_help = false;
                            } else if app.show_sequence {
                                app.show_sequence = false;
                            }
                        }
                        _ => {}
                    }
                }
            }
        }

        if app.should_quit {
            break;
        }

        // Ensure ribbon mesh cache is fresh (rebuilds only when color scheme changes).
        // Must happen outside terminal.draw() since ribbon_mesh() needs &mut self.
        // Only rebuild when in Cartoon mode — Backbone/Wireframe don't use the
        // ribbon mesh, so skipping this preserves the lazy-mesh optimization for
        // large structures that start in a non-Cartoon mode.
        let mesh_was_rebuilt = app.viz_mode == VizMode::Cartoon && app.mesh_is_dirty();
        if app.viz_mode == VizMode::Cartoon {
            app.ribbon_mesh();
        }

        // Always poll the background interface thread, even during skipped
        // frames, so the result is absorbed as soon as it's available.
        let interface_was_pending = app.interface_pending();
        app.poll_background_interface();
        let interface_absorbed = interface_was_pending && !app.interface_pending();

        // Adaptive frame skipping: if the previous draw took longer than the
        // tick rate, skip frames proportionally.  User input always forces a
        // redraw so the UI stays responsive.
        //
        // Do NOT call app.tick() during skipped frames — that would advance
        // auto-rotate without a corresponding render, causing the protein to
        // "jump" when rendering resumes.  Instead we just sleep and let the
        // camera's dt-clamping handle the gap.
        if frames_to_skip > 0 && !had_input {
            frames_to_skip -= 1;
            // Reset the camera's tick timer so the next real tick doesn't see
            // a huge accumulated dt from the skipped frames.
            app.camera.reset_tick_timer();
            std::thread::sleep(tick_rate);
            continue;
        }

        // Nothing on screen changes unless input arrived, an animation is
        // running, or background state was just absorbed.  Redrawing anyway
        // would re-run the whole rasterize + transmit pipeline at every tick,
        // forever, for an image identical to the one already on screen.
        //
        // The one extra case is `settled`: FullHD renders at a reduced
        // resolution while the camera moves, so the frame after it comes to
        // rest must be drawn to replace it with the sharp one.
        let interacting = app.is_interacting();
        let settled = was_interacting && !interacting;
        was_interacting = interacting;

        let must_redraw = had_input
            || interacting
            || settled
            || app.ssh_hd_warning
            || app.needs_clear
            || mesh_was_rebuilt
            || interface_absorbed
            || frame_count < 2;
        if !must_redraw {
            app.tick();
            std::thread::sleep(tick_rate);
            continue;
        }

        // Render
        frame_count += 1;
        if frame_count <= 3 || frame_count % 300 == 0 {
            log!(
                logfile,
                "frame {} render start (render_mode={:?} viz={:?} interface={} last_draw={:?})",
                frame_count,
                app.render_mode,
                app.viz_mode,
                app.show_interface,
                last_draw_duration
            );
        }

        // After a render-mode switch, force ratatui to redraw every cell.
        // Without this, its diff-based rendering may leave stale characters
        // from the previous mode (e.g. braille dots under a FullHD image).
        if app.needs_clear {
            // Delete any Kitty graphics images that may be lingering from
            // a previous FullHD session.  Harmless no-op if there are none.
            let cleanup = render::kitty_png::KittyPngImage::cleanup_escape();
            execute!(terminal.backend_mut(), crossterm::style::Print(&cleanup))?;
            // Drop any shared memory object the terminal never got round to
            // reading, so switching modes cannot leave objects behind.
            render::kitty_shm::unlink_all();
            terminal.clear()?;
            app.needs_clear = false;
        }

        // ...and again after input, because the key that opened the panel or
        // resized it arrived after the first sync: the frame about to be drawn
        // must not be the one frame with a stale layout.
        sync_sequence_viewport(&mut app, (term_cols, term_rows));

        let draw_start = Instant::now();
        terminal.draw(|frame| {
            // If interface is active, split horizontally: sidebar | main
            let main_area = if app.show_interface {
                let horiz = Layout::default()
                    .direction(Direction::Horizontal)
                    .constraints([
                        Constraint::Length(ui::interface_panel::SIDEBAR_WIDTH),
                        Constraint::Min(20),
                    ])
                    .split(frame.area());

                let summary = app.interface_analysis.summary(&app.protein);
                let chain_names = app.chain_names();
                let interaction_counts = app.interface_analysis.interaction_counts();
                ui::interface_panel::render_interface_panel(
                    frame,
                    horiz[0],
                    &summary,
                    app.current_chain,
                    &chain_names,
                    app.show_interactions,
                    interaction_counts,
                );
                horiz[1]
            } else {
                frame.area()
            };

            let sequence_height = ui::sequence_panel::height_for(&app, main_area.height);
            let chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Length(1),               // Header
                    Constraint::Min(3),                  // Viewport
                    Constraint::Length(sequence_height), // Sequence panel
                    Constraint::Length(2),               // Status bar
                    Constraint::Length(1),               // Help bar
                ])
                .split(main_area);

            ui::header::render_header(frame, chunks[0], &app.protein.name);
            ui::viewport::render_viewport(frame, chunks[1], &app);
            if sequence_height > 0 {
                ui::sequence_panel::render_sequence_panel(frame, chunks[2], &app);
            }
            ui::statusbar::render_statusbar(frame, chunks[3], &app);
            ui::helpbar::render_helpbar(frame, chunks[4], &app);

            if app.show_help {
                ui::help_overlay::render_help_overlay(frame, frame.area());
            }
        })?;
        last_draw_duration = draw_start.elapsed();

        // If the draw took longer than two tick periods, skip some frames to
        // let the terminal catch up and avoid saturating the PTY write buffer.
        if last_draw_duration > tick_rate * 2 {
            // Skip 1-3 frames depending on how slow the draw was.
            frames_to_skip = ((last_draw_duration.as_millis() / tick_rate.as_millis()) as u32)
                .saturating_sub(1)
                .min(3);
        }

        app.tick();

        // Sleep for the remainder of the tick period to cap at ~30 FPS.
        // Account for the time already spent drawing so the frame rate stays
        // consistent regardless of render cost.
        let elapsed = draw_start.elapsed();
        if let Some(remaining) = tick_rate.checked_sub(elapsed) {
            std::thread::sleep(remaining);
        }
    }

    // Signal input thread to stop
    quit_flag.store(true, Ordering::Relaxed);

    // Restore terminal
    render::kitty_shm::unlink_all();
    disable_raw_mode()?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen)?;
    terminal.show_cursor()?;

    Ok(())
}

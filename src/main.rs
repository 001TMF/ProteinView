mod app;
mod event;
mod model;
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
use std::sync::atomic::Ordering;
use std::time::{Duration, Instant};

use app::{App, AppConfig, ConnectionType, RenderMode, VizMode};

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

    /// Render mode: braille, halfblock (or hd), fullhd (or pixel)
    #[arg(long = "render", value_name = "MODE")]
    render_mode: Option<String>,

    /// Color scheme: structure, element, chain, plddt, bfactor, rainbow
    #[arg(long, default_value = "structure")]
    color: String,

    /// Visualization mode: cartoon, backbone, wireframe
    #[arg(long, default_value = "cartoon")]
    mode: String,

    /// Fetch structure from RCSB PDB by ID
    #[arg(long)]
    fetch: Option<String>,

    /// Write debug log to file (e.g. --log debug.log)
    #[arg(long)]
    log: Option<String>,

    /// Number of render threads (default: 4)
    #[arg(long, default_value = "4")]
    threads: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Cap rayon thread pool. 4 threads is the sweet spot: the framebuffer
    // only has ~60 tiles (64x64) so more threads hit diminishing returns,
    // and 4 leaves cores free for the terminal emulator and OS.
    let num_threads = cli.threads.max(1);
    match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
    {
        Ok(()) => {}
        Err(e) => eprintln!("Warning: failed to initialize rayon thread pool: {e}"),
    }

    // Determine the file path
    let file_path = if let Some(pdb_id) = &cli.fetch {
        parser::fetch::fetch_pdb(pdb_id)?
    } else if let Some(path) = &cli.file {
        path.clone()
    } else {
        eprintln!("Error: provide a file path or use --fetch <PDB_ID>");
        std::process::exit(1);
    };

    // Load protein structure (dispatch by file extension)
    let lower = file_path.to_lowercase();
    let is_xyz = lower.ends_with(".xyz");
    let protein = if is_xyz {
        parser::xyz::load_xyz(&file_path)?
    } else {
        parser::pdb::load_structure(&file_path)?
    };
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

    // Get terminal dimensions before entering alternate screen
    let (term_cols, term_rows) = crossterm::terminal::size().unwrap_or((80, 24));
    log!(logfile, "terminal size: {}x{}", term_cols, term_rows);

    // Install panic hook that restores the terminal
    let original_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = disable_raw_mode();
        let _ = execute!(io::stderr(), LeaveAlternateScreen);
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

    // Parse CLI color scheme override
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

    // Parse CLI visualization mode override
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

    // XYZ files default to Element coloring + Wireframe mode unless overridden
    let (color_override, viz_mode) = if is_xyz {
        let color = if color_override.is_none() && cli.color == "structure" {
            Some(render::color::ColorSchemeType::Element)
        } else {
            color_override
        };
        let viz = if cli.mode == "cartoon" {
            VizMode::Wireframe
        } else {
            viz_mode
        };
        (color, viz)
    } else {
        (color_override, viz_mode)
    };

    // Create app with actual terminal dimensions for dynamic zoom
    let mut app = App::new(
        protein,
        AppConfig {
            render_mode,
            viz_mode,
            user_explicit_mode,
            color_override,
        },
        term_cols,
        term_rows,
        picker,
    );
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

    loop {
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
                    match key.code {
                        KeyCode::Char('q') => app.should_quit = true,
                        KeyCode::Char('c')
                            if key
                                .modifiers
                                .contains(crossterm::event::KeyModifiers::CONTROL) =>
                        {
                            app.should_quit = true
                        }
                        KeyCode::Char('h') | KeyCode::Left => app.camera.rotate_y(-1.0),
                        KeyCode::Char('l') | KeyCode::Right => app.camera.rotate_y(1.0),
                        KeyCode::Char('j') | KeyCode::Down => app.camera.rotate_x(1.0),
                        KeyCode::Char('k') | KeyCode::Up => app.camera.rotate_x(-1.0),
                        KeyCode::Char('u') => app.camera.rotate_z(-1.0),
                        KeyCode::Char('i') => app.camera.rotate_z(1.0),
                        KeyCode::Char('+') | KeyCode::Char('=') => app.camera.zoom_in(),
                        KeyCode::Char('-') => app.camera.zoom_out(),
                        KeyCode::Char('w') => app.camera.pan(0.0, 1.0),
                        KeyCode::Char('s') => app.camera.pan(0.0, -1.0),
                        KeyCode::Char('a') => app.camera.pan(-1.0, 0.0),
                        KeyCode::Char('d') => app.camera.pan(1.0, 0.0),
                        KeyCode::Char('r') => {
                            let (cols, rows) =
                                crossterm::terminal::size().unwrap_or((term_cols, term_rows));
                            app.camera.reset();
                            app.recalculate_zoom(cols, rows);
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
                        KeyCode::Esc => {
                            if app.show_help {
                                app.show_help = false;
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
        if app.viz_mode == VizMode::Cartoon {
            app.ribbon_mesh();
        }

        // Always poll the background interface thread, even during skipped
        // frames, so the result is absorbed as soon as it's available.
        app.poll_background_interface();

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
            terminal.clear()?;
            app.needs_clear = false;
        }

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

            let chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Length(1), // Header
                    Constraint::Min(3),    // Viewport
                    Constraint::Length(2), // Status bar
                    Constraint::Length(1), // Help bar
                ])
                .split(main_area);

            ui::header::render_header(frame, chunks[0], &app.protein.name);
            ui::viewport::render_viewport(frame, chunks[1], &app);
            ui::statusbar::render_statusbar(frame, chunks[2], &app);
            ui::helpbar::render_helpbar(frame, chunks[3]);

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
    disable_raw_mode()?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen)?;
    terminal.show_cursor()?;

    Ok(())
}

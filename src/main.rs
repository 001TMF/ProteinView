mod app;
mod event;
mod model;
mod parser;
mod render;
mod ui;

use std::io;
use std::time::Duration;
use anyhow::Result;
use clap::Parser;
use crossterm::{
    event::KeyCode,
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::prelude::*;
use ratatui::layout::{Constraint, Direction, Layout};

use app::App;

/// Terminal protein structure viewer
#[derive(Parser)]
#[command(name = "proteinview", version, about = "TUI protein structure viewer")]
struct Cli {
    /// Path to PDB or mmCIF file
    file: Option<String>,

    /// Use HD pixel rendering (sixel/kitty)
    #[arg(long, alias = "pixel")]
    hd: bool,

    /// Color scheme: structure, chain, element, bfactor, rainbow
    #[arg(long, default_value = "structure")]
    color: String,

    /// Visualization mode: backbone, ballandstick, wireframe
    #[arg(long, default_value = "backbone")]
    mode: String,

    /// Fetch structure from RCSB PDB by ID
    #[arg(long)]
    fetch: Option<String>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Determine the file path
    let file_path = if let Some(pdb_id) = &cli.fetch {
        parser::fetch::fetch_pdb(pdb_id)?
    } else if let Some(path) = &cli.file {
        path.clone()
    } else {
        eprintln!("Error: provide a file path or use --fetch <PDB_ID>");
        std::process::exit(1);
    };

    // Load protein structure
    let protein = parser::pdb::load_structure(&file_path)?;
    eprintln!("Loaded: {} ({} chains, {} residues, {} atoms)",
        protein.name,
        protein.chains.len(),
        protein.residue_count(),
        protein.atom_count(),
    );

    // Create app
    let mut app = App::new(protein, cli.hd);

    // Setup terminal
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    // Main loop
    let tick_rate = Duration::from_millis(33); // ~30 FPS
    loop {
        // Render
        terminal.draw(|frame| {
            let chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Length(1),     // Header
                    Constraint::Min(3),        // Viewport
                    Constraint::Length(2),      // Status bar
                    Constraint::Length(1),      // Help bar
                ])
                .split(frame.area());

            ui::header::render_header(frame, chunks[0], &app.protein.name);
            ui::viewport::render_viewport(frame, chunks[1], &app);
            ui::statusbar::render_statusbar(frame, chunks[2], &app);
            ui::helpbar::render_helpbar(frame, chunks[3]);

            if app.show_help {
                ui::help_overlay::render_help_overlay(frame, frame.area());
            }
        })?;

        if app.should_quit { break; }

        // Handle input
        if let Some(key) = event::poll_event(tick_rate)? {
            match key.code {
                KeyCode::Char('q') => app.should_quit = true,
                KeyCode::Char('h') => app.camera.rotate_y(-1.0),
                KeyCode::Char('l') => app.camera.rotate_y(1.0),
                KeyCode::Char('j') => app.camera.rotate_x(1.0),
                KeyCode::Char('k') => app.camera.rotate_x(-1.0),
                KeyCode::Char('u') => app.camera.rotate_z(-1.0),
                KeyCode::Char('i') => app.camera.rotate_z(1.0),
                KeyCode::Char('+') | KeyCode::Char('=') => app.camera.zoom_in(),
                KeyCode::Char('-') => app.camera.zoom_out(),
                KeyCode::Char('w') => app.camera.pan(0.0, 1.0),
                KeyCode::Char('s') => app.camera.pan(0.0, -1.0),
                KeyCode::Char('a') => app.camera.pan(-1.0, 0.0),
                KeyCode::Char('d') => app.camera.pan(1.0, 0.0),
                KeyCode::Char('r') => app.camera.reset(),
                KeyCode::Char('c') => app.cycle_color(),
                KeyCode::Char('v') => app.cycle_viz_mode(),
                KeyCode::Char('m') => app.hd_mode = !app.hd_mode,
                KeyCode::Char('[') => app.prev_chain(),
                KeyCode::Char(']') => app.next_chain(),
                KeyCode::Char(' ') => app.camera.auto_rotate = !app.camera.auto_rotate,
                KeyCode::Char('?') => app.show_help = !app.show_help,
                KeyCode::Esc => { if app.show_help { app.show_help = false; } },
                _ => {}
            }
        }

        app.tick();
    }

    // Restore terminal
    disable_raw_mode()?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen)?;
    terminal.show_cursor()?;

    Ok(())
}

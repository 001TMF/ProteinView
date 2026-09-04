//! The scrollable chain-sequence panel.
//!
//! Every chain in the structure is laid out as one-letter codes, wrapped to the
//! panel width and scrolled as a single list.  The cursor and the selection are
//! drawn straight onto the letters, so picking residues here is what drives the
//! ball-and-stick overlay in the 3D view.

use ratatui::Frame;
use ratatui::layout::Rect;
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph};

use crate::app::App;
use crate::model::protein::{Chain, MoleculeType, SecondaryStructure};
use crate::model::sequence::{GROUP, SeqRow, column_offset, molecule_label, one_letter};
use crate::render::palette::palette;

/// Width of the left gutter: a right-aligned residue number plus a space.
pub const GUTTER: u16 = 7;

/// Color of a picked residue, from the palette's `[selection]` section.
fn selection_color() -> Color {
    let [red, green, blue] = palette().selection.marker.0;
    Color::Rgb(red, green, blue)
}

/// Color of the cursor cell, from the palette's `[selection]` section.
fn cursor_color() -> Color {
    let [red, green, blue] = palette().selection.cursor.0;
    Color::Rgb(red, green, blue)
}

/// Height the panel takes in a layout `total_rows` tall, or 0 when closed.
///
/// The rest of the interface needs seven rows (header, a minimum viewport,
/// status bar, help bar), so the panel never squeezes the 3D view out of
/// existence however far `>` is held down.
pub fn height_for(app: &App, total_rows: u16) -> u16 {
    if !app.show_sequence {
        return 0;
    }
    const CHROME_ROWS: u16 = 7;
    app.seq_panel_height
        .min(total_rows.saturating_sub(CHROME_ROWS))
}

/// Render the sequence panel into `area`.
pub fn render_sequence_panel(frame: &mut Frame, area: Rect, app: &App) {
    if area.height < 2 {
        return;
    }

    let block = Block::default()
        .borders(Borders::TOP)
        .border_style(Style::default().fg(Color::DarkGray))
        .title(" Sequence ")
        .title_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    let inner = block.inner(area);
    frame.render_widget(block, area);

    if inner.height == 0 {
        return;
    }

    // Cursor line: what is under the cursor, and what the selection holds.
    frame.render_widget(
        Paragraph::new(cursor_line(app)),
        Rect::new(inner.x, inner.y, inner.width, 1),
    );

    let rows_area = Rect::new(
        inner.x,
        inner.y + 1,
        inner.width,
        inner.height.saturating_sub(1),
    );
    if rows_area.height == 0 {
        return;
    }

    let layout = app.sequence_layout();
    let mut lines: Vec<Line> = Vec::with_capacity(rows_area.height as usize);
    for offset in 0..rows_area.height as usize {
        let row = app.seq_scroll + offset;
        match layout.rows.get(row) {
            Some(SeqRow::Header(chain_index)) => {
                lines.push(header_line(app, *chain_index));
            }
            Some(SeqRow::Residues { chain, start, len }) => {
                lines.push(residue_line(app, *chain, *start, *len));
            }
            None => lines.push(Line::from("")),
        }
    }

    frame.render_widget(Paragraph::new(lines), rows_area);
}

/// `Lb 245 ARG  helix │ 12 residues in 3 chains │ ball&stick on`
fn cursor_line(app: &App) -> Line<'static> {
    let mut spans = vec![Span::styled(" ", Style::default())];

    match app.seq_cursor_residue() {
        Some((chain, residue)) => {
            let letter = one_letter(&residue.name, chain.molecule_type);
            let insertion = residue.insertion_code.as_deref().unwrap_or("");
            spans.push(Span::styled(
                format!("{} ", chain.id),
                Style::default()
                    .fg(Color::Cyan)
                    .add_modifier(Modifier::BOLD),
            ));
            spans.push(Span::styled(
                format!("{}{} ", residue.seq_num, insertion),
                Style::default().fg(Color::White),
            ));
            spans.push(Span::styled(
                format!("{} ({letter}) ", residue.name),
                Style::default().fg(Color::Gray),
            ));
            if chain.molecule_type == MoleculeType::Protein {
                spans.push(Span::styled(
                    format!("{} ", secondary_label(residue.secondary_structure)),
                    Style::default().fg(Color::DarkGray),
                ));
            }
            if app.selection.contains(app.seq_cursor.0, app.seq_cursor.1) {
                spans.push(Span::styled(
                    "[selected] ",
                    Style::default().fg(selection_color()),
                ));
            }
        }
        None => spans.push(Span::styled(
            "no residues ",
            Style::default().fg(Color::DarkGray),
        )),
    }

    spans.push(Span::styled("│ ", Style::default().fg(Color::DarkGray)));
    if app.selection.is_empty() {
        spans.push(Span::styled(
            "nothing selected ",
            Style::default().fg(Color::DarkGray),
        ));
    } else {
        spans.push(Span::styled(
            format!(
                "{} res in {} chain{} ",
                app.selection.count(),
                app.selection.chain_count(),
                if app.selection.chain_count() == 1 {
                    ""
                } else {
                    "s"
                }
            ),
            Style::default().fg(selection_color()),
        ));
        spans.push(Span::styled("│ ", Style::default().fg(Color::DarkGray)));
        spans.push(Span::styled(
            format!("{} ", app.selection.describe(&app.protein, 4)),
            Style::default().fg(Color::White),
        ));
        spans.push(Span::styled("│ ", Style::default().fg(Color::DarkGray)));
        spans.push(Span::styled(
            if app.show_ball_stick {
                "ball&stick on "
            } else {
                "ball&stick off "
            },
            Style::default().fg(if app.show_ball_stick {
                Color::Rgb(255, 170, 0)
            } else {
                Color::DarkGray
            }),
        ));
    }

    Line::from(spans)
}

fn secondary_label(secondary: SecondaryStructure) -> &'static str {
    match secondary {
        SecondaryStructure::Helix => "helix",
        SecondaryStructure::Sheet => "sheet",
        SecondaryStructure::Turn => "turn",
        SecondaryStructure::Coil => "coil",
    }
}

/// `> Lb  Protein  138 res  1-138`
fn header_line(app: &App, chain_index: usize) -> Line<'static> {
    let Some(chain) = app.protein.chains.get(chain_index) else {
        return Line::from("");
    };
    let range = match (chain.residues.first(), chain.residues.last()) {
        (Some(first), Some(last)) => format!("{}-{}", first.seq_num, last.seq_num),
        _ => "empty".to_string(),
    };
    let marker_style = if chain_index == app.current_chain {
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD)
    } else {
        Style::default().fg(Color::DarkGray)
    };
    let mut spans = vec![
        Span::styled(
            if chain_index == app.current_chain {
                "▸ "
            } else {
                "  "
            },
            marker_style,
        ),
        Span::styled(
            chain.id.clone(),
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(
            format!(
                "  {}  {} res  {range}",
                molecule_label(chain.molecule_type),
                chain.residues.len()
            ),
            Style::default().fg(Color::DarkGray),
        ),
    ];
    if app.selection.chain_has_any(chain_index) {
        spans.push(Span::styled("  ●", Style::default().fg(selection_color())));
    }
    Line::from(spans)
}

/// One wrapped run of residues, with its starting residue number in the gutter.
fn residue_line(app: &App, chain_index: usize, start: usize, len: usize) -> Line<'static> {
    let Some(chain) = app.protein.chains.get(chain_index) else {
        return Line::from("");
    };

    let number = chain
        .residues
        .get(start)
        .map(|residue| residue.seq_num.to_string())
        .unwrap_or_default();
    let mut spans = vec![Span::styled(
        format!("{number:>width$} ", width = GUTTER as usize - 1),
        Style::default().fg(Color::DarkGray),
    )];

    // Consecutive residues that share a style are emitted as one span: a
    // 200-column row of one-letter codes would otherwise cost 200 spans a
    // frame, all with identical styling.
    let mut run = String::new();
    let mut run_style: Option<Style> = None;
    let flush = |run: &mut String, style: &mut Option<Style>, spans: &mut Vec<Span>| {
        if !run.is_empty() {
            spans.push(Span::styled(std::mem::take(run), style.unwrap_or_default()));
        }
    };

    for offset in 0..len {
        let index = start + offset;
        let Some(residue) = chain.residues.get(index) else {
            break;
        };
        // Group separator, styled like the gutter so it never reads as a gap
        // in the sequence itself.
        if offset > 0 && offset % GROUP == 0 {
            flush(&mut run, &mut run_style, &mut spans);
            run_style = None;
            spans.push(Span::styled(" ", Style::default().fg(Color::DarkGray)));
        }

        let style = residue_style(app, chain, chain_index, index);
        if run_style != Some(style) {
            flush(&mut run, &mut run_style, &mut spans);
            run_style = Some(style);
        }
        run.push(one_letter(&residue.name, chain.molecule_type));
    }
    flush(&mut run, &mut run_style, &mut spans);

    Line::from(spans)
}

fn residue_style(app: &App, chain: &Chain, chain_index: usize, index: usize) -> Style {
    let is_cursor = app.show_sequence && app.seq_cursor == (chain_index, index);
    let is_selected = app.selection.contains(chain_index, index);

    if is_cursor {
        return Style::default()
            .bg(cursor_color())
            .fg(Color::Black)
            .add_modifier(Modifier::BOLD);
    }
    if is_selected {
        return Style::default()
            .bg(selection_color())
            .fg(Color::Black)
            .add_modifier(Modifier::BOLD);
    }
    // Otherwise the letter carries the same color the residue has in 3D, so
    // the panel and the structure read as one picture.
    let color = chain
        .residues
        .get(index)
        .map(|residue| app.color_scheme.residue_color(residue, chain))
        .unwrap_or(Color::Gray);
    Style::default().fg(color)
}

/// Screen column of a residue within the row, used by tests and any future
/// mouse support.
#[allow(dead_code)]
pub fn residue_column(column: usize) -> u16 {
    GUTTER + column_offset(column) as u16
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::app::{AppConfig, RenderMode, VizMode};
    use crate::model::protein::{Atom, MoleculeType, Protein, Residue, SecondaryStructure};
    use crate::model::selection::ResidueColorOverrides;
    use ratatui::Terminal;
    use ratatui::backend::TestBackend;
    use ratatui_image::picker::Picker;

    fn protein() -> Protein {
        let residue = |i: usize, name: &str| Residue {
            name: name.to_string(),
            seq_num: i as i32 + 1,
            insertion_code: None,
            atoms: vec![Atom {
                name: "CA".to_string(),
                element: "C".to_string(),
                x: i as f64,
                y: 0.0,
                z: 0.0,
                b_factor: 10.0,
                is_backbone: true,
                is_hetero: false,
            }],
            secondary_structure: SecondaryStructure::Coil,
        };
        Protein {
            name: "panel".to_string(),
            chains: vec![
                Chain {
                    id: "L1".to_string(),
                    molecule_type: MoleculeType::RNA,
                    residues: (0..25).map(|i| residue(i, "G")).collect(),
                },
                Chain {
                    id: "Lb".to_string(),
                    molecule_type: MoleculeType::Protein,
                    residues: (0..12).map(|i| residue(i, "TRP")).collect(),
                },
            ],
            ligands: Vec::new(),
        }
    }

    fn app() -> App {
        let mut app = App::new(
            protein(),
            AppConfig {
                render_mode: RenderMode::Braille,
                viz_mode: VizMode::Backbone,
                user_explicit_mode: true,
                color_override: None,
                residue_colors: ResidueColorOverrides::default(),
            },
            80,
            40,
            Picker::halfblocks(),
        );
        app.show_sequence = true;
        app.set_sequence_viewport(80, 10);
        app
    }

    fn draw(app: &App) -> String {
        let mut terminal = Terminal::new(TestBackend::new(80, 10)).unwrap();
        terminal
            .draw(|frame| render_sequence_panel(frame, Rect::new(0, 0, 80, 10), app))
            .unwrap();
        let buffer = terminal.backend().buffer().clone();
        (0..10)
            .map(|y| {
                (0..80)
                    .map(|x| buffer[(x, y)].symbol().to_string())
                    .collect::<String>()
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn renders_every_chain_with_one_letter_codes() {
        let text = draw(&app());
        assert!(text.contains("L1"), "chain header missing:\n{text}");
        assert!(text.contains("GGGGGGGGGG"), "RNA sequence missing:\n{text}");
        assert!(text.contains("Lb"), "second chain missing:\n{text}");
        assert!(
            text.contains("WWWWWWWWWW"),
            "protein sequence missing:\n{text}"
        );
    }

    #[test]
    fn selection_summary_reports_picked_residues() {
        let mut app = app();
        app.seq_move_horizontal(3, false);
        app.seq_toggle_selection();
        app.seq_move_horizontal(2, true);
        let text = draw(&app);
        assert!(text.contains("3 res in 1 chain"), "{text}");
        assert!(text.contains("L1:4-6"), "{text}");
    }
}

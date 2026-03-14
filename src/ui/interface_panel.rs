use ratatui::layout::Rect;
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph, Wrap};
use ratatui::Frame;

/// Width of the interface sidebar in columns.
pub const SIDEBAR_WIDTH: u16 = 32;

/// Render the interface analysis as a left-edge sidebar (non-occluding).
pub fn render_interface_panel(
    frame: &mut Frame,
    area: Rect,
    summary_lines: &[String],
    focus_chain: usize,
    chain_names: &[String],
    show_interactions: bool,
    interaction_counts: (usize, usize, usize, usize),
) {
    let mut lines: Vec<Line> = Vec::new();

    // Chain roles
    let focus_name = chain_names
        .get(focus_chain)
        .cloned()
        .unwrap_or_else(|| "?".to_string());

    let other_names: Vec<&str> = chain_names
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != focus_chain)
        .map(|(_, name)| name.as_str())
        .collect();
    let other_label = if other_names.is_empty() {
        "-".to_string()
    } else {
        other_names.join(", ")
    };

    lines.push(Line::from(vec![
        Span::styled(" Ab: ", Style::default().fg(Color::DarkGray)),
        Span::styled(
            format!("[{}]", focus_name),
            Style::default()
                .fg(Color::Green)
                .add_modifier(Modifier::BOLD),
        ),
    ]));

    lines.push(Line::from(vec![
        Span::styled(" Ag: ", Style::default().fg(Color::DarkGray)),
        Span::styled(
            format!("[{}]", other_label),
            Style::default()
                .fg(Color::Rgb(255, 165, 0))
                .add_modifier(Modifier::BOLD),
        ),
    ]));

    lines.push(Line::from(Span::styled(
        " [/]: cycle chain",
        Style::default().fg(Color::DarkGray),
    )));

    lines.push(Line::from(""));

    // Contact summary
    for text in summary_lines {
        lines.push(Line::from(Span::styled(
            format!(" {}", text),
            Style::default().fg(Color::White),
        )));
    }

    lines.push(Line::from(""));

    // Visual key
    lines.push(Line::from(vec![
        Span::styled(" ", Style::default()),
        Span::styled("\u{2588}", Style::default().fg(Color::Rgb(0, 255, 100))),
        Span::styled(" Ab interface", Style::default().fg(Color::DarkGray)),
    ]));
    lines.push(Line::from(vec![
        Span::styled(" ", Style::default()),
        Span::styled("\u{2588}", Style::default().fg(Color::Rgb(40, 100, 60))),
        Span::styled(" Ab other", Style::default().fg(Color::DarkGray)),
    ]));
    lines.push(Line::from(vec![
        Span::styled(" ", Style::default()),
        Span::styled("\u{2588}", Style::default().fg(Color::Rgb(255, 165, 0))),
        Span::styled(" Ag interface", Style::default().fg(Color::DarkGray)),
    ]));
    lines.push(Line::from(vec![
        Span::styled(" ", Style::default()),
        Span::styled("\u{2588}", Style::default().fg(Color::Rgb(100, 80, 60))),
        Span::styled(" Ag other", Style::default().fg(Color::DarkGray)),
    ]));

    if summary_lines.iter().any(|l| l.contains("Ligand") || l.contains("Ion")) {
        lines.push(Line::from(""));
        lines.push(Line::from(vec![
            Span::styled(" ", Style::default()),
            Span::styled("\u{2588}", Style::default().fg(Color::Rgb(255, 0, 255))),
            Span::styled(" Ligand", Style::default().fg(Color::DarkGray)),
        ]));
        lines.push(Line::from(vec![
            Span::styled(" ", Style::default()),
            Span::styled("\u{2588}", Style::default().fg(Color::Rgb(0, 255, 255))),
            Span::styled(" Ion", Style::default().fg(Color::DarkGray)),
        ]));
    }

    // Interaction section
    lines.push(Line::from(""));
    if show_interactions {
        let (hbonds, salt_bridges, hydrophobic, other) = interaction_counts;
        lines.push(Line::from(Span::styled(
            " Interactions",
            Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD),
        )));
        lines.push(Line::from(vec![
            Span::styled(" ", Style::default()),
            Span::styled("\u{2588}", Style::default().fg(Color::Rgb(0, 220, 255))),
            Span::styled(format!(" H-bond: {}", hbonds), Style::default().fg(Color::White)),
        ]));
        lines.push(Line::from(vec![
            Span::styled(" ", Style::default()),
            Span::styled("\u{2588}", Style::default().fg(Color::Rgb(255, 80, 80))),
            Span::styled(format!(" Salt bridge: {}", salt_bridges), Style::default().fg(Color::White)),
        ]));
        lines.push(Line::from(vec![
            Span::styled(" ", Style::default()),
            Span::styled("\u{2588}", Style::default().fg(Color::Rgb(220, 200, 60))),
            Span::styled(format!(" Hydrophobic: {}", hydrophobic), Style::default().fg(Color::White)),
        ]));
        lines.push(Line::from(vec![
            Span::styled(" ", Style::default()),
            Span::styled("\u{2588}", Style::default().fg(Color::Rgb(160, 160, 160))),
            Span::styled(format!(" Other: {}", other), Style::default().fg(Color::White)),
        ]));
    } else {
        lines.push(Line::from(Span::styled(
            " [Shift+I]: show interactions",
            Style::default().fg(Color::DarkGray),
        )));
    }

    let panel = Paragraph::new(lines)
        .block(
            Block::default()
                .borders(Borders::RIGHT)
                .border_style(Style::default().fg(Color::DarkGray))
                .title(" Interface ")
                .title_style(Style::default().fg(Color::Cyan)),
        )
        .wrap(Wrap { trim: false });

    frame.render_widget(panel, area);
}

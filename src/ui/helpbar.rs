use ratatui::Frame;
use ratatui::layout::Rect;
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;

use crate::app::App;

/// Render the keybinding hints bar at the bottom.
///
/// The hints change with the sequence panel: while it is open the arrow keys
/// drive the cursor rather than the camera, and saying so here is cheaper than
/// making the user open the help overlay to find out.
pub fn render_helpbar(frame: &mut Frame, area: Rect, app: &App) {
    if app.show_sequence {
        frame.render_widget(Paragraph::new(sequence_hints()), area);
        return;
    }
    let help = Paragraph::new(Line::from(vec![
        Span::styled("╰── ", Style::default().fg(Color::DarkGray)),
        Span::styled(
            "h/l",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": rotY  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "j/k",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": rotX  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "+/-",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": zoom  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "wasd",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": pan  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "c",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": color  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "v",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": mode  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "f",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": interface  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "I",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": interactions  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "?",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": help  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "g",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": ligands  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "S",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": sequence  ", Style::default().fg(Color::Gray)),
        Span::styled(
            "q",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(": quit ", Style::default().fg(Color::Gray)),
        Span::styled("──╯", Style::default().fg(Color::DarkGray)),
    ]));
    frame.render_widget(help, area);
}

/// Hints shown while the sequence panel has the arrow keys.
fn sequence_hints() -> Line<'static> {
    let key = |text: &'static str| {
        Span::styled(
            text,
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        )
    };
    let label = |text: &'static str| Span::styled(text, Style::default().fg(Color::Gray));

    Line::from(vec![
        Span::styled("╰── ", Style::default().fg(Color::DarkGray)),
        key("←→↑↓"),
        label(": cursor  "),
        key("shift+←→"),
        label(": range  "),
        key("↵"),
        label(": pick  "),
        key("A"),
        label(": chain  "),
        key("x"),
        label(": clear  "),
        key("b"),
        label(": ball&stick  "),
        key("z"),
        label(": centre  "),
        key("[ ]"),
        label(": chain  "),
        key("hjkl"),
        label(": rotate  "),
        key("<>"),
        label(": size  "),
        key("S/esc"),
        label(": close "),
        Span::styled("──╯", Style::default().fg(Color::DarkGray)),
    ])
}

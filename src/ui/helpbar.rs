use ratatui::layout::Rect;
use ratatui::style::{Color, Style, Modifier};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

/// Render the keybinding hints bar at the bottom
pub fn render_helpbar(frame: &mut Frame, area: Rect) {
    let help = Paragraph::new(Line::from(vec![
        Span::styled("╰── ", Style::default().fg(Color::DarkGray)),
        Span::styled("h/l", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)),
        Span::styled(": rotY  ", Style::default().fg(Color::Gray)),
        Span::styled("j/k", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)),
        Span::styled(": rotX  ", Style::default().fg(Color::Gray)),
        Span::styled("+/-", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)),
        Span::styled(": zoom  ", Style::default().fg(Color::Gray)),
        Span::styled("wasd", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)),
        Span::styled(": pan  ", Style::default().fg(Color::Gray)),
        Span::styled("c", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)),
        Span::styled(": color  ", Style::default().fg(Color::Gray)),
        Span::styled("v", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)),
        Span::styled(": mode  ", Style::default().fg(Color::Gray)),
        Span::styled("?", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)),
        Span::styled(": help  ", Style::default().fg(Color::Gray)),
        Span::styled("q", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)),
        Span::styled(": quit ", Style::default().fg(Color::Gray)),
        Span::styled("──╯", Style::default().fg(Color::DarkGray)),
    ]));
    frame.render_widget(help, area);
}

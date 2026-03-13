use ratatui::layout::Rect;
use ratatui::style::{Color, Style, Modifier};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Clear, Paragraph, Wrap};
use ratatui::Frame;

/// Render a centered help overlay
pub fn render_help_overlay(frame: &mut Frame, area: Rect) {
    // Center the popup (60x20)
    let popup_width = 60u16.min(area.width - 4);
    let popup_height = 21u16.min(area.height - 4);
    let x = (area.width - popup_width) / 2;
    let y = (area.height - popup_height) / 2;
    let popup_area = Rect::new(x, y, popup_width, popup_height);

    frame.render_widget(Clear, popup_area);

    let help_text = vec![
        Line::from(Span::styled("  Keybindings", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD))),
        Line::from(""),
        Line::from(vec![Span::styled("  h / l      ", Style::default().fg(Color::Yellow)), Span::raw("Rotate Y-axis")]),
        Line::from(vec![Span::styled("  j / k      ", Style::default().fg(Color::Yellow)), Span::raw("Rotate X-axis")]),
        Line::from(vec![Span::styled("  u / i      ", Style::default().fg(Color::Yellow)), Span::raw("Rotate Z-axis (roll)")]),
        Line::from(vec![Span::styled("  + / -      ", Style::default().fg(Color::Yellow)), Span::raw("Zoom in / out")]),
        Line::from(vec![Span::styled("  w/a/s/d    ", Style::default().fg(Color::Yellow)), Span::raw("Pan up/left/down/right")]),
        Line::from(vec![Span::styled("  r          ", Style::default().fg(Color::Yellow)), Span::raw("Reset view")]),
        Line::from(vec![Span::styled("  c          ", Style::default().fg(Color::Yellow)), Span::raw("Cycle color scheme")]),
        Line::from(vec![Span::styled("  v          ", Style::default().fg(Color::Yellow)), Span::raw("Cycle viz mode")]),
        Line::from(vec![Span::styled("  m          ", Style::default().fg(Color::Yellow)), Span::raw("Toggle Braille / HD")]),
        Line::from(vec![Span::styled("  M          ", Style::default().fg(Color::Yellow)), Span::raw("Toggle HD / FullHD (Sixel/Kitty)")]),
        Line::from(vec![Span::styled("  [ / ]      ", Style::default().fg(Color::Yellow)), Span::raw("Prev / next chain")]),
        Line::from(vec![Span::styled("  g          ", Style::default().fg(Color::Yellow)), Span::raw("Toggle ligand visibility")]),
        Line::from(vec![Span::styled("  Space      ", Style::default().fg(Color::Yellow)), Span::raw("Toggle auto-rotation")]),
        Line::from(vec![Span::styled("  ?          ", Style::default().fg(Color::Yellow)), Span::raw("Toggle this help")]),
        Line::from(vec![Span::styled("  q          ", Style::default().fg(Color::Yellow)), Span::raw("Quit")]),
        Line::from(""),
        Line::from(Span::styled("  Press ? or Esc to close", Style::default().fg(Color::DarkGray))),
    ];

    let help = Paragraph::new(help_text)
        .block(Block::default()
            .borders(Borders::ALL)
            .border_style(Style::default().fg(Color::Cyan))
            .title(" Help "))
        .wrap(Wrap { trim: false });

    frame.render_widget(help, popup_area);
}

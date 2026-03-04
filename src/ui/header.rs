use ratatui::layout::Rect;
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

/// Render the retro-styled header with protein name
pub fn render_header(frame: &mut Frame, area: Rect, protein_name: &str) {
    let title_text = format!(" ProteinView ─── {} ", protein_name);
    let fill = " ─".repeat((area.width as usize).saturating_sub(title_text.len() + 4) / 2);
    let header = Paragraph::new(Line::from(vec![
        Span::styled("╭─── ", Style::default().fg(Color::DarkGray)),
        Span::styled(
            "ProteinView",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(" ─── ", Style::default().fg(Color::DarkGray)),
        Span::styled(protein_name, Style::default().fg(Color::Yellow)),
        Span::styled(fill, Style::default().fg(Color::DarkGray)),
        Span::styled("╮", Style::default().fg(Color::DarkGray)),
    ]));
    frame.render_widget(header, area);
}

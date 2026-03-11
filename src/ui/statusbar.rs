use ratatui::Frame;
use ratatui::layout::Rect;
use ratatui::style::{Color, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;

use crate::app::App;

/// Render the status bar showing current mode and info
pub fn render_statusbar(frame: &mut Frame, area: Rect, app: &App) {
    let chain_info = if let Some(chain) = app.protein.chains.get(app.current_chain) {
        format!("Chain {} ", chain.id)
    } else {
        "No chains ".to_string()
    };

    let res_count = app.protein.residue_count();
    let render_mode = if app.hd_mode { "HD" } else { "Braille" };

    let status = Paragraph::new(Line::from(vec![
        Span::styled("├", Style::default().fg(Color::DarkGray)),
        Span::styled(
            "─".repeat(area.width as usize - 2),
            Style::default().fg(Color::DarkGray),
        ),
        Span::styled("┤", Style::default().fg(Color::DarkGray)),
    ]));
    frame.render_widget(status, area);

    // Render the actual status info on the next line if area has height > 1
    if area.height > 1 {
        let info_area = Rect::new(area.x, area.y + 1, area.width, 1);
        let info = Paragraph::new(Line::from(vec![
            Span::styled("│ ", Style::default().fg(Color::DarkGray)),
            Span::styled(&chain_info, Style::default().fg(Color::Cyan)),
            Span::styled("│ ", Style::default().fg(Color::DarkGray)),
            Span::styled(
                format!("{} res ", res_count),
                Style::default().fg(Color::White),
            ),
            Span::styled("│ ", Style::default().fg(Color::DarkGray)),
            Span::styled(app.viz_mode.name(), Style::default().fg(Color::Green)),
            Span::styled(" │ ", Style::default().fg(Color::DarkGray)),
            Span::styled(
                app.color_scheme.scheme_type.name(),
                Style::default().fg(Color::Yellow),
            ),
            Span::styled(" │ ", Style::default().fg(Color::DarkGray)),
            Span::styled(render_mode, Style::default().fg(Color::Magenta)),
            Span::raw(" "),
        ]));
        frame.render_widget(info, info_area);
    }
}

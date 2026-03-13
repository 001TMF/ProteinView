use ratatui::layout::Rect;
use ratatui::style::{Color, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::app::{App, ConnectionType, RenderMode, VizMode};

/// Render the status bar showing current mode and info
pub fn render_statusbar(frame: &mut Frame, area: Rect, app: &App) {
    let chain_info = if let Some(chain) = app.protein.chains.get(app.current_chain) {
        format!("Chain {} ", chain.id)
    } else {
        "No chains ".to_string()
    };

    let res_count = app.protein.residue_count();
    let render_mode_name = app.render_mode.name();

    let border_fill = (area.width as usize).saturating_sub(2);
    let status = Paragraph::new(Line::from(vec![
        Span::styled("\u{251c}", Style::default().fg(Color::DarkGray)),
        Span::styled("\u{2500}".repeat(border_fill), Style::default().fg(Color::DarkGray)),
        Span::styled("\u{2524}", Style::default().fg(Color::DarkGray)),
    ]));
    frame.render_widget(status, area);

    // Render the actual status info on the next line if area has height > 1
    if area.height > 1 {
        let info_area = Rect::new(area.x, area.y + 1, area.width, 1);
        // Build spans dynamically
        let mut spans = vec![
            Span::styled("\u{2502} ", Style::default().fg(Color::DarkGray)),
            Span::styled(&chain_info, Style::default().fg(Color::Cyan)),
            Span::styled("\u{2502} ", Style::default().fg(Color::DarkGray)),
            Span::styled(format!("{} res ", res_count), Style::default().fg(Color::White)),
        ];

        if app.protein.ligand_count() > 0 {
            spans.push(Span::styled("\u{2502} ", Style::default().fg(Color::DarkGray)));
            spans.push(Span::styled(
                format!("{} ligands ", app.protein.ligand_count()),
                Style::default().fg(Color::Rgb(255, 0, 255)),
            ));
        }

        // In Braille mode, Cartoon is silently degraded to Backbone rendering
        // (the basic braille renderer has no triangle mesh support).  Show this
        // honestly so the user isn't confused.
        let viz_label = if app.render_mode == RenderMode::Braille && app.viz_mode == VizMode::Cartoon {
            "Backbone*"
        } else {
            app.viz_mode.name()
        };

        spans.extend_from_slice(&[
            Span::styled("\u{2502} ", Style::default().fg(Color::DarkGray)),
            Span::styled(viz_label, Style::default().fg(Color::Green)),
            Span::styled(" \u{2502} ", Style::default().fg(Color::DarkGray)),
            Span::styled(app.color_scheme.scheme_type.name(), Style::default().fg(Color::Yellow)),
            Span::styled(" \u{2502} ", Style::default().fg(Color::DarkGray)),
            Span::styled(render_mode_name, Style::default().fg(Color::Magenta)),
            Span::raw(" "),
        ]);

        // Show SSH indicator
        if app.connection_type == ConnectionType::Ssh {
            spans.push(Span::styled("\u{2502} ", Style::default().fg(Color::DarkGray)));
            spans.push(Span::styled("SSH", Style::default().fg(Color::Rgb(255, 165, 0))));
            spans.push(Span::raw(" "));
        }

        // Show SSH HD warning
        if app.ssh_hd_warning {
            spans.push(Span::styled("\u{2502} ", Style::default().fg(Color::DarkGray)));
            spans.push(Span::styled(
                "\u{26a0} Full HD over SSH may be slow",
                Style::default().fg(Color::Rgb(255, 200, 0)),
            ));
            spans.push(Span::raw(" "));
        }

        let info = Paragraph::new(Line::from(spans));
        frame.render_widget(info, info_area);
    }
}

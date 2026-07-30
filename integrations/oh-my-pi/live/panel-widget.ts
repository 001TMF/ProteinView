import type { Component, ImageDimensions, ImageOptions, ImageTheme, TUI } from "@oh-my-pi/pi-tui";
import { Image, truncateToWidth, visibleWidth } from "@oh-my-pi/pi-tui";
import type { PanelSnapshot } from "./panel-client.ts";

const MIN_EXPANDED_ROWS = 18;
const MIN_EXPANDED_COLUMNS = 44;
const EDITOR_RESERVE_ROWS = 10;
const PANEL_HEIGHT_FRACTION = 0.45;
const PANEL_CHROME_ROWS = 2;
const IMAGE_KEY = "proteinview-live-panel";
// biome-ignore lint/complexity/useRegexLiterals: Constructor strings avoid literal control escapes rejected by noControlCharactersInRegex.
const ANSI_ESCAPE = new RegExp("\\x1b(?:\\[[0-?]*[ -/]*[@-~]|\\][^\\x07]*(?:\\x07|\\x1b\\\\))", "g");
// biome-ignore lint/complexity/useRegexLiterals: Constructor strings avoid literal control escapes rejected by noControlCharactersInRegex.
const DISPLAY_CONTROL_CHARACTERS = new RegExp("[\\u0000-\\u001f\\u007f-\\u009f]", "g");

export interface ProteinViewPanelTheme {
	fg(color: "accent" | "dim" | "muted" | "error", text: string): string;
	bold(text: string): string;
}

export interface ProteinViewPanelWidgetOptions {
	client: ProteinViewPanelWidgetClient;
	identifier: string;
	tui: TUI;
	theme: ProteinViewPanelTheme;
	imageFactory?: ProteinViewPanelImageFactory;
}

export interface ProteinViewPanelWidgetClient {
	readonly snapshot: PanelSnapshot | undefined;
	subscribe(listener: (snapshot: PanelSnapshot) => void): () => void;
}

export interface ProteinViewPanelImage extends Component {
	setData(base64Data: string, mimeType: string, dimensions?: ImageDimensions): void;
}

export type ProteinViewPanelImageFactory = (
	base64Data: string,
	mimeType: string,
	theme: ImageTheme,
	options: ImageOptions,
	dimensions: ImageDimensions,
) => ProteinViewPanelImage;

function fitLine(line: string, width: number): string {
	if (width <= 0) return "";
	const truncated = truncateToWidth(line, width);
	const remaining = Math.max(0, width - visibleWidth(truncated));
	return truncated + " ".repeat(remaining);
}

function maxPanelRows(terminalRows: number): number {
	return Math.max(0, Math.min(Math.floor(terminalRows * PANEL_HEIGHT_FRACTION), terminalRows - EDITOR_RESERVE_ROWS));
}

function safeLabel(value: string, fallback: string): string {
	const cleaned = value.replace(ANSI_ESCAPE, "").replace(DISPLAY_CONTROL_CHARACTERS, " ").trim().slice(0, 120);
	return cleaned || fallback;
}

export class ProteinViewPanelWidget implements Component {
	#client: ProteinViewPanelWidgetClient;
	#identifier: string;
	#tui: TUI;
	#theme: ProteinViewPanelTheme;
	#imageFactory: ProteinViewPanelImageFactory;
	#snapshot: PanelSnapshot;
	#image: ProteinViewPanelImage | undefined;
	#imageHeight = 0;
	#focused = false;
	#error: string | undefined;
	#unsubscribe: () => void;
	#disposed = false;
	#cachedWidth = 0;
	#cachedTerminalRows = 0;
	#cachedRevision = 0;
	#cachedFocused = false;
	#cachedError: string | undefined;
	#cachedLines: readonly string[] | undefined;

	constructor(options: ProteinViewPanelWidgetOptions) {
		const snapshot = options.client.snapshot;
		if (snapshot === undefined) {
			throw new Error("ProteinView panel widget requires a ready frame");
		}
		this.#client = options.client;
		this.#identifier = safeLabel(options.identifier, "structure");
		this.#tui = options.tui;
		this.#theme = options.theme;
		this.#imageFactory =
			options.imageFactory ??
			((base64Data, mimeType, theme, imageOptions, dimensions) =>
				new Image(base64Data, mimeType, theme, imageOptions, dimensions));
		this.#snapshot = snapshot;
		this.#unsubscribe = this.#client.subscribe(next => {
			this.#snapshot = next;
			if (this.#image !== undefined) {
				this.#image.setData(next.png.toString("base64"), "image/png", {
					widthPx: next.frame.width,
					heightPx: next.frame.height,
				});
			}
			this.invalidate();
			this.#tui.requestRender();
		});
	}

	setFocused(focused: boolean): void {
		if (this.#focused === focused) return;
		this.#focused = focused;
		this.invalidate();
		this.#tui.requestRender();
	}

	setError(error: string | undefined): void {
		if (this.#error === error) return;
		this.#error = error;
		this.invalidate();
		this.#tui.requestRender();
	}

	invalidate(): void {
		this.#cachedLines = undefined;
		this.#image?.invalidate?.();
	}

	render(width: number): readonly string[] {
		const terminalRows = Math.max(1, this.#tui.terminal.rows);
		if (
			this.#cachedLines !== undefined &&
			this.#cachedWidth === width &&
			this.#cachedTerminalRows === terminalRows &&
			this.#cachedRevision === this.#snapshot.revision &&
			this.#cachedFocused === this.#focused &&
			this.#cachedError === this.#error
		) {
			return this.#cachedLines;
		}

		const budget = maxPanelRows(terminalRows);
		const collapsed =
			terminalRows < MIN_EXPANDED_ROWS || width < MIN_EXPANDED_COLUMNS || budget < PANEL_CHROME_ROWS + 2;
		const header = fitLine(this.#headerText(collapsed), width);
		const hint = fitLine(this.#hintText(collapsed), width);
		let lines: readonly string[];

		if (collapsed) {
			lines = [header, hint];
		} else {
			const nextImageHeight = Math.max(1, budget - PANEL_CHROME_ROWS);
			if (this.#image === undefined || this.#imageHeight !== nextImageHeight) {
				this.#imageHeight = nextImageHeight;
				this.#image = this.#imageFactory(
					this.#snapshot.png.toString("base64"),
					"image/png",
					{ fallbackColor: text => this.#theme.fg("muted", text) },
					{
						maxHeightCells: nextImageHeight,
						filename: `${this.#identifier}.png`,
						budget: this.#tui.imageBudget,
						imageKey: IMAGE_KEY,
					},
					{
						widthPx: this.#snapshot.frame.width,
						heightPx: this.#snapshot.frame.height,
					},
				);
			}
			const imageLines = this.#image.render(width);
			lines = [header, ...imageLines.slice(0, nextImageHeight), hint];
		}

		this.#cachedWidth = width;
		this.#cachedTerminalRows = terminalRows;
		this.#cachedRevision = this.#snapshot.revision;
		this.#cachedFocused = this.#focused;
		this.#cachedError = this.#error;
		this.#cachedLines = lines;
		return lines;
	}

	dispose(): void {
		if (this.#disposed) return;
		this.#disposed = true;
		this.#unsubscribe();
	}

	#headerText(collapsed: boolean): string {
		const state = this.#snapshot.state;
		const chain = safeLabel(state.presentation.current_chain_id ?? "—", "—");
		const focus = this.#focused ? " · focused" : "";
		const collapsedLabel = collapsed ? " · compact" : "";
		const label =
			` ProteinView · ${this.#identifier} · FullHD source 1920×1080 · ` +
			`${state.presentation.viz_mode}/${state.presentation.effective_color} · chain ${chain}` +
			`${focus}${collapsedLabel}`;
		if (this.#error !== undefined) {
			return this.#theme.fg("error", `${label} · ${this.#error}`);
		}
		return this.#focused ? this.#theme.fg("accent", this.#theme.bold(label)) : this.#theme.fg("muted", label);
	}

	#hintText(collapsed: boolean): string {
		if (collapsed) {
			return this.#theme.fg("dim", " F8 focus/expand when terminal is larger · q close");
		}
		return this.#theme.fg(
			"dim",
			" F8 focus · Esc/Tab chat · ←↑↓→ rotate · +/- zoom · r reset · f fit · c color · v view · [/] chain · i interface · n contacts · l ligands · a spin · q close",
		);
	}
}

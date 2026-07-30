import type { ExtensionAPI, ExtensionContext } from "@oh-my-pi/pi-coding-agent";
import { matchesKey, replaceTabs } from "@oh-my-pi/pi-tui";
import {
	type ProteinViewLiveOpenRequest,
	type ProteinViewLiveResidueColor,
	registerProteinViewLiveOpenHandler,
} from "../live/open-bridge.ts";
import {
	type PanelClientOptions,
	type PanelCommand,
	type PanelCommandResult,
	type PanelState,
	ProteinViewPanelClient,
} from "../live/panel-client.ts";
import { ProteinViewPanelWidget } from "../live/panel-widget.ts";

const WIDGET_KEY = "proteinview-live";
const STATUS_KEY = "proteinview-live-focus";
const ROTATION_STEP = 0.12;
const AUTO_ROTATE_INTERVAL_MS = 100;
const MAX_RESIDUE_COLORS = 256;
// biome-ignore lint/complexity/useRegexLiterals: Constructor strings avoid literal control escapes rejected by noControlCharactersInRegex.
const VISIBLE_CONTROL_CHARACTERS = new RegExp("[\\u0000-\\u001f\\u007f-\\u009f]", "g");

export interface ProteinViewPanelClientLike {
	readonly snapshot: ProteinViewPanelClient["snapshot"];
	readonly state: PanelState | undefined;
	readonly closed: boolean;
	readonly failure: Error | undefined;
	subscribe(listener: Parameters<ProteinViewPanelClient["subscribe"]>[0]): () => void;
	subscribeFailure(listener: Parameters<ProteinViewPanelClient["subscribeFailure"]>[0]): () => void;
	command(command: PanelCommand): Promise<PanelCommandResult>;
	close(): Promise<void>;
	abort(error?: Error): Promise<void>;
}

export interface ProteinViewExtensionDependencies {
	startClient(options: PanelClientOptions): Promise<ProteinViewPanelClientLike>;
}

export type ProteinViewControlDetails = { closed: true } | PanelCommandResult;

const DEFAULT_DEPENDENCIES: ProteinViewExtensionDependencies = {
	startClient: options => ProteinViewPanelClient.start(options),
};

function cancelled(signal: AbortSignal | undefined): Error | undefined {
	if (!signal?.aborted) return undefined;
	return signal.reason instanceof Error ? signal.reason : new Error("ProteinView control cancelled");
}

function visibleError(error: unknown): string {
	const message = error instanceof Error ? error.message : String(error);
	return replaceTabs(message).replace(VISIBLE_CONTROL_CHARACTERS, " ").trim().slice(0, 240);
}

function requireField<T>(value: T | undefined, label: string): T {
	if (value === undefined) throw new Error(`${label} is required for this action`);
	return value;
}

export class ProteinViewPanelHost {
	#dependencies: ProteinViewExtensionDependencies;
	#client: ProteinViewPanelClientLike | undefined;
	#context: ExtensionContext | undefined;
	#widget: ProteinViewPanelWidget | undefined;
	#inputUnsubscribe: (() => void) | undefined;
	#failureUnsubscribe: (() => void) | undefined;
	#focused = false;
	#autoRotateInFlight = false;
	#lifecycle: Promise<void> = Promise.resolve();

	constructor(dependencies: ProteinViewExtensionDependencies = DEFAULT_DEPENDENCIES) {
		this.#dependencies = dependencies;
	}

	get active(): boolean {
		return this.#client !== undefined && !this.#client.closed;
	}

	get state(): PanelState | undefined {
		return this.#client?.state;
	}

	open(request: ProteinViewLiveOpenRequest, context: ExtensionContext): Promise<void> {
		const next = this.#lifecycle.then(() => this.#open(request, context));
		this.#lifecycle = next.then(
			() => undefined,
			() => undefined,
		);
		return next;
	}

	close(): Promise<void> {
		const next = this.#lifecycle.then(() => this.#close());
		this.#lifecycle = next.then(
			() => undefined,
			() => undefined,
		);
		return next;
	}

	toggleFocus(context?: ExtensionContext): void {
		const uiContext = context ?? this.#context;
		if (!this.active || uiContext === undefined) {
			uiContext?.ui.notify("No live ProteinView panel is open", "info");
			return;
		}
		this.#setFocused(!this.#focused);
	}

	handleTerminalInput(data: string): { consume?: boolean } | undefined {
		if (matchesKey(data, "f8")) {
			this.toggleFocus();
			return { consume: true };
		}
		if (!this.#focused) return undefined;
		if (matchesKey(data, "escape") || matchesKey(data, "tab")) {
			this.#setFocused(false);
			return { consume: true };
		}
		if (matchesKey(data, "q")) {
			void this.close();
			return { consume: true };
		}

		let command: PanelCommand | undefined;
		if (matchesKey(data, "left")) {
			command = { command: "rotate", axis: "y", delta: -ROTATION_STEP };
		} else if (matchesKey(data, "right")) {
			command = { command: "rotate", axis: "y", delta: ROTATION_STEP };
		} else if (matchesKey(data, "up")) {
			command = { command: "rotate", axis: "x", delta: -ROTATION_STEP };
		} else if (matchesKey(data, "down")) {
			command = { command: "rotate", axis: "x", delta: ROTATION_STEP };
		} else if (matchesKey(data, "+") || matchesKey(data, "=")) {
			command = { command: "zoom", direction: "in" };
		} else if (matchesKey(data, "-")) {
			command = { command: "zoom", direction: "out" };
		} else if (matchesKey(data, "r")) {
			command = { command: "reset" };
		} else if (matchesKey(data, "f")) {
			command = { command: "fit" };
		} else if (matchesKey(data, "c")) {
			command = { command: "cycle_color" };
		} else if (matchesKey(data, "v")) {
			command = { command: "cycle_viz" };
		} else if (matchesKey(data, "[")) {
			command = { command: "select_chain", direction: "prev" };
		} else if (matchesKey(data, "]")) {
			command = { command: "select_chain", direction: "next" };
		} else if (matchesKey(data, "i")) {
			command = {
				command: "set_interface",
				enabled: !(this.state?.presentation.interface ?? false),
			};
		} else if (matchesKey(data, "n")) {
			command = {
				command: "set_interactions",
				enabled: !(this.state?.presentation.interactions ?? false),
			};
		} else if (matchesKey(data, "l")) {
			command = {
				command: "set_ligands",
				enabled: !(this.state?.presentation.ligands ?? true),
			};
		} else if (matchesKey(data, "a")) {
			command = {
				command: "set_auto_rotate",
				enabled: !(this.state?.camera.auto_rotate ?? false),
			};
		}
		if (command === undefined) return { consume: true };
		this.#dispatch(command);
		return { consume: true };
	}

	async control(command: PanelCommand, signal?: AbortSignal): Promise<PanelCommandResult> {
		const cancellation = cancelled(signal);
		if (cancellation !== undefined) throw cancellation;
		const client = this.#client;
		if (client === undefined || client.closed) {
			throw new Error("No live ProteinView panel is open");
		}
		const onAbort = (): void => {
			void client.abort(cancelled(signal) ?? new Error("ProteinView control cancelled"));
		};
		signal?.addEventListener("abort", onAbort, { once: true });
		try {
			const result = await client.command(command);
			const afterCommand = cancelled(signal);
			if (afterCommand !== undefined) {
				await client.abort(afterCommand);
				throw afterCommand;
			}
			this.#widget?.setError(undefined);
			return result;
		} finally {
			signal?.removeEventListener("abort", onAbort);
		}
	}

	autoRotateTick(): void {
		const client = this.#client;
		if (client === undefined || client.closed || !client.state?.camera.auto_rotate || this.#autoRotateInFlight) {
			return;
		}
		this.#autoRotateInFlight = true;
		void client
			.command({ command: "render" })
			.catch(error => {
				this.#widget?.setError(visibleError(error));
			})
			.finally(() => {
				this.#autoRotateInFlight = false;
			});
	}

	async #open(request: ProteinViewLiveOpenRequest, context: ExtensionContext): Promise<void> {
		await this.#close();
		const client = await this.#dependencies.startClient({
			binary: request.binary,
			cwd: request.cwd,
			signal: request.signal,
			request,
		});
		const onAbort = (): void => {
			void client.abort(cancelled(request.signal) ?? new Error("ProteinView panel launch cancelled"));
		};
		request.signal?.addEventListener("abort", onAbort, { once: true });
		try {
			const beforePresentation = cancelled(request.signal);
			if (beforePresentation !== undefined) throw beforePresentation;
			await this.#applyInitialPresentation(client, request);
			const afterPresentation = cancelled(request.signal);
			if (afterPresentation !== undefined) throw afterPresentation;
		} catch (error) {
			const cancellation = cancelled(request.signal);
			if (cancellation !== undefined) {
				await client.abort(cancellation);
			} else {
				await client.close();
			}
			throw error;
		} finally {
			request.signal?.removeEventListener("abort", onAbort);
		}

		this.#client = client;
		this.#context = context;
		this.#focused = false;
		const failure = client.failure;
		if (failure !== undefined) {
			await client.close();
			throw failure;
		}
		this.#failureUnsubscribe = client.subscribeFailure(error => {
			if (this.#client !== client) return;
			this.#focused = false;
			this.#autoRotateInFlight = false;
			this.#inputUnsubscribe?.();
			this.#inputUnsubscribe = undefined;
			this.#widget?.dispose();
			this.#widget = undefined;
			context.ui.setStatus(STATUS_KEY, undefined);
			context.ui.setWidget(WIDGET_KEY, undefined, { placement: "aboveEditor" });
			context.ui.notify(`ProteinView live panel stopped: ${visibleError(error)}`, "error");
		});
		context.ui.setWidget(
			WIDGET_KEY,
			(tui, theme) => {
				const widget = new ProteinViewPanelWidget({
					client,
					identifier: request.identifier,
					tui,
					theme,
				});
				widget.setFocused(this.#focused);
				this.#widget = widget;
				return widget;
			},
			{ placement: "aboveEditor" },
		);
		this.#inputUnsubscribe = context.ui.onTerminalInput(data => this.handleTerminalInput(data));
	}

	async #applyInitialPresentation(
		client: ProteinViewPanelClientLike,
		request: ProteinViewLiveOpenRequest,
	): Promise<void> {
		if (!request.showLigands) {
			await client.command({ command: "set_ligands", enabled: false });
		}
		if (request.interfaceChain === undefined) {
			if (request.showInteractions) {
				throw new Error("ProteinView interactions require an interface chain");
			}
			return;
		}
		await client.command({ command: "set_chain", chain: request.interfaceChain });
		await client.command({ command: "set_interface", enabled: true });
		if (request.showInteractions) {
			await client.command({ command: "set_interactions", enabled: true });
		}
	}

	async #close(): Promise<void> {
		const client = this.#client;
		const context = this.#context;
		this.#client = undefined;
		this.#context = undefined;
		this.#focused = false;
		this.#autoRotateInFlight = false;
		this.#failureUnsubscribe?.();
		this.#failureUnsubscribe = undefined;
		this.#inputUnsubscribe?.();
		this.#inputUnsubscribe = undefined;
		context?.ui.setStatus(STATUS_KEY, undefined);
		context?.ui.setWidget(WIDGET_KEY, undefined, { placement: "aboveEditor" });
		this.#widget?.dispose();
		this.#widget = undefined;
		await client?.close();
	}

	#setFocused(focused: boolean): void {
		this.#focused = focused;
		this.#widget?.setFocused(focused);
		this.#context?.ui.setStatus(STATUS_KEY, focused ? "ProteinView focused · Esc/Tab returns to chat" : undefined);
	}

	#dispatch(command: PanelCommand): void {
		void this.control(command).catch(error => {
			this.#widget?.setError(visibleError(error));
		});
	}
}

export function installProteinViewExtension(
	pi: ExtensionAPI,
	dependencies: ProteinViewExtensionDependencies = DEFAULT_DEPENDENCIES,
): ProteinViewPanelHost {
	const host = new ProteinViewPanelHost(dependencies);
	let unregisterOpenHandler: (() => void) | undefined;
	let autoRotateTimer: Timer | undefined;
	let timerContext: ExtensionContext | undefined;

	const deactivate = async (): Promise<void> => {
		unregisterOpenHandler?.();
		unregisterOpenHandler = undefined;
		if (autoRotateTimer !== undefined && timerContext !== undefined) {
			timerContext.clearTimer(autoRotateTimer);
		}
		autoRotateTimer = undefined;
		timerContext = undefined;
		await host.close();
	};

	const activate = async (context: ExtensionContext): Promise<void> => {
		await deactivate();
		if (!context.hasUI) return;
		unregisterOpenHandler = registerProteinViewLiveOpenHandler(request => host.open(request, context));
		timerContext = context;
		autoRotateTimer = context.setInterval(() => host.autoRotateTick(), AUTO_ROTATE_INTERVAL_MS);
	};

	pi.setLabel("ProteinView Live");
	pi.on("session_start", (_event, context) => activate(context));
	pi.on("session_switch", (_event, context) => activate(context));
	pi.on("session_branch", (_event, context) => activate(context));
	pi.on("session_tree", (_event, context) => activate(context));
	pi.on("session_shutdown", () => deactivate());

	pi.registerShortcut("f8", {
		description: "Focus or release the live ProteinView panel",
		handler: context => host.toggleFocus(context),
	});

	const z = pi.zod;
	const actions = [
		"rotate",
		"pan",
		"zoom",
		"reset",
		"fit",
		"set_color",
		"cycle_color",
		"set_viz",
		"cycle_viz",
		"select_chain",
		"set_chain",
		"set_interface",
		"set_interactions",
		"set_ligands",
		"set_auto_rotate",
		"set_residue_colors",
		"clear_residue_colors",
		"close",
		"get_state",
	] as const;
	const parameters = z.object({
		action: z.enum(actions),
		axis: z.enum(["x", "y", "z"] as const).optional(),
		delta: z.number().finite().optional(),
		dx: z.number().finite().optional(),
		dy: z.number().finite().optional(),
		factor: z.number().finite().positive().optional(),
		direction: z.enum(["in", "out", "prev", "next"] as const).optional(),
		color: z.enum(["structure", "element", "chain", "plddt", "bfactor", "rainbow"] as const).optional(),
		mode: z.enum(["cartoon", "backbone", "wireframe"] as const).optional(),
		chain: z
			.string()
			.regex(/^[!-~]{1,32}$/)
			.describe("Exact case-sensitive chain ID")
			.optional(),
		enabled: z.boolean().optional(),
		residues: z
			.array(
				z.object({
					chain: z
						.string()
						.regex(/^[!-~]{1,32}$/)
						.describe("Exact case-sensitive chain ID"),
					residue_number: z.number().int().min(-2_147_483_648).max(2_147_483_647),
					insertion_code: z
						.string()
						.regex(/^[!-~]{1,32}$/)
						.optional(),
					color: z
						.string()
						.regex(/^[0-9A-F]{6}$/)
						.describe("Uppercase RGB color in RRGGBB form"),
				}),
			)
			.max(MAX_RESIDUE_COLORS)
			.optional(),
	});

	pi.registerTool<typeof parameters, ProteinViewControlDetails>({
		name: "proteinview_control",
		label: "ProteinView control",
		description:
			"Control the currently open live ProteinView panel. Use exact chain/residue identifiers and uppercase RRGGBB colors for residue overrides.",
		parameters,
		strict: true,
		approval: "read",
		loadMode: "essential",
		async execute(_toolCallId, params, signal) {
			const cancellation = cancelled(signal);
			if (cancellation !== undefined) throw cancellation;
			if (params.action === "close") {
				await host.close();
				return {
					content: [{ type: "text", text: "Closed the live ProteinView panel." }],
					details: { closed: true },
				};
			}

			let command: PanelCommand;
			switch (params.action) {
				case "rotate":
					command = {
						command: "rotate",
						axis: requireField(params.axis, "axis"),
						delta: requireField(params.delta, "delta"),
					};
					break;
				case "pan":
					command = {
						command: "pan",
						dx: requireField(params.dx, "dx"),
						dy: requireField(params.dy, "dy"),
					};
					break;
				case "zoom": {
					const factor = params.factor;
					const direction = params.direction === "in" || params.direction === "out" ? params.direction : undefined;
					if ((factor === undefined) === (direction === undefined)) {
						throw new Error("zoom requires exactly one of factor or direction=in|out");
					}
					command =
						factor === undefined
							? { command: "zoom", direction: direction as "in" | "out" }
							: { command: "zoom", factor };
					break;
				}
				case "reset":
					command = { command: "reset" };
					break;
				case "fit":
					command = { command: "fit" };
					break;
				case "set_color":
					command = {
						command: "set_color",
						color: requireField(params.color, "color"),
					};
					break;
				case "cycle_color":
					command = { command: "cycle_color" };
					break;
				case "set_viz":
					command = {
						command: "set_viz",
						mode: requireField(params.mode, "mode"),
					};
					break;
				case "cycle_viz":
					command = { command: "cycle_viz" };
					break;
				case "select_chain": {
					const direction = requireField(params.direction, "direction");
					if (direction !== "prev" && direction !== "next") {
						throw new Error("select_chain direction must be prev or next");
					}
					command = { command: "select_chain", direction };
					break;
				}
				case "set_chain":
					command = {
						command: "set_chain",
						chain: requireField(params.chain, "chain"),
					};
					break;
				case "set_interface":
					command = {
						command: "set_interface",
						enabled: requireField(params.enabled, "enabled"),
					};
					break;
				case "set_interactions":
					command = {
						command: "set_interactions",
						enabled: requireField(params.enabled, "enabled"),
					};
					break;
				case "set_ligands":
					command = {
						command: "set_ligands",
						enabled: requireField(params.enabled, "enabled"),
					};
					break;
				case "set_auto_rotate":
					command = {
						command: "set_auto_rotate",
						enabled: requireField(params.enabled, "enabled"),
					};
					break;
				case "set_residue_colors": {
					const residues = requireField(params.residues, "residues");
					command = {
						command: "set_residue_colors",
						residues: residues.map(
							(residue): ProteinViewLiveResidueColor => ({
								chain: residue.chain,
								residue_number: residue.residue_number,
								...(residue.insertion_code === undefined ? {} : { insertion_code: residue.insertion_code }),
								color: residue.color,
							}),
						),
					};
					break;
				}
				case "clear_residue_colors":
					command = { command: "set_residue_colors", residues: [] };
					break;
				case "get_state":
					command = { command: "get_state" };
					break;
			}
			const result = await host.control(command, signal);
			return {
				content: [
					{
						type: "text",
						text:
							`ProteinView ${params.action} complete at revision ${result.revision}. ` +
							`Mode ${result.state.presentation.viz_mode}, color ${result.state.presentation.effective_color}.`,
					},
				],
				details: result,
			};
		},
	});

	return host;
}

export default function proteinViewExtension(pi: ExtensionAPI): void {
	installProteinViewExtension(pi);
}

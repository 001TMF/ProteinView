import * as fs from "node:fs/promises";
import * as os from "node:os";
import * as path from "node:path";
import { FULLHD_HEIGHT, FULLHD_WIDTH, MAX_PNG_BYTES, validateFullHdPng } from "../tools/index.ts";
import type { ProteinViewLiveOpenRequest, ProteinViewLiveResidueColor } from "./open-bridge.ts";

export const PANEL_PROTOCOL = "proteinview-panel";
export const PANEL_PROTOCOL_VERSION = 1;
export const PANEL_REQUEST_LIMIT_BYTES = 64 * 1024;
export const PANEL_RESPONSE_LIMIT_BYTES = 256 * 1024;
export const PANEL_STDERR_LIMIT_BYTES = 16 * 1024;
export const PANEL_READY_TIMEOUT_MS = 30_000;
export const PANEL_COMMAND_TIMEOUT_MS = 30_000;
export const PANEL_SHUTDOWN_TIMEOUT_MS = 2_000;

// biome-ignore lint/complexity/useRegexLiterals: Constructor strings avoid literal control escapes rejected by noControlCharactersInRegex.
const ANSI_ESCAPE = new RegExp("\\x1b(?:\\[[0-?]*[ -/]*[@-~]|\\][^\\x07]*(?:\\x07|\\x1b\\\\))", "g");
// biome-ignore lint/complexity/useRegexLiterals: Constructor strings avoid literal control escapes rejected by noControlCharactersInRegex.
const CONTROL_CHARACTERS = new RegExp("[\\u0000-\\u001f\\u007f-\\u009f]", "g");

export type PanelColor = "structure" | "element" | "chain" | "plddt" | "bfactor" | "rainbow";
export type PanelVizMode = "cartoon" | "backbone" | "wireframe";
export type PanelAxis = "x" | "y" | "z";
export type PanelDirection = "prev" | "next";

export type PanelCommand =
	| { command: "render" }
	| { command: "rotate"; axis: PanelAxis; delta: number }
	| { command: "pan"; dx: number; dy: number }
	| { command: "zoom"; factor: number }
	| { command: "zoom"; direction: "in" | "out" }
	| { command: "reset" }
	| { command: "fit" }
	| { command: "set_color"; color: PanelColor }
	| { command: "cycle_color" }
	| { command: "set_viz"; mode: PanelVizMode }
	| { command: "cycle_viz" }
	| { command: "select_chain"; direction: PanelDirection }
	| { command: "set_chain"; chain: string }
	| { command: "set_interface"; enabled: boolean }
	| { command: "set_interactions"; enabled: boolean }
	| { command: "set_ligands"; enabled: boolean }
	| { command: "set_auto_rotate"; enabled: boolean }
	| { command: "set_residue_colors"; residues: ProteinViewLiveResidueColor[] }
	| { command: "get_state" };

type InternalPanelCommand = PanelCommand | { command: "shutdown" };

export interface PanelStructureState {
	name: string;
	chain_count: number;
	residue_count: number;
	atom_count: number;
	ligand_count: number;
	chains: string[];
}

export interface PanelCameraState {
	rot_x: number;
	rot_y: number;
	rot_z: number;
	zoom: number;
	pan_x: number;
	pan_y: number;
	auto_rotate: boolean;
}

export interface PanelViewportState {
	width: number;
	height: number;
}

export interface PanelResidueColorState extends ProteinViewLiveResidueColor {
	insertion_code?: string;
	residue_name: string;
}

export interface PanelPresentationState {
	viz_mode: PanelVizMode;
	color: PanelColor;
	effective_color: PanelColor | "interface";
	current_chain_index: number;
	current_chain_id?: string;
	interface: boolean;
	interactions: boolean;
	ligands: boolean;
	residue_colors: PanelResidueColorState[];
}

export interface PanelState {
	structure: PanelStructureState;
	camera: PanelCameraState;
	viewport: PanelViewportState;
	presentation: PanelPresentationState;
}

export interface PanelFrame {
	path: string;
	mime_type: "image/png";
	width: number;
	height: number;
}

export interface PanelSnapshot {
	revision: number;
	state: PanelState;
	frame: PanelFrame;
	png: Buffer;
}

export interface PanelCommandResult {
	revision: number;
	state: PanelState;
	frame: PanelFrame;
}

export interface PanelSpawnOptions {
	cwd: string;
}

export interface PanelTransport {
	stdout: ReadableStream<Uint8Array>;
	stderr: ReadableStream<Uint8Array>;
	exited: Promise<number>;
	write(data: string): Promise<void>;
	closeInput(): void;
	kill(): void;
}

export interface PanelRuntime {
	createPrivateDir(): Promise<string>;
	removePrivateDir(directory: string): Promise<void>;
	readFile(filePath: string): Promise<Buffer>;
	spawn(argv: string[], options: PanelSpawnOptions): PanelTransport;
}

export interface PanelClientOptions {
	binary?: string;
	cwd: string;
	signal?: AbortSignal;
	request: ProteinViewLiveOpenRequest;
	runtime?: PanelRuntime;
	readyTimeoutMs?: number;
	commandTimeoutMs?: number;
	shutdownTimeoutMs?: number;
}

interface PendingResponse {
	id: string;
	command: InternalPanelCommand["command"];
	resolve: (result: PanelCommandResult) => void;
	reject: (error: Error) => void;
	timer: Timer;
}

interface ReadyRecord {
	revision: number;
	state: PanelState;
	frame: PanelFrame;
}

interface SuccessResponse {
	id: string;
	revision: number;
	state: PanelState;
	frame: PanelFrame;
}

interface ErrorResponse {
	id: string;
	revision: number;
	code: string;
	message: string;
}

type ParsedResponse = SuccessResponse | ErrorResponse;
type SnapshotListener = (snapshot: PanelSnapshot) => void;
type FailureListener = (error: Error) => void;

function record(value: unknown, label: string): Record<string, unknown> {
	if (typeof value !== "object" || value === null || Array.isArray(value)) {
		throw new Error(`${label} must be an object`);
	}
	return value as Record<string, unknown>;
}

function exactKeys(
	value: Record<string, unknown>,
	required: readonly string[],
	optional: readonly string[],
	label: string,
): void {
	const allowed = new Set([...required, ...optional]);
	for (const key of required) {
		if (!(key in value)) throw new Error(`${label} is missing ${key}`);
	}
	for (const key of Object.keys(value)) {
		if (!allowed.has(key)) throw new Error(`${label} contains unsupported field ${key}`);
	}
}

function stringValue(value: unknown, label: string, allowEmpty = false): string {
	if (typeof value !== "string" || (!allowEmpty && value.length === 0)) {
		throw new Error(`${label} must be ${allowEmpty ? "a string" : "a non-empty string"}`);
	}
	return value;
}

function finiteNumber(value: unknown, label: string): number {
	if (typeof value !== "number" || !Number.isFinite(value)) {
		throw new Error(`${label} must be a finite number`);
	}
	return value;
}

function nonNegativeInteger(value: unknown, label: string): number {
	if (typeof value !== "number" || !Number.isSafeInteger(value) || value < 0) {
		throw new Error(`${label} must be a non-negative safe integer`);
	}
	return value;
}

function booleanValue(value: unknown, label: string): boolean {
	if (typeof value !== "boolean") throw new Error(`${label} must be a boolean`);
	return value;
}

function oneOf<T extends string>(value: unknown, allowed: readonly T[], label: string): T {
	if (typeof value !== "string" || !allowed.includes(value as T)) {
		throw new Error(`${label} has unsupported value ${String(value)}`);
	}
	return value as T;
}

function nullableString(value: unknown, label: string): string | undefined {
	if (value === null || value === undefined) return undefined;
	return stringValue(value, label);
}

function parseState(value: unknown): PanelState {
	const state = record(value, "panel state");
	exactKeys(state, ["structure", "camera", "viewport", "presentation"], [], "panel state");

	const structure = record(state.structure, "panel structure");
	exactKeys(
		structure,
		["name", "chain_count", "residue_count", "atom_count", "ligand_count", "chains"],
		[],
		"panel structure",
	);
	if (!Array.isArray(structure.chains)) {
		throw new Error("panel structure chains must be an array");
	}
	const chains = structure.chains.map((chain, index) => stringValue(chain, `panel structure chains[${index}]`));

	const camera = record(state.camera, "panel camera");
	exactKeys(camera, ["rot_x", "rot_y", "rot_z", "zoom", "pan_x", "pan_y", "auto_rotate"], [], "panel camera");

	const viewport = record(state.viewport, "panel viewport");
	exactKeys(viewport, ["width", "height"], [], "panel viewport");

	const presentation = record(state.presentation, "panel presentation");
	exactKeys(
		presentation,
		[
			"viz_mode",
			"color",
			"effective_color",
			"current_chain_index",
			"current_chain_id",
			"interface",
			"interactions",
			"ligands",
			"residue_colors",
		],
		[],
		"panel presentation",
	);
	if (!Array.isArray(presentation.residue_colors)) {
		throw new Error("panel presentation residue_colors must be an array");
	}
	const residueColors = presentation.residue_colors.map((entry, index) => {
		const residue = record(entry, `panel residue_colors[${index}]`);
		exactKeys(
			residue,
			["chain", "residue_number", "insertion_code", "residue_name", "color"],
			[],
			`panel residue_colors[${index}]`,
		);
		const residueNumber = finiteNumber(residue.residue_number, `panel residue_colors[${index}].residue_number`);
		if (!Number.isInteger(residueNumber)) {
			throw new Error(`panel residue_colors[${index}].residue_number must be an integer`);
		}
		const color = stringValue(residue.color, `panel residue_colors[${index}].color`);
		if (!/^[0-9A-F]{6}$/.test(color)) {
			throw new Error(`panel residue_colors[${index}].color must be uppercase RRGGBB`);
		}
		const insertionCode = nullableString(residue.insertion_code, `panel residue_colors[${index}].insertion_code`);
		return {
			chain: stringValue(residue.chain, `panel residue_colors[${index}].chain`),
			residue_number: residueNumber,
			...(insertionCode === undefined ? {} : { insertion_code: insertionCode }),
			residue_name: stringValue(residue.residue_name, `panel residue_colors[${index}].residue_name`),
			color,
		};
	});

	const currentChainId = nullableString(presentation.current_chain_id, "panel presentation current_chain_id");
	return {
		structure: {
			name: stringValue(structure.name, "panel structure name", true),
			chain_count: nonNegativeInteger(structure.chain_count, "panel structure chain_count"),
			residue_count: nonNegativeInteger(structure.residue_count, "panel structure residue_count"),
			atom_count: nonNegativeInteger(structure.atom_count, "panel structure atom_count"),
			ligand_count: nonNegativeInteger(structure.ligand_count, "panel structure ligand_count"),
			chains,
		},
		camera: {
			rot_x: finiteNumber(camera.rot_x, "panel camera rot_x"),
			rot_y: finiteNumber(camera.rot_y, "panel camera rot_y"),
			rot_z: finiteNumber(camera.rot_z, "panel camera rot_z"),
			zoom: finiteNumber(camera.zoom, "panel camera zoom"),
			pan_x: finiteNumber(camera.pan_x, "panel camera pan_x"),
			pan_y: finiteNumber(camera.pan_y, "panel camera pan_y"),
			auto_rotate: booleanValue(camera.auto_rotate, "panel camera auto_rotate"),
		},
		viewport: {
			width: nonNegativeInteger(viewport.width, "panel viewport width"),
			height: nonNegativeInteger(viewport.height, "panel viewport height"),
		},
		presentation: {
			viz_mode: oneOf(
				presentation.viz_mode,
				["cartoon", "backbone", "wireframe"] as const,
				"panel presentation viz_mode",
			),
			color: oneOf(
				presentation.color,
				["structure", "element", "chain", "plddt", "bfactor", "rainbow"] as const,
				"panel presentation color",
			),
			effective_color: oneOf(
				presentation.effective_color,
				["structure", "element", "chain", "plddt", "bfactor", "rainbow", "interface"] as const,
				"panel presentation effective_color",
			),
			current_chain_index: nonNegativeInteger(
				presentation.current_chain_index,
				"panel presentation current_chain_index",
			),
			...(currentChainId === undefined ? {} : { current_chain_id: currentChainId }),
			interface: booleanValue(presentation.interface, "panel presentation interface"),
			interactions: booleanValue(presentation.interactions, "panel presentation interactions"),
			ligands: booleanValue(presentation.ligands, "panel presentation ligands"),
			residue_colors: residueColors,
		},
	};
}

function parseFrame(value: unknown, expectedPath: string): PanelFrame {
	const frame = record(value, "panel frame");
	exactKeys(frame, ["path", "mime_type", "width", "height"], [], "panel frame");
	const framePath = stringValue(frame.path, "panel frame path");
	if (path.resolve(framePath) !== path.resolve(expectedPath)) {
		throw new Error("panel frame path does not match the private output path");
	}
	const mimeType = stringValue(frame.mime_type, "panel frame mime_type");
	if (mimeType !== "image/png") throw new Error("panel frame mime_type must be image/png");
	const width = nonNegativeInteger(frame.width, "panel frame width");
	const height = nonNegativeInteger(frame.height, "panel frame height");
	if (width !== FULLHD_WIDTH || height !== FULLHD_HEIGHT) {
		throw new Error(`panel frame dimensions were ${width}x${height}; expected ${FULLHD_WIDTH}x${FULLHD_HEIGHT}`);
	}
	return { path: framePath, mime_type: "image/png", width, height };
}

function parseReady(value: unknown, expectedPath: string): ReadyRecord {
	const ready = record(value, "panel ready record");
	exactKeys(ready, ["type", "protocol", "version", "revision", "state", "frame", "limits"], [], "panel ready record");
	if (ready.type !== "ready") throw new Error("first panel record must be ready");
	if (ready.protocol !== PANEL_PROTOCOL) throw new Error("unexpected panel protocol name");
	if (ready.version !== PANEL_PROTOCOL_VERSION) {
		throw new Error(`unsupported panel protocol version ${String(ready.version)}`);
	}
	const revision = nonNegativeInteger(ready.revision, "panel ready revision");
	if (revision !== 1) throw new Error(`initial panel revision must be 1, received ${revision}`);
	const limits = record(ready.limits, "panel limits");
	exactKeys(limits, ["min_dimension", "max_dimension", "max_pixels", "max_request_bytes"], [], "panel limits");
	if (limits.max_request_bytes !== PANEL_REQUEST_LIMIT_BYTES) {
		throw new Error("panel request byte limit does not match protocol v1");
	}
	return {
		revision,
		state: parseState(ready.state),
		frame: parseFrame(ready.frame, expectedPath),
	};
}

function parseResponse(value: unknown, expectedPath: string): ParsedResponse {
	const response = record(value, "panel response");
	const ok = booleanValue(response.ok, "panel response ok");
	exactKeys(
		response,
		ok
			? ["type", "protocol", "version", "id", "ok", "revision", "state", "frame"]
			: ["type", "protocol", "version", "id", "ok", "revision", "error"],
		[],
		"panel response",
	);
	if (response.type !== "response") throw new Error("panel emitted a non-response record");
	if (response.protocol !== PANEL_PROTOCOL) throw new Error("unexpected panel protocol name");
	if (response.version !== PANEL_PROTOCOL_VERSION) {
		throw new Error(`unsupported panel protocol version ${String(response.version)}`);
	}
	const id = stringValue(response.id, "panel response id");
	const revision = nonNegativeInteger(response.revision, "panel response revision");
	if (ok) {
		return {
			id,
			revision,
			state: parseState(response.state),
			frame: parseFrame(response.frame, expectedPath),
		};
	}
	const error = record(response.error, "panel response error");
	exactKeys(error, ["code", "message"], [], "panel response error");
	return {
		id,
		revision,
		code: stringValue(error.code, "panel response error code"),
		message: stringValue(error.message, "panel response error message"),
	};
}

function isErrorResponse(response: ParsedResponse): response is ErrorResponse {
	return "code" in response;
}

function sanitizeStderr(value: string): string {
	return value.replace(ANSI_ESCAPE, "").replace(CONTROL_CHARACTERS, " ").trim();
}

function processError(message: string, stderr: string): Error {
	const suffix = sanitizeStderr(stderr);
	return new Error(suffix ? `${message}: ${suffix}` : message);
}

class BoundedByteTail {
	#bytes = Buffer.alloc(0);

	constructor(readonly limit: number) {}

	append(value: Uint8Array): void {
		if (value.byteLength >= this.limit) {
			this.#bytes = Buffer.from(value.subarray(value.byteLength - this.limit));
			return;
		}
		const combined = Buffer.concat([this.#bytes, Buffer.from(value)]);
		this.#bytes = combined.byteLength <= this.limit ? combined : combined.subarray(combined.byteLength - this.limit);
	}

	text(): string {
		return this.#bytes.toString("utf8");
	}
}

async function consumeCappedLines(
	stream: ReadableStream<Uint8Array>,
	limit: number,
	onLine: (line: string) => Promise<void>,
): Promise<void> {
	const reader = stream.getReader();
	const decoder = new TextDecoder("utf-8", { fatal: true });
	let pending = Buffer.alloc(0);
	try {
		while (true) {
			const { done, value } = await reader.read();
			if (done) break;
			if (value === undefined) continue;
			const chunk = Buffer.from(value);
			let offset = 0;
			while (offset < chunk.byteLength) {
				const newline = chunk.indexOf(0x0a, offset);
				const end = newline === -1 ? chunk.byteLength : newline;
				const segment = chunk.subarray(offset, end);
				if (pending.byteLength + segment.byteLength > limit) {
					throw new Error(`panel NDJSON record exceeds the ${limit}-byte response limit`);
				}
				pending = pending.byteLength === 0 ? Buffer.from(segment) : Buffer.concat([pending, segment]);
				if (newline === -1) break;
				if (pending.at(-1) === 0x0d) pending = pending.subarray(0, -1);
				if (pending.byteLength === 0) throw new Error("panel emitted a blank NDJSON record");
				const line = decoder.decode(pending);
				pending = Buffer.alloc(0);
				await onLine(line);
				offset = newline + 1;
			}
		}
		if (pending.byteLength !== 0) {
			throw new Error("panel stdout ended with an unterminated NDJSON record");
		}
	} finally {
		reader.releaseLock();
	}
}

async function drainStderr(stream: ReadableStream<Uint8Array>, tail: BoundedByteTail): Promise<void> {
	const reader = stream.getReader();
	try {
		while (true) {
			const { done, value } = await reader.read();
			if (done) break;
			if (value !== undefined) tail.append(value);
		}
	} finally {
		reader.releaseLock();
	}
}

function requestLine(id: string, command: InternalPanelCommand): string {
	const line = `${JSON.stringify({
		version: PANEL_PROTOCOL_VERSION,
		id,
		...command,
	})}\n`;
	if (Buffer.byteLength(line, "utf8") > PANEL_REQUEST_LIMIT_BYTES) {
		throw new Error(`panel request exceeds the ${PANEL_REQUEST_LIMIT_BYTES}-byte protocol limit`);
	}
	return line;
}

function commandRenders(command: InternalPanelCommand["command"]): boolean {
	return command !== "get_state" && command !== "shutdown";
}

function deadline<T>(promise: Promise<T>, ms: number, onTimeout: () => void, message: string): Promise<T> {
	const timeout = Promise.withResolvers<T>();
	const timer = setTimeout(() => {
		onTimeout();
		timeout.reject(new Error(message));
	}, ms);
	return Promise.race([promise, timeout.promise]).finally(() => clearTimeout(timer));
}

export class PanelProtocolError extends Error {
	readonly code: string;

	constructor(code: string, message: string) {
		super(message);
		this.name = "PanelProtocolError";
		this.code = code;
	}
}

export class ProteinViewPanelClient {
	#transport: PanelTransport;
	#runtime: PanelRuntime;
	#privateDir: string;
	#outputPath: string;
	#commandTimeoutMs: number;
	#shutdownTimeoutMs: number;
	#ready = Promise.withResolvers<ReadyRecord>();
	#pending: PendingResponse | undefined;
	#queue: Promise<void> = Promise.resolve();
	#listeners = new Set<SnapshotListener>();
	#failureListeners = new Set<FailureListener>();
	#stderr = new BoundedByteTail(PANEL_STDERR_LIMIT_BYTES);
	#stdoutLoop: Promise<void>;
	#stderrLoop: Promise<void>;
	#revision = 0;
	#nextId = 1;
	#snapshot: PanelSnapshot | undefined;
	#fault: Error | undefined;
	#closing = false;
	#closed = false;
	#expectedExit = false;
	#cleanupPromise: Promise<void> | undefined;

	static async start(options: PanelClientOptions): Promise<ProteinViewPanelClient> {
		if (options.signal?.aborted) {
			throw options.signal.reason instanceof Error
				? options.signal.reason
				: new Error("ProteinView panel launch cancelled");
		}
		const runtime = options.runtime ?? createDefaultPanelRuntime();
		const privateDir = await runtime.createPrivateDir();
		const outputPath = path.join(privateDir, "frame.png");
		const binary = options.binary ?? options.request.binary ?? (Bun.env.PROTEINVIEW_BIN?.trim() || "proteinview");
		const argv = buildPanelArgv(binary, outputPath, options.request);
		let transport: PanelTransport;
		try {
			transport = runtime.spawn(argv, { cwd: options.cwd });
		} catch (error) {
			await runtime.removePrivateDir(privateDir);
			throw processError(
				`unable to launch ProteinView panel server: ${error instanceof Error ? error.message : String(error)}`,
				"",
			);
		}
		const client = new ProteinViewPanelClient(
			transport,
			runtime,
			privateDir,
			outputPath,
			options.commandTimeoutMs ?? PANEL_COMMAND_TIMEOUT_MS,
			options.shutdownTimeoutMs ?? PANEL_SHUTDOWN_TIMEOUT_MS,
		);
		try {
			const cancelled = Promise.withResolvers<ReadyRecord>();
			const onAbort = (): void => {
				const error =
					options.signal?.reason instanceof Error
						? options.signal.reason
						: new Error("ProteinView panel launch cancelled");
				client.#fail(error, true);
				cancelled.reject(error);
			};
			options.signal?.addEventListener("abort", onAbort, { once: true });
			let ready: ReadyRecord;
			try {
				ready = await deadline(
					Promise.race([client.#ready.promise, cancelled.promise]),
					options.readyTimeoutMs ?? PANEL_READY_TIMEOUT_MS,
					() => client.#fail(new Error("ProteinView panel ready handshake timed out"), true),
					"ProteinView panel ready handshake timed out",
				);
				await client.#acceptFrame(ready.revision, ready.state, ready.frame);
			} finally {
				options.signal?.removeEventListener("abort", onAbort);
			}
			return client;
		} catch (error) {
			await client.close();
			throw error;
		}
	}

	constructor(
		transport: PanelTransport,
		runtime: PanelRuntime,
		privateDir: string,
		outputPath: string,
		commandTimeoutMs: number,
		shutdownTimeoutMs: number,
	) {
		this.#transport = transport;
		this.#runtime = runtime;
		this.#privateDir = privateDir;
		this.#outputPath = outputPath;
		this.#commandTimeoutMs = commandTimeoutMs;
		this.#shutdownTimeoutMs = shutdownTimeoutMs;
		this.#stdoutLoop = this.#runStdout();
		this.#stderrLoop = drainStderr(transport.stderr, this.#stderr).catch(() => {});
		void transport.exited.then(exitCode => {
			if (this.#expectedExit || this.#closing || this.#closed) return;
			this.#fail(
				processError(`ProteinView panel server exited unexpectedly with code ${exitCode}`, this.#stderr.text()),
				false,
			);
		});
	}

	get snapshot(): PanelSnapshot | undefined {
		return this.#snapshot;
	}

	get state(): PanelState | undefined {
		return this.#snapshot?.state;
	}

	get closed(): boolean {
		return this.#closed || this.#fault !== undefined;
	}

	get failure(): Error | undefined {
		return this.#fault;
	}

	subscribe(listener: SnapshotListener): () => void {
		this.#listeners.add(listener);
		const current = this.#snapshot;
		if (current !== undefined) listener(current);
		return () => {
			this.#listeners.delete(listener);
		};
	}

	subscribeFailure(listener: FailureListener): () => void {
		this.#failureListeners.add(listener);
		const failure = this.#fault;
		if (failure !== undefined) listener(failure);
		return () => {
			this.#failureListeners.delete(listener);
		};
	}

	command(command: PanelCommand): Promise<PanelCommandResult> {
		if (this.#fault !== undefined) {
			return Promise.reject(this.#fault);
		}
		if (this.#closing || this.#closed) {
			return Promise.reject(new Error("ProteinView panel is closed"));
		}
		return this.#enqueue(command, false);
	}

	async close(): Promise<void> {
		if (this.#cleanupPromise !== undefined) return this.#cleanupPromise;
		this.#cleanupPromise = this.#close();
		return this.#cleanupPromise;
	}

	abort(error: Error = new Error("ProteinView panel cancelled")): Promise<void> {
		this.#fail(error, true);
		return this.close();
	}

	async #close(): Promise<void> {
		if (this.#closed) return;
		this.#closing = true;
		try {
			if (this.#fault === undefined && this.#revision > 0) {
				const shutdownTimeoutMessage = "ProteinView panel shutdown timed out";
				try {
					await deadline(
						(async () => {
							await this.#enqueue({ command: "shutdown" }, true);
							this.#expectedExit = true;
							await this.#transport.exited;
						})(),
						this.#shutdownTimeoutMs,
						() => this.#fail(new Error(shutdownTimeoutMessage), true),
						shutdownTimeoutMessage,
					);
				} catch {
					this.#expectedExit = true;
					this.#transport.kill();
				}
			} else {
				this.#expectedExit = true;
				this.#transport.kill();
			}
			this.#transport.closeInput();
			await deadline(
				this.#transport.exited.then(() => undefined),
				this.#shutdownTimeoutMs,
				() => this.#transport.kill(),
				"ProteinView panel process did not exit",
			).catch(() => {});
			await Promise.allSettled([this.#stdoutLoop, this.#stderrLoop]);
		} finally {
			this.#closed = true;
			this.#listeners.clear();
			this.#failureListeners.clear();
			await this.#runtime.removePrivateDir(this.#privateDir);
		}
	}

	#enqueue(command: InternalPanelCommand, allowClosing: boolean): Promise<PanelCommandResult> {
		const result = this.#queue.then(() => this.#send(command, allowClosing));
		this.#queue = result.then(
			() => undefined,
			() => undefined,
		);
		return result;
	}

	async #send(command: InternalPanelCommand, allowClosing: boolean): Promise<PanelCommandResult> {
		if (this.#fault !== undefined) throw this.#fault;
		if (this.#closed || (this.#closing && !allowClosing)) {
			throw new Error("ProteinView panel is closed");
		}
		if (this.#pending !== undefined) {
			throw new Error("ProteinView panel attempted concurrent protocol requests");
		}
		const id = `omp-${this.#nextId}`;
		this.#nextId += 1;
		const response = Promise.withResolvers<PanelCommandResult>();
		const timer = setTimeout(() => {
			const error = new Error(
				`ProteinView panel command ${command.command} timed out after ${this.#commandTimeoutMs}ms`,
			);
			this.#fail(error, true);
		}, this.#commandTimeoutMs);
		this.#pending = {
			id,
			command: command.command,
			resolve: response.resolve,
			reject: response.reject,
			timer,
		};
		try {
			await this.#transport.write(requestLine(id, command));
		} catch (error) {
			const writeError = processError(
				`failed to write ProteinView panel command: ${error instanceof Error ? error.message : String(error)}`,
				this.#stderr.text(),
			);
			this.#fail(writeError, true);
		}
		return response.promise;
	}

	async #runStdout(): Promise<void> {
		try {
			await consumeCappedLines(this.#transport.stdout, PANEL_RESPONSE_LIMIT_BYTES, async line => {
				let value: unknown;
				try {
					value = JSON.parse(line);
				} catch (error) {
					throw new Error(
						`ProteinView panel emitted invalid JSON: ${error instanceof Error ? error.message : String(error)}`,
					);
				}
				if (this.#revision === 0) {
					const ready = parseReady(value, this.#outputPath);
					this.#revision = ready.revision;
					this.#ready.resolve(ready);
					return;
				}
				const pending = this.#pending;
				if (pending === undefined) {
					throw new Error("ProteinView panel emitted an unsolicited response");
				}
				const response = parseResponse(value, this.#outputPath);
				if (response.id !== pending.id) {
					throw new Error(`ProteinView panel response id ${response.id} did not match ${pending.id}`);
				}
				if (isErrorResponse(response)) {
					if (response.revision !== this.#revision) {
						throw new Error("failed panel command changed the frame revision");
					}
					clearTimeout(pending.timer);
					this.#pending = undefined;
					pending.reject(new PanelProtocolError(response.code, response.message));
					return;
				}
				const expectedRevision = commandRenders(pending.command) ? this.#revision + 1 : this.#revision;
				if (response.revision !== expectedRevision) {
					throw new Error(
						`ProteinView panel revision ${response.revision} did not match expected ${expectedRevision}`,
					);
				}
				if (response.revision > this.#revision) {
					await this.#acceptFrame(response.revision, response.state, response.frame, pending);
				} else if (this.#snapshot !== undefined) {
					this.#snapshot = {
						...this.#snapshot,
						state: response.state,
						frame: response.frame,
					};
				}
				if (this.#pending !== pending || this.#fault !== undefined) {
					throw this.#fault ?? new Error("ProteinView panel command was cancelled");
				}
				clearTimeout(pending.timer);
				this.#pending = undefined;
				pending.resolve({
					revision: response.revision,
					state: response.state,
					frame: response.frame,
				});
			});
			if (!this.#closing && !this.#closed) {
				throw processError("ProteinView panel stdout closed unexpectedly", this.#stderr.text());
			}
		} catch (error) {
			if (!this.#closing && !this.#closed) {
				this.#fail(error instanceof Error ? error : new Error(String(error)), true);
			}
		}
	}

	async #acceptFrame(
		revision: number,
		state: PanelState,
		frame: PanelFrame,
		pending?: PendingResponse,
	): Promise<void> {
		const bytes = await this.#runtime.readFile(this.#outputPath);
		if (bytes.byteLength > MAX_PNG_BYTES) {
			throw new Error(`ProteinView panel PNG exceeds the ${MAX_PNG_BYTES}-byte safety limit`);
		}
		validateFullHdPng(bytes);
		if (this.#fault !== undefined || (pending !== undefined && this.#pending !== pending)) {
			throw this.#fault ?? new Error("ProteinView panel command was cancelled");
		}
		this.#revision = revision;
		const snapshot = { revision, state, frame, png: bytes };
		this.#snapshot = snapshot;
		for (const listener of this.#listeners) listener(snapshot);
	}

	#fail(error: Error, kill: boolean): void {
		if (this.#fault !== undefined) return;
		this.#fault = error;
		this.#ready.reject(error);
		const pending = this.#pending;
		if (pending !== undefined) {
			clearTimeout(pending.timer);
			this.#pending = undefined;
			pending.reject(error);
		}
		for (const listener of this.#failureListeners) {
			try {
				listener(error);
			} catch {
				// A UI failure listener must not prevent process teardown.
			}
		}
		this.#failureListeners.clear();
		if (kill) {
			this.#expectedExit = true;
			this.#transport.kill();
		}
		void this.close();
	}
}

export function buildPanelArgv(binary: string, outputPath: string, request: ProteinViewLiveOpenRequest): string[] {
	const xyzDefault = path.extname(request.inputPath).toLowerCase() === ".xyz" && request.explicitColor !== true;
	const color = xyzDefault ? "element" : request.color;
	const args = [
		binary,
		request.inputPath,
		"--panel-server",
		"--output",
		outputPath,
		"--panel-width",
		String(FULLHD_WIDTH),
		"--panel-height",
		String(FULLHD_HEIGHT),
		"--color",
		color,
	];
	if (request.mode !== undefined) args.push("--mode", request.mode);
	for (const residue of request.residueColors) {
		const insertion = residue.insertion_code ? `[${residue.insertion_code}]` : "";
		args.push("--residue-color", `${residue.chain}:${residue.residue_number}${insertion}=${residue.color}`);
	}
	return args;
}

export function createDefaultPanelRuntime(): PanelRuntime {
	return {
		createPrivateDir: () => fs.mkdtemp(path.join(os.tmpdir(), "omp-proteinview-live-")),
		removePrivateDir: directory => fs.rm(directory, { recursive: true, force: true }),
		readFile: async filePath => {
			const file = Bun.file(filePath);
			if (file.size > MAX_PNG_BYTES) {
				throw new Error(`ProteinView panel PNG exceeds the ${MAX_PNG_BYTES}-byte safety limit`);
			}
			return Buffer.from(await file.arrayBuffer());
		},
		spawn: (argv, options) => {
			const process = Bun.spawn({
				cmd: argv,
				cwd: options.cwd,
				env: Bun.env,
				stdin: "pipe",
				stdout: "pipe",
				stderr: "pipe",
				windowsHide: true,
			});
			const input = process.stdin;
			return {
				stdout: process.stdout as ReadableStream<Uint8Array>,
				stderr: process.stderr as ReadableStream<Uint8Array>,
				exited: process.exited,
				write: async data => {
					input.write(data);
					await input.flush();
				},
				closeInput: () => {
					input.end();
				},
				kill: () => {
					try {
						process.kill();
					} catch {
						// The process already exited.
					}
				},
			};
		},
	};
}

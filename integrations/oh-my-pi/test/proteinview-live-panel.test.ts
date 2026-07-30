import { afterEach, describe, expect, it } from "bun:test";
import { deflateSync } from "node:zlib";
import type { ExtensionContext } from "@oh-my-pi/pi-coding-agent";
import type { Component, ImageDimensions, ImageOptions, ImageTheme, TUI } from "@oh-my-pi/pi-tui";
import { ProteinViewPanelHost } from "../extensions/index.ts";
import type { ProteinViewLiveOpenRequest } from "../live/open-bridge.ts";
import {
	buildPanelArgv,
	PANEL_PROTOCOL,
	PANEL_PROTOCOL_VERSION,
	type PanelCommand,
	type PanelRuntime,
	type PanelSnapshot,
	type PanelState,
	type PanelTransport,
	ProteinViewPanelClient,
} from "../live/panel-client.ts";
import {
	type ProteinViewPanelImage,
	type ProteinViewPanelImageFactory,
	ProteinViewPanelWidget,
} from "../live/panel-widget.ts";
import { FULLHD_HEIGHT, FULLHD_WIDTH } from "../tools/index.ts";

const encoder = new TextEncoder();
const cleanups: Array<() => Promise<void>> = [];

afterEach(async () => {
	await Promise.all(cleanups.splice(0).map(cleanup => cleanup()));
});

function crc32(bytes: Buffer): number {
	let crc = 0xffffffff;
	for (const byte of bytes) {
		crc ^= byte;
		for (let bit = 0; bit < 8; bit += 1) {
			crc = (crc >>> 1) ^ (crc & 1 ? 0xedb88320 : 0);
		}
	}
	return (crc ^ 0xffffffff) >>> 0;
}

function pngChunk(type: string, data: Buffer): Buffer {
	const typeBytes = Buffer.from(type, "ascii");
	const chunk = Buffer.alloc(12 + data.length);
	chunk.writeUInt32BE(data.length, 0);
	typeBytes.copy(chunk, 4);
	data.copy(chunk, 8);
	chunk.writeUInt32BE(crc32(Buffer.concat([typeBytes, data])), 8 + data.length);
	return chunk;
}

function fullHdPng(fill = 0): Buffer {
	const signature = Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]);
	const ihdr = Buffer.alloc(13);
	ihdr.writeUInt32BE(FULLHD_WIDTH, 0);
	ihdr.writeUInt32BE(FULLHD_HEIGHT, 4);
	ihdr[8] = 1;
	const rowBytes = Math.ceil(FULLHD_WIDTH / 8);
	const pixels = Buffer.alloc((rowBytes + 1) * FULLHD_HEIGHT, fill);
	return Buffer.concat([
		signature,
		pngChunk("IHDR", ihdr),
		pngChunk("IDAT", deflateSync(pixels)),
		pngChunk("IEND", Buffer.alloc(0)),
	]);
}

function initialState(): PanelState {
	return {
		structure: {
			name: "Test protein",
			chain_count: 2,
			residue_count: 4,
			atom_count: 20,
			ligand_count: 1,
			chains: ["A", "B"],
		},
		camera: {
			rot_x: 0,
			rot_y: 0,
			rot_z: 0,
			zoom: 1,
			pan_x: 0,
			pan_y: 0,
			auto_rotate: false,
		},
		viewport: { width: FULLHD_WIDTH, height: FULLHD_HEIGHT },
		presentation: {
			viz_mode: "cartoon",
			color: "chain",
			effective_color: "chain",
			current_chain_index: 0,
			current_chain_id: "A",
			interface: false,
			interactions: false,
			ligands: true,
			residue_colors: [],
		},
	};
}

function request(): ProteinViewLiveOpenRequest {
	return {
		inputPath: "/workspace/input.pdb",
		identifier: "input.pdb",
		cwd: "/workspace",
		mode: "cartoon",
		color: "chain",
		showInteractions: false,
		showLigands: true,
		residueColors: [],
	};
}

class FakePanelProcess implements PanelTransport {
	#stdout = new TransformStream<Uint8Array>();
	#stderr = new TransformStream<Uint8Array>();
	#stdoutWriter = this.#stdout.writable.getWriter();
	#stderrWriter = this.#stderr.writable.getWriter();
	#exit = Promise.withResolvers<number>();
	#outputPath: string;
	#revision = 1;
	#state = initialState();
	#closedStreams = false;
	requests: Array<Record<string, unknown>> = [];
	killed = false;
	closeInputCalled = false;
	autoRespond = true;
	heldRequests: Array<Record<string, unknown>> = [];

	readonly stdout = this.#stdout.readable;
	readonly stderr = this.#stderr.readable;
	readonly exited = this.#exit.promise;

	constructor(outputPath: string) {
		this.#outputPath = outputPath;
	}

	async ready(): Promise<void> {
		await this.#emit({
			type: "ready",
			protocol: PANEL_PROTOCOL,
			version: PANEL_PROTOCOL_VERSION,
			revision: this.#revision,
			state: this.#state,
			frame: this.#frame(),
			limits: {
				min_dimension: 64,
				max_dimension: 8192,
				max_pixels: 67_108_864,
				max_request_bytes: 65_536,
			},
		});
	}

	async write(data: string): Promise<void> {
		const parsed = JSON.parse(data) as Record<string, unknown>;
		this.requests.push(parsed);
		if (!this.autoRespond) {
			this.heldRequests.push(parsed);
			return;
		}
		await this.respond(parsed);
	}

	closeInput(): void {
		this.closeInputCalled = true;
	}

	kill(): void {
		this.killed = true;
		this.#exit.resolve(137);
		void this.#closeStreams();
	}

	async exitUnexpectedly(code = 1): Promise<void> {
		this.#exit.resolve(code);
		await this.#closeStreams();
	}

	async respond(requestRecord: Record<string, unknown>): Promise<void> {
		const command = requestRecord.command as PanelCommand["command"] | "shutdown";
		if (command !== "get_state" && command !== "shutdown") {
			this.#revision += 1;
			this.#apply(command, requestRecord);
		}
		if (command === "shutdown") {
			// Deliberately publish exit before the final response is consumed. A
			// graceful close must not misclassify this as an unexpected death.
			this.#exit.resolve(0);
		}
		await this.#emit({
			type: "response",
			protocol: PANEL_PROTOCOL,
			version: PANEL_PROTOCOL_VERSION,
			id: requestRecord.id,
			ok: true,
			revision: this.#revision,
			state: this.#state,
			frame: this.#frame(),
		});
		if (command === "shutdown") await this.#closeStreams();
	}

	async emitRaw(value: string): Promise<void> {
		await this.#stdoutWriter.write(encoder.encode(value));
	}

	async emitStderr(value: string): Promise<void> {
		await this.#stderrWriter.write(encoder.encode(value));
	}

	async #emit(value: unknown): Promise<void> {
		await this.#stdoutWriter.write(encoder.encode(`${JSON.stringify(value)}\n`));
	}

	#frame() {
		return {
			path: this.#outputPath,
			mime_type: "image/png",
			width: FULLHD_WIDTH,
			height: FULLHD_HEIGHT,
		};
	}

	#apply(command: PanelCommand["command"], values: Record<string, unknown>): void {
		if (command === "set_auto_rotate") {
			this.#state = {
				...this.#state,
				camera: {
					...this.#state.camera,
					auto_rotate: values.enabled === true,
				},
			};
		}
		if (command === "cycle_viz") {
			this.#state = {
				...this.#state,
				presentation: {
					...this.#state.presentation,
					viz_mode: "backbone",
				},
			};
		}
	}

	async #closeStreams(): Promise<void> {
		if (this.#closedStreams) return;
		this.#closedStreams = true;
		await Promise.allSettled([this.#stdoutWriter.close(), this.#stderrWriter.close()]);
	}
}

class FakeRuntime implements PanelRuntime {
	transport: FakePanelProcess | undefined;
	argv: string[] = [];
	removed: string[] = [];
	readCount = 0;
	png = fullHdPng();

	async createPrivateDir(): Promise<string> {
		return "/private/omp-proteinview-live-test";
	}

	async removePrivateDir(directory: string): Promise<void> {
		this.removed.push(directory);
	}

	async readFile(): Promise<Buffer> {
		this.readCount += 1;
		return this.png;
	}

	spawn(argv: string[]): PanelTransport {
		this.argv = argv;
		const outputIndex = argv.indexOf("--output");
		const outputPath = argv[outputIndex + 1] as string;
		this.transport = new FakePanelProcess(outputPath);
		queueMicrotask(() => {
			void this.transport?.ready();
		});
		return this.transport;
	}
}

async function startFakeClient(runtime = new FakeRuntime()): Promise<{
	client: ProteinViewPanelClient;
	runtime: FakeRuntime;
	process: FakePanelProcess;
}> {
	const client = await ProteinViewPanelClient.start({
		cwd: "/workspace",
		request: request(),
		runtime,
		readyTimeoutMs: 2_000,
		commandTimeoutMs: 2_000,
		shutdownTimeoutMs: 2_000,
	});
	const process = runtime.transport;
	if (process === undefined) throw new Error("fake process was not spawned");
	cleanups.push(() => client.close());
	return { client, runtime, process };
}

describe("ProteinView live panel client", () => {
	it("launches strict FullHD panel argv and reads the ready frame", async () => {
		const { client, runtime } = await startFakeClient();

		expect(runtime.argv).toEqual([
			"proteinview",
			"/workspace/input.pdb",
			"--panel-server",
			"--output",
			"/private/omp-proteinview-live-test/frame.png",
			"--panel-width",
			"1920",
			"--panel-height",
			"1080",
			"--color",
			"chain",
			"--mode",
			"cartoon",
		]);
		expect(client.snapshot?.revision).toBe(1);
		expect(client.snapshot?.frame).toMatchObject({ width: 1920, height: 1080 });
		expect(runtime.readCount).toBe(1);
	});

	it("serializes commands and reads exactly one PNG for each new revision", async () => {
		const { client, runtime, process } = await startFakeClient();
		process.autoRespond = false;

		const first = client.command({ command: "cycle_viz" });
		const second = client.command({
			command: "set_auto_rotate",
			enabled: true,
		});
		await Bun.sleep(0);
		expect(process.requests.map(entry => entry.command)).toEqual(["cycle_viz"]);

		const firstRequest = process.heldRequests.shift();
		if (firstRequest === undefined) throw new Error("first request was not held");
		await process.respond(firstRequest);
		await first;
		await Bun.sleep(0);
		expect(process.requests.map(entry => entry.command)).toEqual(["cycle_viz", "set_auto_rotate"]);

		const secondRequest = process.heldRequests.shift();
		if (secondRequest === undefined) throw new Error("second request was not held");
		await process.respond(secondRequest);
		const result = await second;
		expect(result.state.camera.auto_rotate).toBe(true);
		expect(runtime.readCount).toBe(3);
		process.autoRespond = true;
	});

	it("gracefully shuts down without killing and removes its private directory", async () => {
		const { client, runtime, process } = await startFakeClient();
		await client.close();

		expect(process.requests.at(-1)?.command).toBe("shutdown");
		expect(process.killed).toBe(false);
		expect(process.closeInputCalled).toBe(true);
		expect(runtime.removed).toEqual(["/private/omp-proteinview-live-test"]);
	});

	it("applies the shutdown timeout while the shutdown response is still pending", async () => {
		const runtime = new FakeRuntime();
		const client = await ProteinViewPanelClient.start({
			cwd: "/workspace",
			request: request(),
			runtime,
			readyTimeoutMs: 2_000,
			commandTimeoutMs: 1_000,
			shutdownTimeoutMs: 10,
		});
		const process = runtime.transport;
		if (process === undefined) throw new Error("fake process was not spawned");
		process.autoRespond = false;

		const outcome = await Promise.race([
			client.close().then(() => "closed"),
			Bun.sleep(250).then(() => "outer-timeout"),
		]);

		expect(outcome).toBe("closed");
		expect(process.requests.at(-1)?.command).toBe("shutdown");
		expect(process.killed).toBe(true);
		expect(process.closeInputCalled).toBe(true);
		expect(runtime.removed).toEqual(["/private/omp-proteinview-live-test"]);
	});

	it("preempts a pending user command when close reaches its shutdown deadline", async () => {
		const runtime = new FakeRuntime();
		const client = await ProteinViewPanelClient.start({
			cwd: "/workspace",
			request: request(),
			runtime,
			readyTimeoutMs: 2_000,
			commandTimeoutMs: 1_000,
			shutdownTimeoutMs: 10,
		});
		const process = runtime.transport;
		if (process === undefined) throw new Error("fake process was not spawned");
		process.autoRespond = false;
		const commandError = client.command({ command: "cycle_viz" }).catch(error => error);
		await Bun.sleep(0);

		const outcome = await Promise.race([
			client.close().then(() => "closed"),
			Bun.sleep(250).then(() => "outer-timeout"),
		]);
		const error = await commandError;

		expect(outcome).toBe("closed");
		expect(error).toBeInstanceOf(Error);
		expect((error as Error).message).toContain("shutdown timed out");
		expect(process.requests.map(entry => entry.command)).toEqual(["cycle_viz"]);
		expect(process.killed).toBe(true);
		expect(process.closeInputCalled).toBe(true);
		expect(runtime.removed).toEqual(["/private/omp-proteinview-live-test"]);
	});

	it("aborts promptly after ready while initial presentation is rendering", async () => {
		const { client, process } = await startFakeClient();
		process.autoRespond = false;
		const controller = new AbortController();
		const openRequest = request();
		openRequest.showLigands = false;
		openRequest.signal = controller.signal;
		const host = new ProteinViewPanelHost({
			startClient: async () => client,
		});

		const opening = host.open(openRequest, {} as ExtensionContext);
		await Bun.sleep(0);
		expect(process.requests.at(-1)?.command).toBe("set_ligands");

		controller.abort(new Error("cancel initial presentation"));
		await expect(opening).rejects.toThrow("cancel initial presentation");
		expect(process.killed).toBe(true);
		expect(client.closed).toBe(true);
	});

	it("rejects unsolicited or malformed protocol output and kills the process", async () => {
		const { client, process } = await startFakeClient();
		await process.emitRaw('{"type":"response"}\n');
		await Bun.sleep(0);

		await expect(client.command({ command: "get_state" })).rejects.toThrow("unsolicited response");
		expect(process.killed).toBe(true);
	});

	it("rejects the pending command promptly when its replacement PNG is corrupt", async () => {
		const { client, runtime, process } = await startFakeClient();
		runtime.png = Buffer.from("not a PNG", "utf8");

		await expect(client.command({ command: "cycle_viz" })).rejects.toThrow("PNG");
		expect(process.killed).toBe(true);
	});

	it("sets an exact initial interface chain without cycling through intermediate chains", async () => {
		const { client, process } = await startFakeClient();
		const host = new ProteinViewPanelHost({
			startClient: async () => client,
		});
		const openRequest = request();
		openRequest.interfaceChain = "B";
		const context = {
			ui: {
				setWidget: () => {},
				onTerminalInput: () => () => {},
				setStatus: () => {},
			},
		} as unknown as ExtensionContext;

		await host.open(openRequest, context);
		expect(process.requests.slice(0, 2).map(entry => entry.command)).toEqual(["set_chain", "set_interface"]);
		expect(process.requests[0]?.chain).toBe("B");
		await host.close();
	});

	it("deactivates and unmounts after an unexpected exit without queuing auto-rotate work", async () => {
		const { client, process } = await startFakeClient();
		await client.command({ command: "set_auto_rotate", enabled: true });
		const host = new ProteinViewPanelHost({
			startClient: async () => client,
		});
		const widgets: unknown[] = [];
		const notices: string[] = [];
		const context = {
			ui: {
				setWidget: (_key: string, widget: unknown) => widgets.push(widget),
				onTerminalInput: () => () => {},
				setStatus: () => {},
				notify: (message: string) => notices.push(message),
			},
		} as unknown as ExtensionContext;
		await host.open(request(), context);
		const requestsBeforeExit = process.requests.length;

		await process.emitStderr("render failed\nnext\tline\r");
		await process.exitUnexpectedly(9);
		await Bun.sleep(0);
		host.autoRotateTick();
		await Bun.sleep(0);

		expect(host.active).toBe(false);
		expect(client.closed).toBe(true);
		expect(widgets.at(-1)).toBeUndefined();
		expect(notices).toEqual([
			"ProteinView live panel stopped: ProteinView panel server exited unexpectedly with code 9: render failed next line",
		]);
		expect(process.requests).toHaveLength(requestsBeforeExit);
		await host.close();
	});

	it("omits mode when auto selection is requested and appends exact residue colors", () => {
		const autoRequest = request();
		autoRequest.mode = undefined;
		autoRequest.residueColors = [
			{
				chain: "A",
				residue_number: 42,
				insertion_code: "B",
				color: "FF00AA",
			},
		];
		const argv = buildPanelArgv("pv", "/tmp/frame.png", autoRequest);

		expect(argv).not.toContain("--mode");
		expect(argv.slice(-2)).toEqual(["--residue-color", "A:42[B]=FF00AA"]);
	});

	it("preserves ProteinView XYZ defaults when color and mode were not explicit", () => {
		const xyzRequest = request();
		xyzRequest.inputPath = "/private/input.xyz";
		xyzRequest.mode = undefined;
		xyzRequest.color = "structure";
		xyzRequest.explicitColor = false;
		const argv = buildPanelArgv("pv", "/tmp/frame.png", xyzRequest);

		expect(argv).not.toContain("--mode");
		expect(argv.slice(-2)).toEqual(["--color", "element"]);
	});
});

class FakeWidgetClient {
	#listeners = new Set<(snapshot: PanelSnapshot) => void>();
	snapshot: PanelSnapshot | undefined;

	constructor(snapshot: PanelSnapshot) {
		this.snapshot = snapshot;
	}

	subscribe(listener: (snapshot: PanelSnapshot) => void): () => void {
		this.#listeners.add(listener);
		listener(this.snapshot as PanelSnapshot);
		return () => {
			this.#listeners.delete(listener);
		};
	}

	publish(snapshot: PanelSnapshot): void {
		this.snapshot = snapshot;
		for (const listener of this.#listeners) listener(snapshot);
	}
}

class FakeImage implements ProteinViewPanelImage {
	setDataCalls = 0;
	lastData = "";

	setData(base64Data: string): void {
		this.setDataCalls += 1;
		this.lastData = base64Data;
	}

	render(): readonly string[] {
		return ["[image]"];
	}

	invalidate(): void {}
}

function snapshot(revision: number, state = initialState()): PanelSnapshot {
	return {
		revision,
		state,
		frame: {
			path: "/private/frame.png",
			mime_type: "image/png",
			width: FULLHD_WIDTH,
			height: FULLHD_HEIGHT,
		},
		png: fullHdPng(revision),
	};
}

describe("ProteinView live panel widget", () => {
	it("updates pixels in place and recreates the same-key image only when height constraints change", () => {
		const client = new FakeWidgetClient(snapshot(1));
		const fakeTerminal = { rows: 40 };
		const fakeTui = {
			terminal: fakeTerminal,
			imageBudget: {},
			requestRender: () => {},
		} as unknown as TUI;
		const images: Array<{ image: FakeImage; options: ImageOptions }> = [];
		const imageFactory: ProteinViewPanelImageFactory = (
			_base64Data: string,
			_mimeType: string,
			_theme: ImageTheme,
			options: ImageOptions,
			_dimensions: ImageDimensions,
		): Component & ProteinViewPanelImage => {
			const image = new FakeImage();
			images.push({ image, options });
			return image;
		};
		const widget = new ProteinViewPanelWidget({
			client,
			identifier: "4HHB",
			tui: fakeTui,
			theme: {
				fg: (_color, text) => text,
				bold: text => text,
			},
			imageFactory,
		});

		expect(widget.render(100).length).toBeGreaterThan(2);
		expect(images).toHaveLength(1);
		expect(images[0]?.options.imageKey).toBe("proteinview-live-panel");
		client.publish(snapshot(2));
		widget.render(100);
		expect(images).toHaveLength(1);
		expect(images[0]?.image.setDataCalls).toBe(1);

		fakeTerminal.rows = 50;
		widget.render(100);
		expect(images).toHaveLength(2);
		expect(images[1]?.options.imageKey).toBe("proteinview-live-panel");
		widget.dispose();
	});

	it("collapses to two full-width rows on a tiny terminal", () => {
		const client = new FakeWidgetClient(snapshot(1));
		const fakeTui = {
			terminal: { rows: 14 },
			imageBudget: {},
			requestRender: () => {},
		} as unknown as TUI;
		const widget = new ProteinViewPanelWidget({
			client,
			identifier: "4HHB",
			tui: fakeTui,
			theme: {
				fg: (_color, text) => text,
				bold: text => text,
			},
			imageFactory: () => new FakeImage(),
		});

		const lines = widget.render(60);
		expect(lines).toHaveLength(2);
		expect(lines.every(line => line.length === 60)).toBe(true);
		expect(lines[1]).toContain("larger");
		widget.dispose();
	});

	it("strips terminal escapes and control characters from structure and chain labels", () => {
		const maliciousState = initialState();
		maliciousState.presentation.current_chain_id = "\u001b]52;c;payload\u0007A\nB";
		const client = new FakeWidgetClient(snapshot(1, maliciousState));
		const fakeTui = {
			terminal: { rows: 40 },
			imageBudget: {},
			requestRender: () => {},
		} as unknown as TUI;
		const widget = new ProteinViewPanelWidget({
			client,
			identifier: "\u001b[31mevil\u001b[0m\nstructure",
			tui: fakeTui,
			theme: {
				fg: (_color, text) => text,
				bold: text => text,
			},
			imageFactory: () => new FakeImage(),
		});

		const header = widget.render(100)[0] as string;
		expect(header).not.toContain("\u001b");
		expect(header).not.toContain("\u0007");
		expect(header).not.toContain("\n");
		expect(header).toContain("evil structure");
		expect(header).toContain("chain A B");
		widget.dispose();
	});
});

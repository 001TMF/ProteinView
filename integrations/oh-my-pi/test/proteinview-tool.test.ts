import { afterEach, describe, expect, it } from "bun:test";
import { existsSync } from "node:fs";
import { mkdtemp, readFile, rm, symlink, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { deflateSync } from "node:zlib";
import type { AgentToolResult, CustomToolAPI, CustomToolContext } from "@oh-my-pi/pi-coding-agent";
import { type ProteinViewLiveOpenRequest, registerProteinViewLiveOpenHandler } from "../live/open-bridge.ts";
import {
	createProteinViewTool,
	FULLHD_HEIGHT,
	FULLHD_WIDTH,
	MAX_SOURCE_BYTES,
	type ProteinViewDetails,
	type ProteinViewParams,
	RENDER_TIMEOUT_MS,
	validateFullHdPng,
} from "../tools/index.ts";

const workspaces: string[] = [];
const liveHandlerCleanups: Array<() => void> = [];

afterEach(async () => {
	for (const cleanup of liveHandlerCleanups.splice(0)) cleanup();
	await Promise.all(workspaces.splice(0).map(path => rm(path, { recursive: true, force: true })));
});

function chainableSchema(extra: Record<string, unknown> = {}) {
	const schema = {
		...extra,
		optional: () => schema,
		describe: () => schema,
		regex: () => schema,
		min: () => schema,
		max: () => schema,
		int: () => schema,
	};
	return schema;
}

const fakeZod = {
	enum: (values: readonly string[]) => chainableSchema({ values }),
	literal: (value: string) => chainableSchema({ value }),
	string: () => chainableSchema(),
	number: () => chainableSchema(),
	boolean: () => chainableSchema(),
	array: (item: unknown) => chainableSchema({ item }),
	object: (shape: Record<string, unknown>) => chainableSchema({ shape }),
	discriminatedUnion: (key: string, options: unknown[]) => chainableSchema({ key, options }),
};

async function workspace(): Promise<string> {
	const path = await mkdtemp(join(tmpdir(), "proteinview-plugin-test-"));
	workspaces.push(path);
	return path;
}

type ExecResult = { stdout: string; stderr: string; code: number; killed: boolean };
type ExecMock = (
	command: string,
	args: string[],
	options?: { cwd?: string; signal?: AbortSignal; timeout?: number },
) => Promise<ExecResult>;

function mockApi(cwd: string, exec: ExecMock): CustomToolAPI {
	return {
		cwd,
		exec,
		zod: fakeZod,
	} as unknown as CustomToolAPI;
}

function mockContext(cwd: string, fetchImpl?: NonNullable<CustomToolContext["fetch"]>): CustomToolContext {
	return {
		sessionManager: { getCwd: () => cwd },
		fetch: fetchImpl,
	} as unknown as CustomToolContext;
}

function crc32(bytes: Buffer): number {
	let crc = 0xffffffff;
	for (const byte of bytes) {
		crc ^= byte;
		for (let bit = 0; bit < 8; bit++) {
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

function makePng(width = FULLHD_WIDTH, height = FULLHD_HEIGHT): Buffer {
	const signature = Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]);
	const ihdr = Buffer.alloc(13);
	ihdr.writeUInt32BE(width, 0);
	ihdr.writeUInt32BE(height, 4);
	ihdr[8] = 1;
	ihdr[9] = 0;
	ihdr[10] = 0;
	ihdr[11] = 0;
	ihdr[12] = 0;

	const rowBytes = Math.ceil(width / 8);
	const pixels = Buffer.alloc((rowBytes + 1) * height);
	const idat = deflateSync(pixels);
	return Buffer.concat([signature, pngChunk("IHDR", ihdr), pngChunk("IDAT", idat), pngChunk("IEND", Buffer.alloc(0))]);
}

interface ProteinViewToolExecutable {
	execute(
		toolCallId: string,
		params: ProteinViewParams,
		onUpdate: unknown,
		ctx: CustomToolContext,
		signal?: AbortSignal,
	): Promise<AgentToolResult<ProteinViewDetails>>;
}

async function execute(
	tool: ProteinViewToolExecutable,
	params: ProteinViewParams,
	ctx: CustomToolContext,
	signal?: AbortSignal,
) {
	const result = await tool.execute("call-1", params, undefined, ctx, signal);
	if (result.details === undefined) {
		throw new Error("ProteinView test tool result is missing details");
	}
	const content = result.content.map(block => {
		if (block.type !== "text") {
			throw new Error("ProteinView test tool result unexpectedly contains image content");
		}
		return block;
	});
	return { content, details: result.details };
}

describe("ProteinView Oh My Pi tool", () => {
	it("is strict, exec-approved, essential, and named for natural protein rendering requests", async () => {
		const cwd = await workspace();
		const tool = createProteinViewTool(
			mockApi(cwd, async () => ({ stdout: "", stderr: "", code: 0, killed: false })),
		);

		expect(tool.name).toBe("render_protein");
		expect(tool.label).toBe("ProteinView");
		expect(tool.strict).toBe(true);
		expect(tool.approval).toBe("exec");
		expect(tool.loadMode).toBe("essential");
		expect(tool.description).toContain("show me 1UBQ PDB protein");
	});

	it("fetches a normalized PDB ID and returns only a details image using exact snapshot argv", async () => {
		const cwd = await workspace();
		const controller = new AbortController();
		let fetchedUrl = "";
		let fetchSignal: AbortSignal | undefined;
		let captured:
			| {
					command: string;
					args: string[];
					options?: { cwd?: string; signal?: AbortSignal; timeout?: number };
					input: Buffer;
			  }
			| undefined;
		let privateDir = "";

		const fetchImpl: NonNullable<CustomToolContext["fetch"]> = async (input, init) => {
			fetchedUrl = input.toString();
			fetchSignal = init?.signal as AbortSignal;
			return new Response("data_1UBQ\n#\n", {
				status: 200,
				headers: { "content-type": "chemical/x-mmcif", "content-length": "14" },
			});
		};
		const api = mockApi(cwd, async (command, args, options) => {
			privateDir = dirname(args[2] as string);
			const input = await readFile(args[0] as string);
			await writeFile(args[2] as string, makePng());
			captured = { command, args, options, input };
			return { stdout: "", stderr: "Rendered FullHD PNG", code: 0, killed: false };
		});
		const tool = createProteinViewTool(api);

		const result = await execute(
			tool,
			{
				source: "pdb",
				pdb_id: "1ubq",
				mode: "cartoon",
				interface_chain: "A",
				show_interactions: true,
				show_ligands: false,
				residue_colors: [
					{ chain: "A", residue_number: 42, insertion_code: "a", color: "ff0000" },
					{ chain: "B", residue_number: -1, color: "00ff00" },
				],
			},
			mockContext(cwd, fetchImpl),
			controller.signal,
		);

		expect(fetchedUrl).toBe("https://files.rcsb.org/download/1UBQ.cif");
		expect(fetchSignal).toBeInstanceOf(AbortSignal);
		expect(captured?.command).toBe("proteinview");
		const inputPath = captured?.args[0] as string;
		const outputPath = captured?.args[2] as string;
		expect(captured?.args).toEqual([
			inputPath,
			"--snapshot",
			outputPath,
			"--snapshot-width",
			"1920",
			"--snapshot-height",
			"1080",
			"--mode",
			"cartoon",
			"--snapshot-interface-chain",
			"A",
			"--snapshot-interactions",
			"--snapshot-hide-ligands",
			"--residue-color",
			"A:42[A]=FF0000",
			"--residue-color",
			"B:-1=00FF00",
		]);
		expect(captured?.options).toEqual({
			cwd,
			signal: controller.signal,
			timeout: RENDER_TIMEOUT_MS,
		});
		expect(captured?.input.toString()).toBe("data_1UBQ\n#\n");
		expect(result.content).toEqual([
			{
				type: "text",
				text: "Rendered PDB 1UBQ with ProteinView at 1920x1080 (cartoon, interface, interface chain A with interactions, ligands hidden, 2 exact residue colors).",
			},
		]);
		expect(result.details.images).toHaveLength(1);
		expect(result.details.images[0]?.mimeType).toBe("image/png");
		expect(result.details.interfaceChain).toBe("A");
		expect(result.details.color).toBe("interface");
		expect(result.details.showInteractions).toBe(true);
		expect(result.details.showLigands).toBe(false);
		expect(result.details.live).toBe(false);
		expect(result.details.residueColors).toEqual([
			{ chain: "A", residue_number: 42, insertion_code: "A", color: "FF0000" },
			{ chain: "B", residue_number: -1, color: "00FF00" },
		]);
		expect(() => validateFullHdPng(Buffer.from(result.details.images[0]?.data ?? "", "base64"))).not.toThrow();
		expect(existsSync(privateDir)).toBe(false);
	});

	it("opens the live extension before removing its private input without rendering a transcript snapshot", async () => {
		const cwd = await workspace();
		const source = join(cwd, "input.pdb");
		await writeFile(source, "HEADER LIVE\nATOM\n");
		let liveRequest: ProteinViewLiveOpenRequest | undefined;
		let inputExistedDuringOpen = false;
		let execCalls = 0;
		liveHandlerCleanups.push(
			registerProteinViewLiveOpenHandler(async request => {
				liveRequest = request;
				inputExistedDuringOpen = existsSync(request.inputPath);
			}),
		);
		const tool = createProteinViewTool(
			mockApi(cwd, async () => {
				execCalls += 1;
				throw new Error("snapshot execution is forbidden after a successful live open");
			}),
		);

		const result = await execute(tool, { source: "file", path: "input.pdb", color: "chain" }, mockContext(cwd));

		expect(inputExistedDuringOpen).toBe(true);
		expect(liveRequest).toMatchObject({
			identifier: "input.pdb",
			cwd,
			mode: undefined,
			color: "chain",
			explicitColor: true,
			showInteractions: false,
			showLigands: true,
			residueColors: [],
		});
		expect(existsSync(liveRequest?.inputPath ?? "")).toBe(false);
		expect(execCalls).toBe(0);
		expect(result.details.live).toBe(true);
		expect(result.details.images).toHaveLength(0);
		expect(result.content[0]?.text).toContain("persistent ProteinView panel above the editor");
		expect(result.content[0]?.text).not.toContain("Rendered");
	});

	it("renders one FullHD snapshot when live-panel startup fails", async () => {
		const cwd = await workspace();
		const source = join(cwd, "input.pdb");
		await writeFile(source, "HEADER FALLBACK\nATOM\n");
		let execCalls = 0;
		liveHandlerCleanups.push(
			registerProteinViewLiveOpenHandler(async () => {
				throw new Error("panel unavailable");
			}),
		);
		const tool = createProteinViewTool(
			mockApi(cwd, async (_command, args) => {
				execCalls += 1;
				await writeFile(args[2] as string, makePng());
				return { stdout: "", stderr: "", code: 0, killed: false };
			}),
		);

		const result = await execute(tool, { source: "file", path: "input.pdb", color: "chain" }, mockContext(cwd));

		expect(execCalls).toBe(1);
		expect(result.details.live).toBe(false);
		expect(result.details.liveError).toBe("panel unavailable");
		expect(result.details.images).toHaveLength(1);
		expect(result.content[0]?.text).toContain(
			"Live panel unavailable (panel unavailable); showing the FullHD snapshot.",
		);
	});

	it("reports and renders ProteinView XYZ defaults consistently", async () => {
		const cwd = await workspace();
		await writeFile(join(cwd, "molecule.xyz"), "1\nmolecule\nC 0 0 0\n");
		let renderedArgs: string[] = [];
		const tool = createProteinViewTool(
			mockApi(cwd, async (_command, args) => {
				renderedArgs = args;
				await writeFile(args[2] as string, makePng());
				return { stdout: "", stderr: "", code: 0, killed: false };
			}),
		);

		const result = await execute(tool, { source: "file", path: "molecule.xyz" }, mockContext(cwd));

		expect(renderedArgs).not.toContain("--mode");
		expect(renderedArgs).not.toContain("--color");
		expect(result.details.mode).toBe("wireframe");
		expect(result.details.color).toBe("element");
		expect(result.content[0]?.text).toContain("(wireframe, element");
	});

	it("snapshots a workspace-local file through a private copy", async () => {
		const cwd = await workspace();
		const source = join(cwd, "input.pdb");
		await writeFile(source, "HEADER TEST\nATOM\n");
		let privateInput = "";
		let copiedBytes = "";
		let privateDir = "";
		let renderedArgs: string[] = [];

		const tool = createProteinViewTool(
			mockApi(cwd, async (_command, args) => {
				renderedArgs = args;
				privateInput = args[0] as string;
				privateDir = dirname(privateInput);
				copiedBytes = (await readFile(privateInput)).toString();
				await writeFile(args[2] as string, makePng());
				return { stdout: "", stderr: "", code: 0, killed: false };
			}),
		);
		const result = await execute(tool, { source: "file", path: "input.pdb", color: "chain" }, mockContext(cwd));

		expect(privateInput).not.toBe(source);
		expect(copiedBytes).toBe("HEADER TEST\nATOM\n");
		expect(result.details.identifier).toBe("input.pdb");
		expect(result.details.mode).toBe("auto");
		expect(result.details.color).toBe("chain");
		expect(result.details.showInteractions).toBe(false);
		expect(result.details.showLigands).toBe(true);
		expect(renderedArgs).not.toContain("--mode");
		expect(renderedArgs.slice(-2)).toEqual(["--color", "chain"]);
		expect(existsSync(privateDir)).toBe(false);
	});

	it("sanitizes local filenames before showing them in the live terminal UI", async () => {
		const cwd = await workspace();
		const unsafeName = "evil\u001b]0;owned\u0007\n.pdb";
		await writeFile(join(cwd, unsafeName), "HEADER TEST\n");
		let opened: ProteinViewLiveOpenRequest | undefined;
		liveHandlerCleanups.push(
			registerProteinViewLiveOpenHandler(async request => {
				opened = request;
			}),
		);
		const tool = createProteinViewTool(
			mockApi(cwd, async () => {
				throw new Error("live success must not invoke snapshot rendering");
			}),
		);

		const result = await execute(tool, { source: "file", path: unsafeName }, mockContext(cwd));

		expect(opened?.identifier).toBe("evil .pdb");
		expect(result.details.identifier).toBe("evil .pdb");
		expect(result.content[0]?.text).not.toContain("\u001b");
		expect(result.content[0]?.text).not.toContain("\n");
	});

	it("rejects malformed PDB IDs before fetch or process execution", async () => {
		const cwd = await workspace();
		let fetched = false;
		let executed = false;
		const tool = createProteinViewTool(
			mockApi(cwd, async () => {
				executed = true;
				return { stdout: "", stderr: "", code: 0, killed: false };
			}),
		);

		await expect(
			execute(
				tool,
				{ source: "pdb", pdb_id: "1UBQ;touch /tmp/nope" },
				mockContext(cwd, async () => {
					fetched = true;
					return new Response();
				}),
			),
		).rejects.toThrow("PDB ID must be exactly four");
		expect(fetched).toBe(false);
		expect(executed).toBe(false);
	});

	it("rejects malformed and duplicate exact residue colors before execution", async () => {
		const cwd = await workspace();
		const source = join(cwd, "input.pdb");
		await writeFile(source, "HEADER TEST\n");
		let executions = 0;
		const tool = createProteinViewTool(
			mockApi(cwd, async () => {
				executions += 1;
				return { stdout: "", stderr: "", code: 0, killed: false };
			}),
		);

		await expect(
			execute(
				tool,
				{
					source: "file",
					path: "input.pdb",
					residue_colors: [{ chain: "A", residue_number: 1, color: "#FF0000" }],
				},
				mockContext(cwd),
			),
		).rejects.toThrow("six hexadecimal digits");
		await expect(
			execute(
				tool,
				{
					source: "file",
					path: "input.pdb",
					residue_colors: [
						{ chain: "A", residue_number: 1, insertion_code: "a", color: "FF0000" },
						{ chain: "A", residue_number: 1, insertion_code: "A", color: "00FF00" },
					],
				},
				mockContext(cwd),
			),
		).rejects.toThrow("duplicate target");
		expect(executions).toBe(0);
	});

	it("rejects paths outside the workspace and symlinks that escape it", async () => {
		const cwd = await workspace();
		const outside = await workspace();
		const outsideFile = join(outside, "outside.pdb");
		await writeFile(outsideFile, "HEADER OUTSIDE\n");
		await symlink(outsideFile, join(cwd, "escape.pdb"));
		let executions = 0;
		const tool = createProteinViewTool(
			mockApi(cwd, async () => {
				executions += 1;
				return { stdout: "", stderr: "", code: 0, killed: false };
			}),
		);

		await expect(execute(tool, { source: "file", path: outsideFile }, mockContext(cwd))).rejects.toThrow(
			"inside the current workspace",
		);
		await expect(execute(tool, { source: "file", path: "escape.pdb" }, mockContext(cwd))).rejects.toThrow(
			"inside the current workspace",
		);
		expect(executions).toBe(0);
	});

	it("rejects an oversized declared RCSB response before execution", async () => {
		const cwd = await workspace();
		let executions = 0;
		const tool = createProteinViewTool(
			mockApi(cwd, async () => {
				executions += 1;
				return { stdout: "", stderr: "", code: 0, killed: false };
			}),
		);

		await expect(
			execute(
				tool,
				{ source: "pdb", pdb_id: "1UBQ" },
				mockContext(
					cwd,
					async () =>
						new Response("small", {
							headers: { "content-length": String(MAX_SOURCE_BYTES + 1) },
						}),
				),
			),
		).rejects.toThrow("safety limit");
		expect(executions).toBe(0);
	});

	it("propagates cancellation through api.exec and removes temporary files", async () => {
		const cwd = await workspace();
		const source = join(cwd, "input.pdb");
		await writeFile(source, "HEADER TEST\n");
		const controller = new AbortController();
		let privateDir = "";
		let forwardedSignal: AbortSignal | undefined;
		const cancellation = new Error("cancelled by test");
		const tool = createProteinViewTool(
			mockApi(cwd, async (_command, args, options) => {
				privateDir = dirname(args[2] as string);
				forwardedSignal = options?.signal;
				controller.abort(cancellation);
				return { stdout: "", stderr: "", code: 0, killed: true };
			}),
		);

		await expect(
			execute(tool, { source: "file", path: "input.pdb" }, mockContext(cwd), controller.signal),
		).rejects.toThrow("cancelled by test");
		expect(forwardedSignal).toBe(controller.signal);
		expect(existsSync(privateDir)).toBe(false);
	});

	it("rejects wrong PNG dimensions and removes the private directory", async () => {
		const cwd = await workspace();
		const source = join(cwd, "input.pdb");
		await writeFile(source, "HEADER TEST\n");
		let privateDir = "";
		const tool = createProteinViewTool(
			mockApi(cwd, async (_command, args) => {
				privateDir = dirname(args[2] as string);
				await writeFile(args[2] as string, makePng(800, 600));
				return { stdout: "", stderr: "", code: 0, killed: false };
			}),
		);

		await expect(execute(tool, { source: "file", path: "input.pdb" }, mockContext(cwd))).rejects.toThrow(
			"expected 1920x1080",
		);
		expect(existsSync(privateDir)).toBe(false);
	});

	it("rejects a truncated or checksum-corrupt PNG", () => {
		const valid = makePng();
		expect(() => validateFullHdPng(valid.subarray(0, 33))).toThrow("valid PNG signature");
		const corrupt = Buffer.from(valid);
		corrupt[29] ^= 0xff;
		expect(() => validateFullHdPng(corrupt)).toThrow("invalid checksum");
	});

	it("reports render timeouts and missing snapshot support without returning an image", async () => {
		const cwd = await workspace();
		const source = join(cwd, "input.pdb");
		await writeFile(source, "HEADER TEST\n");
		const timedOutTool = createProteinViewTool(
			mockApi(cwd, async () => ({ stdout: "", stderr: "", code: 0, killed: true })),
		);
		await expect(execute(timedOutTool, { source: "file", path: "input.pdb" }, mockContext(cwd))).rejects.toThrow(
			"exceeded 60 seconds",
		);

		const oldBinaryTool = createProteinViewTool(
			mockApi(cwd, async () => ({
				stdout: "",
				stderr: "unexpected argument '--snapshot'",
				code: 2,
				killed: false,
			})),
		);
		await expect(execute(oldBinaryTool, { source: "file", path: "input.pdb" }, mockContext(cwd))).rejects.toThrow(
			"supports --snapshot",
		);
	});
});

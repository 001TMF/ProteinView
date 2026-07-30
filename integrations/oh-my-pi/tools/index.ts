import { constants as fsConstants } from "node:fs";
import type { FileHandle } from "node:fs/promises";
import { mkdtemp, open, realpath, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { basename, extname, isAbsolute, join, relative, resolve, sep } from "node:path";
import type { CustomToolAPI, CustomToolContext, CustomToolFactory, ExecResult } from "@oh-my-pi/pi-coding-agent";
import { requestProteinViewLiveOpen } from "../live/open-bridge.ts";

export const FULLHD_WIDTH = 1920;
export const FULLHD_HEIGHT = 1080;
export const MAX_SOURCE_BYTES = 64 * 1024 * 1024;
export const MAX_PNG_BYTES = 24 * 1024 * 1024;
export const FETCH_TIMEOUT_MS = 30_000;
export const RENDER_TIMEOUT_MS = 60_000;
export const MAX_RESIDUE_COLORS = 256;

const PDB_ID_PATTERN = /^[1-9][A-Za-z0-9]{3}$/;
const ALLOWED_EXTENSIONS = new Set([".pdb", ".cif", ".mmcif", ".xyz"]);
const MODES = ["cartoon", "backbone", "wireframe"] as const;
const COLORS = ["structure", "element", "chain", "plddt", "bfactor", "rainbow"] as const;
const PNG_SIGNATURE = Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]);
const CRC32_TABLE = Uint32Array.from({ length: 256 }, (_, index) => {
	let value = index;
	for (let bit = 0; bit < 8; bit++) {
		value = (value >>> 1) ^ (value & 1 ? 0xedb88320 : 0);
	}
	return value >>> 0;
});
// biome-ignore lint/complexity/useRegexLiterals: Constructor strings avoid literal control escapes rejected by noControlCharactersInRegex.
const ANSI_ESCAPE = new RegExp("\\x1b(?:\\[[0-?]*[ -/]*[@-~]|\\][^\\x07]*(?:\\x07|\\x1b\\\\))", "g");
// biome-ignore lint/complexity/useRegexLiterals: Constructor strings avoid literal control escapes rejected by noControlCharactersInRegex.
const CONTROL_CHARACTERS = new RegExp("[\\u0000-\\u0008\\u000b\\u000c\\u000e-\\u001f\\u007f]", "g");
// Includes line breaks, tabs, DEL, and C1 controls because this text is later
// rendered inside a terminal widget rather than only shown as process output.
// biome-ignore lint/complexity/useRegexLiterals: Constructor strings avoid literal control escapes rejected by noControlCharactersInRegex.
const DISPLAY_CONTROL_CHARACTERS = new RegExp("[\\u0000-\\u001f\\u007f-\\u009f]", "g");

type RenderMode = (typeof MODES)[number];
type ColorScheme = (typeof COLORS)[number];

export interface ProteinViewResidueColor {
	chain: string;
	residue_number: number;
	insertion_code?: string;
	color: string;
}

export interface PdbParams {
	source: "pdb";
	pdb_id: string;
	mode?: RenderMode;
	color?: ColorScheme;
	interface_chain?: string;
	show_interactions?: boolean;
	show_ligands?: boolean;
	residue_colors?: ProteinViewResidueColor[];
}

export interface FileParams {
	source: "file";
	path: string;
	mode?: RenderMode;
	color?: ColorScheme;
	interface_chain?: string;
	show_interactions?: boolean;
	show_ligands?: boolean;
	residue_colors?: ProteinViewResidueColor[];
}

export type ProteinViewParams = PdbParams | FileParams;

interface InlineImage {
	data: string;
	mimeType: "image/png";
}

export interface ProteinViewDetails {
	source: "pdb" | "file";
	identifier: string;
	width: number;
	height: number;
	mode: RenderMode | "auto";
	color: ColorScheme | "interface";
	interfaceChain?: string;
	showInteractions: boolean;
	showLigands: boolean;
	residueColors: ProteinViewResidueColor[];
	live: boolean;
	liveError?: string;
	images: InlineImage[];
}

interface DeadlineSignal {
	signal: AbortSignal;
	didTimeout(): boolean;
	dispose(): void;
}

function abortError(message = "Protein rendering cancelled"): Error {
	const error = new Error(message);
	error.name = "AbortError";
	return error;
}

function throwIfAborted(signal?: AbortSignal): void {
	if (!signal?.aborted) return;
	throw signal.reason instanceof Error ? signal.reason : abortError();
}

function createDeadlineSignal(parent: AbortSignal | undefined, timeoutMs: number): DeadlineSignal {
	const controller = new AbortController();
	let timedOut = false;

	const forwardAbort = () => {
		controller.abort(parent?.reason instanceof Error ? parent.reason : abortError());
	};
	if (parent?.aborted) {
		forwardAbort();
	} else {
		parent?.addEventListener("abort", forwardAbort, { once: true });
	}

	const timer = setTimeout(() => {
		timedOut = true;
		controller.abort(abortError(`RCSB fetch exceeded ${Math.ceil(timeoutMs / 1000)} seconds`));
	}, timeoutMs);

	return {
		signal: controller.signal,
		didTimeout: () => timedOut,
		dispose: () => {
			clearTimeout(timer);
			parent?.removeEventListener("abort", forwardAbort);
		},
	};
}

function normalizePdbId(value: unknown): string {
	if (typeof value !== "string") {
		throw new Error("PDB ID must be exactly four alphanumeric characters and begin with 1-9");
	}
	const candidate = value.trim();
	if (!PDB_ID_PATTERN.test(candidate)) {
		throw new Error("PDB ID must be exactly four alphanumeric characters and begin with 1-9");
	}
	return candidate.toUpperCase();
}

function normalizeChoice<T extends string>(value: unknown, allowed: readonly T[], fallback: T, label: string): T {
	if (value === undefined) return fallback;
	if (typeof value === "string" && allowed.includes(value as T)) return value as T;
	throw new Error(`Unsupported ProteinView ${label}: ${String(value)}`);
}

function normalizeInterfaceChain(value: unknown): string | undefined {
	if (value === undefined) return undefined;
	if (typeof value !== "string") throw new Error("ProteinView interface chain must be a string");
	const candidate = value.trim();
	const hasControlCharacter = [...candidate].some(character => {
		const codePoint = character.codePointAt(0) ?? 0;
		return codePoint <= 0x1f || codePoint === 0x7f;
	});
	if (!candidate || candidate.length > 32 || hasControlCharacter) {
		throw new Error("ProteinView interface chain must contain 1-32 printable characters");
	}
	return candidate;
}

function normalizeResidueIdentifier(value: unknown, label: string): string {
	if (typeof value !== "string") throw new Error(`ProteinView residue ${label} must be a string`);
	const candidate = value.trim();
	if (
		!candidate ||
		Buffer.byteLength(candidate, "utf8") > 32 ||
		![...candidate].every(character => {
			const codePoint = character.codePointAt(0) ?? 0;
			return codePoint >= 0x21 && codePoint <= 0x7e;
		})
	) {
		throw new Error(`ProteinView residue ${label} must contain 1-32 printable non-whitespace ASCII bytes`);
	}
	return candidate;
}

function normalizeResidueColors(value: unknown): ProteinViewResidueColor[] {
	if (value === undefined) return [];
	if (!Array.isArray(value)) throw new Error("ProteinView residue_colors must be an array");
	if (value.length > MAX_RESIDUE_COLORS) {
		throw new Error(`ProteinView supports at most ${MAX_RESIDUE_COLORS} exact residue colors`);
	}

	const seen = new Set<string>();
	return value.map((entry, index) => {
		if (typeof entry !== "object" || entry === null || Array.isArray(entry)) {
			throw new Error(`ProteinView residue_colors[${index}] must be an object`);
		}
		const record = entry as Record<string, unknown>;
		const chain = normalizeResidueIdentifier(record.chain, "chain");
		const residueNumber = record.residue_number;
		if (
			typeof residueNumber !== "number" ||
			!Number.isInteger(residueNumber) ||
			residueNumber < -2_147_483_648 ||
			residueNumber > 2_147_483_647
		) {
			throw new Error(`ProteinView residue_colors[${index}].residue_number must be a signed 32-bit integer`);
		}
		const insertionCode =
			record.insertion_code === undefined
				? undefined
				: normalizeResidueIdentifier(record.insertion_code, "insertion code").toUpperCase();
		if (typeof record.color !== "string" || !/^[0-9A-Fa-f]{6}$/.test(record.color)) {
			throw new Error(`ProteinView residue_colors[${index}].color must be exactly six hexadecimal digits`);
		}
		const duplicateKey = JSON.stringify([chain, residueNumber, insertionCode ?? null]);
		if (seen.has(duplicateKey)) {
			throw new Error(
				`ProteinView residue_colors contains duplicate target ${chain}:${residueNumber}${insertionCode ?? ""}`,
			);
		}
		seen.add(duplicateKey);
		return {
			chain,
			residue_number: residueNumber,
			...(insertionCode === undefined ? {} : { insertion_code: insertionCode }),
			color: record.color.toUpperCase(),
		};
	});
}

function residueColorArg(residue: ProteinViewResidueColor): string {
	const insertion = residue.insertion_code ? `[${residue.insertion_code}]` : "";
	return `${residue.chain}:${residue.residue_number}${insertion}=${residue.color}`;
}

function normalizeBoolean(value: unknown, fallback: boolean, label: string): boolean {
	if (value === undefined) return fallback;
	if (typeof value === "boolean") return value;
	throw new Error(`ProteinView ${label} must be true or false`);
}

function isWithinWorkspace(workspace: string, candidate: string): boolean {
	const rel = relative(workspace, candidate);
	return rel !== "" && rel !== ".." && !rel.startsWith(`..${sep}`) && !isAbsolute(rel);
}

async function readWorkspaceStructure(
	cwd: string,
	requestedPath: string,
): Promise<{
	bytes: Buffer;
	extension: string;
	identifier: string;
}> {
	if (typeof requestedPath !== "string" || !requestedPath.trim()) {
		throw new Error("Protein structure path cannot be empty");
	}

	const workspace = await realpath(cwd);
	const requested = resolve(workspace, requestedPath);
	let canonical: string;
	try {
		canonical = await realpath(requested);
	} catch {
		throw new Error(`Protein structure file does not exist: ${requestedPath}`);
	}

	if (!isWithinWorkspace(workspace, canonical)) {
		throw new Error("Protein structure paths must resolve to a regular file inside the current workspace");
	}

	const extension = extname(canonical).toLowerCase();
	if (!ALLOWED_EXTENSIONS.has(extension)) {
		throw new Error("Protein structure must use a .pdb, .cif, .mmcif, or .xyz extension");
	}

	const noFollow = typeof fsConstants.O_NOFOLLOW === "number" ? fsConstants.O_NOFOLLOW : 0;
	let handle: FileHandle;
	try {
		handle = await open(canonical, fsConstants.O_RDONLY | noFollow);
	} catch {
		throw new Error(`Unable to open protein structure file: ${requestedPath}`);
	}

	try {
		const fileStat = await handle.stat();
		if (!fileStat.isFile()) {
			throw new Error("Protein structure path must identify a regular file");
		}
		if (fileStat.size <= 0) throw new Error("Protein structure file is empty");
		if (fileStat.size > MAX_SOURCE_BYTES) {
			throw new Error(`Protein structure exceeds the ${MAX_SOURCE_BYTES}-byte safety limit`);
		}

		const bytes = await handle.readFile();
		if (bytes.length <= 0) throw new Error("Protein structure file is empty");
		if (bytes.length > MAX_SOURCE_BYTES) {
			throw new Error(`Protein structure exceeds the ${MAX_SOURCE_BYTES}-byte safety limit`);
		}
		return { bytes, extension, identifier: basename(canonical) };
	} finally {
		await handle.close();
	}
}

async function readResponseBody(response: Response, signal: AbortSignal): Promise<Buffer> {
	const declaredLength = response.headers.get("content-length");
	if (declaredLength) {
		const parsedLength = Number(declaredLength);
		if (Number.isFinite(parsedLength) && parsedLength > MAX_SOURCE_BYTES) {
			throw new Error(`RCSB structure exceeds the ${MAX_SOURCE_BYTES}-byte safety limit`);
		}
	}
	if (!response.body) throw new Error("RCSB returned an empty response body");

	const reader = response.body.getReader();
	const chunks: Uint8Array[] = [];
	let total = 0;
	try {
		while (true) {
			throwIfAborted(signal);
			const { done, value } = await reader.read();
			if (done) break;
			if (!value) continue;
			total += value.byteLength;
			if (total > MAX_SOURCE_BYTES) {
				await reader.cancel("Protein structure response exceeded safety limit");
				throw new Error(`RCSB structure exceeds the ${MAX_SOURCE_BYTES}-byte safety limit`);
			}
			chunks.push(value);
		}
	} finally {
		reader.releaseLock();
	}

	if (total === 0) throw new Error("RCSB returned an empty structure");
	const bytes = Buffer.allocUnsafe(total);
	let offset = 0;
	for (const chunk of chunks) {
		bytes.set(chunk, offset);
		offset += chunk.byteLength;
	}
	return bytes;
}

async function fetchPdb(
	fetchImpl: NonNullable<CustomToolContext["fetch"]>,
	pdbId: string,
	parentSignal?: AbortSignal,
): Promise<Buffer> {
	const deadline = createDeadlineSignal(parentSignal, FETCH_TIMEOUT_MS);
	try {
		const url = `https://files.rcsb.org/download/${pdbId}.cif`;
		let response: Response;
		try {
			response = await fetchImpl(url, {
				method: "GET",
				redirect: "error",
				signal: deadline.signal,
				headers: { accept: "chemical/x-mmcif,text/plain;q=0.9" },
			});
		} catch (error) {
			if (parentSignal?.aborted) throwIfAborted(parentSignal);
			if (deadline.didTimeout()) {
				throw new Error(`RCSB fetch for ${pdbId} exceeded ${FETCH_TIMEOUT_MS / 1000} seconds`);
			}
			throw new Error(
				`Unable to fetch PDB ${pdbId} from RCSB: ${error instanceof Error ? error.message : String(error)}`,
			);
		}

		if (!response.ok) {
			throw new Error(`Unable to fetch PDB ${pdbId} from RCSB: HTTP ${response.status}`);
		}
		const contentType = response.headers.get("content-type")?.toLowerCase() ?? "";
		if (contentType.includes("text/html")) {
			throw new Error(`RCSB returned HTML instead of structure data for PDB ${pdbId}`);
		}
		return await readResponseBody(response, deadline.signal);
	} finally {
		deadline.dispose();
	}
}

function sanitizeProcessOutput(value: string): string {
	return value.replace(ANSI_ESCAPE, "").replace(CONTROL_CHARACTERS, " ").trim().slice(0, 4000);
}

function sanitizeDisplayIdentifier(value: string): string {
	const cleaned = value.replace(ANSI_ESCAPE, "").replace(DISPLAY_CONTROL_CHARACTERS, " ").trim().slice(0, 120);
	return cleaned || "structure";
}

function crc32(bytes: Buffer): number {
	let crc = 0xffffffff;
	for (const byte of bytes) {
		crc = (crc >>> 8) ^ (CRC32_TABLE[(crc ^ byte) & 0xff] as number);
	}
	return (crc ^ 0xffffffff) >>> 0;
}

export function validateFullHdPng(bytes: Buffer): void {
	if (bytes.length > MAX_PNG_BYTES) {
		throw new Error(`ProteinView PNG exceeds the ${MAX_PNG_BYTES}-byte safety limit`);
	}
	if (bytes.length < 57 || !bytes.subarray(0, PNG_SIGNATURE.length).equals(PNG_SIGNATURE)) {
		throw new Error("ProteinView did not produce a valid PNG signature");
	}

	let offset = PNG_SIGNATURE.length;
	let width = 0;
	let height = 0;
	let chunkIndex = 0;
	let sawImageData = false;
	let sawEnd = false;
	while (offset + 12 <= bytes.length) {
		const length = bytes.readUInt32BE(offset);
		const chunkEnd = offset + 12 + length;
		if (chunkEnd > bytes.length) {
			throw new Error("ProteinView PNG contains a truncated chunk");
		}
		const typeStart = offset + 4;
		const dataStart = offset + 8;
		const dataEnd = dataStart + length;
		const chunkType = bytes.toString("ascii", typeStart, dataStart);
		const expectedCrc = bytes.readUInt32BE(dataEnd);
		if (crc32(bytes.subarray(typeStart, dataEnd)) !== expectedCrc) {
			throw new Error("ProteinView PNG contains an invalid checksum");
		}

		if (chunkIndex === 0) {
			if (chunkType !== "IHDR" || length !== 13) {
				throw new Error("ProteinView PNG is missing a valid IHDR chunk");
			}
			width = bytes.readUInt32BE(dataStart);
			height = bytes.readUInt32BE(dataStart + 4);
		} else if (chunkType === "IDAT") {
			sawImageData ||= length > 0;
		} else if (chunkType === "IEND") {
			if (length !== 0 || chunkEnd !== bytes.length) {
				throw new Error("ProteinView PNG contains an invalid ending");
			}
			sawEnd = true;
			offset = chunkEnd;
			break;
		}
		offset = chunkEnd;
		chunkIndex += 1;
	}

	if (!sawImageData || !sawEnd || offset !== bytes.length) {
		throw new Error("ProteinView PNG is missing complete image data");
	}
	if (width !== FULLHD_WIDTH || height !== FULLHD_HEIGHT) {
		throw new Error(`ProteinView PNG dimensions were ${width}x${height}; expected ${FULLHD_WIDTH}x${FULLHD_HEIGHT}`);
	}
}

async function readAndValidatePng(outputPath: string): Promise<Buffer> {
	const noFollow = typeof fsConstants.O_NOFOLLOW === "number" ? fsConstants.O_NOFOLLOW : 0;
	let handle: FileHandle;
	try {
		handle = await open(outputPath, fsConstants.O_RDONLY | noFollow);
	} catch {
		throw new Error("ProteinView exited successfully but did not create a snapshot");
	}

	try {
		const outputStat = await handle.stat();
		if (!outputStat.isFile()) throw new Error("ProteinView snapshot output is not a regular file");
		if (outputStat.size <= 0) throw new Error("ProteinView produced an empty PNG");
		if (outputStat.size > MAX_PNG_BYTES) {
			throw new Error(`ProteinView PNG exceeds the ${MAX_PNG_BYTES}-byte safety limit`);
		}
		const bytes = await handle.readFile();
		validateFullHdPng(bytes);
		return bytes;
	} finally {
		await handle.close();
	}
}

function missingBinaryError(): Error {
	return new Error(
		"ProteinView binary not found. Install the live-panel/snapshot-capable ProteinView build or set PROTEINVIEW_BIN to its executable path.",
	);
}

export function createProteinViewTool(api: CustomToolAPI) {
	const z = api.zod;
	const styleFields = {
		mode: z.enum(MODES).optional().describe("Protein visualization style; omit to auto-select for large structures"),
		color: z.enum(COLORS).optional().describe("Protein color scheme; defaults to structure"),
		interface_chain: z.string().min(1).max(32).optional().describe("Focus-chain ID for interface-residue coloring"),
		show_interactions: z
			.boolean()
			.optional()
			.describe("Overlay classified interface interaction lines; requires interface_chain"),
		show_ligands: z.boolean().optional().describe("Show ligands and ions; defaults to true"),
		residue_colors: z
			.array(
				z.object({
					chain: z.string().min(1).max(32).describe("Exact case-sensitive polymer chain ID"),
					residue_number: z.number().int().describe("Exact signed PDB/mmCIF residue sequence number"),
					insertion_code: z
						.string()
						.min(1)
						.max(32)
						.optional()
						.describe("Exact residue insertion code; omit only for a blank insertion code"),
					color: z
						.string()
						.regex(/^[0-9A-Fa-f]{6}$/)
						.describe("Exact RGB color as six hexadecimal digits, for example FF0000"),
				}),
			)
			.max(MAX_RESIDUE_COLORS)
			.optional()
			.describe(
				"Exact polymer residues to color; these overrides take precedence over every regular and interface scheme",
			),
	};
	const parameters = z.discriminatedUnion("source", [
		z.object({
			source: z.literal("pdb"),
			pdb_id: z
				.string()
				.regex(PDB_ID_PATTERN)
				.describe("Four-character RCSB Protein Data Bank identifier, for example 1UBQ"),
			...styleFields,
		}),
		z.object({
			source: z.literal("file"),
			path: z.string().min(1).describe("Path to a .pdb, .cif, .mmcif, or .xyz file inside the current workspace"),
			...styleFields,
		}),
	]);

	return {
		name: "render_protein",
		label: "ProteinView",
		strict: true,
		loadMode: "essential" as const,
		approval: "exec" as const,
		description:
			'Open a 3D protein structure in ProteinView\'s persistent FullHD live panel, with an inline FullHD snapshot fallback when the panel is unavailable. Use this whenever the user asks to show, view, visualize, render, or inspect a PDB/Protein Data Bank ID or a local PDB, mmCIF, CIF, or XYZ structure. For example, “show me 1UBQ PDB protein” means source="pdb" and pdb_id="1UBQ".',
		parameters,
		async execute(
			_toolCallId: string,
			rawParams: ProteinViewParams,
			_onUpdate: unknown,
			ctx: CustomToolContext,
			signal?: AbortSignal,
		) {
			throwIfAborted(signal);
			const mode =
				rawParams.mode === undefined ? undefined : normalizeChoice(rawParams.mode, MODES, "cartoon", "mode");
			const color = normalizeChoice(rawParams.color, COLORS, "structure", "color scheme");
			const interfaceChain = normalizeInterfaceChain(rawParams.interface_chain);
			const showInteractions = normalizeBoolean(rawParams.show_interactions, false, "show_interactions");
			const showLigands = normalizeBoolean(rawParams.show_ligands, true, "show_ligands");
			const residueColors = normalizeResidueColors(rawParams.residue_colors);
			if (showInteractions && interfaceChain === undefined) {
				throw new Error("ProteinView show_interactions requires interface_chain");
			}
			if (interfaceChain !== undefined && rawParams.color !== undefined) {
				throw new Error("ProteinView interface_chain uses the interface palette and cannot use color");
			}
			let pdbId: string | undefined;
			if (rawParams.source === "pdb") {
				pdbId = normalizePdbId(rawParams.pdb_id);
			} else if (rawParams.source !== "file") {
				throw new Error(`Unsupported protein source: ${String((rawParams as { source?: unknown }).source)}`);
			}

			const cwd = ctx.sessionManager.getCwd() || api.cwd;
			const privateDir = await mkdtemp(join(tmpdir(), "omp-proteinview-"));
			try {
				throwIfAborted(signal);
				let sourceBytes: Buffer;
				let extension: string;
				let identifier: string;

				if (rawParams.source === "pdb") {
					const canonicalId = pdbId as string;
					sourceBytes = await fetchPdb(ctx.fetch ?? fetch, canonicalId, signal);
					extension = ".cif";
					identifier = canonicalId;
				} else {
					const local = await readWorkspaceStructure(cwd, rawParams.path);
					sourceBytes = local.bytes;
					extension = local.extension;
					identifier = sanitizeDisplayIdentifier(local.identifier);
				}
				throwIfAborted(signal);
				const usesXyzDefaults = extension === ".xyz";
				const resolvedColor = usesXyzDefaults && rawParams.color === undefined ? "element" : color;
				const resolvedMode = usesXyzDefaults && mode === undefined ? "wireframe" : (mode ?? "auto");
				const effectiveColor = interfaceChain === undefined ? resolvedColor : "interface";

				const privateInput = join(privateDir, `structure${extension}`);
				await writeFile(privateInput, sourceBytes, { flag: "wx", mode: 0o600 });
				throwIfAborted(signal);

				const command = process.env.PROTEINVIEW_BIN?.trim() || "proteinview";
				let live = false;
				let liveError: string | undefined;
				try {
					live = await requestProteinViewLiveOpen({
						inputPath: privateInput,
						identifier,
						cwd,
						binary: command,
						signal,
						mode,
						color: resolvedColor,
						explicitColor: rawParams.color !== undefined,
						interfaceChain,
						showInteractions,
						showLigands,
						residueColors,
					});
				} catch (error) {
					if (signal?.aborted) throwIfAborted(signal);
					liveError =
						sanitizeProcessOutput(error instanceof Error ? error.message : String(error)) ||
						"unknown live-panel startup error";
				}

				const images: InlineImage[] = [];
				if (!live) {
					const privateOutput = join(privateDir, "proteinview-fullhd.png");
					const args = [
						privateInput,
						"--snapshot",
						privateOutput,
						"--snapshot-width",
						String(FULLHD_WIDTH),
						"--snapshot-height",
						String(FULLHD_HEIGHT),
					];
					if (mode !== undefined) args.push("--mode", mode);
					if (rawParams.color !== undefined) args.push("--color", color);
					if (interfaceChain !== undefined) {
						args.push("--snapshot-interface-chain", interfaceChain);
					}
					if (showInteractions) args.push("--snapshot-interactions");
					if (!showLigands) args.push("--snapshot-hide-ligands");
					for (const residueColor of residueColors) {
						args.push("--residue-color", residueColorArg(residueColor));
					}

					let execution: ExecResult;
					try {
						execution = await api.exec(command, args, {
							cwd,
							signal,
							timeout: RENDER_TIMEOUT_MS,
						});
					} catch (error) {
						if (signal?.aborted) throwIfAborted(signal);
						if (
							typeof error === "object" &&
							error !== null &&
							"code" in error &&
							(error as { code?: unknown }).code === "ENOENT"
						) {
							throw missingBinaryError();
						}
						throw new Error(
							`Unable to start ProteinView: ${error instanceof Error ? error.message : String(error)}`,
						);
					}

					if (execution.killed) {
						if (signal?.aborted) throwIfAborted(signal);
						throw new Error(`ProteinView snapshot exceeded ${RENDER_TIMEOUT_MS / 1000} seconds`);
					}
					throwIfAborted(signal);
					if (execution.code !== 0) {
						const diagnostic =
							sanitizeProcessOutput(execution.stderr) ||
							sanitizeProcessOutput(execution.stdout) ||
							`exit code ${execution.code}`;
						throw new Error(
							`ProteinView snapshot failed: ${diagnostic}. Ensure the installed ProteinView supports --snapshot.`,
						);
					}

					const png = await readAndValidatePng(privateOutput);
					throwIfAborted(signal);
					images.push({
						data: png.toString("base64"),
						mimeType: "image/png",
					});
				}

				const details: ProteinViewDetails = {
					source: rawParams.source,
					identifier,
					width: FULLHD_WIDTH,
					height: FULLHD_HEIGHT,
					mode: resolvedMode,
					color: effectiveColor,
					interfaceChain,
					showInteractions,
					showLigands,
					residueColors,
					live,
					...(liveError === undefined ? {} : { liveError }),
					images,
				};
				const subject = rawParams.source === "pdb" ? `PDB ${identifier}` : identifier;
				const presentation = `${resolvedMode}, ${effectiveColor}${interfaceChain ? `, interface chain ${interfaceChain}${showInteractions ? " with interactions" : ""}` : ""}${showLigands ? ", ligands shown" : ", ligands hidden"}${residueColors.length > 0 ? `, ${residueColors.length} exact residue color${residueColors.length === 1 ? "" : "s"}` : ""}`;
				const snapshotSummary = `Rendered ${subject} with ProteinView at ${FULLHD_WIDTH}x${FULLHD_HEIGHT} (${presentation}).`;
				const text = live
					? `Opened ${subject} in the persistent ProteinView panel above the editor at FullHD ${FULLHD_WIDTH}x${FULLHD_HEIGHT} (${presentation}).`
					: liveError === undefined
						? snapshotSummary
						: `${snapshotSummary} Live panel unavailable (${liveError}); showing the FullHD snapshot.`;

				return {
					content: [
						{
							type: "text" as const,
							text,
						},
					],
					details,
				};
			} finally {
				await rm(privateDir, { recursive: true, force: true });
			}
		},
	};
}

const factory: CustomToolFactory = api => createProteinViewTool(api);

export default factory;

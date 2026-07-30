export interface ProteinViewLiveResidueColor {
	chain: string;
	residue_number: number;
	insertion_code?: string;
	color: string;
}

export interface ProteinViewLiveOpenRequest {
	inputPath: string;
	identifier: string;
	cwd: string;
	binary?: string;
	signal?: AbortSignal;
	mode?: "cartoon" | "backbone" | "wireframe";
	color: "structure" | "element" | "chain" | "plddt" | "bfactor" | "rainbow";
	/** Distinguishes an explicit `structure` request from the generic default. */
	explicitColor?: boolean;
	interfaceChain?: string;
	showInteractions: boolean;
	showLigands: boolean;
	residueColors: ProteinViewLiveResidueColor[];
}

export type ProteinViewLiveOpenHandler = (request: ProteinViewLiveOpenRequest) => Promise<void>;

let activeHandler: ProteinViewLiveOpenHandler | undefined;

/**
 * Connect the snapshot tool to the interactive extension without coupling
 * either module to the other's lifecycle. Only the active OMP session owns a
 * handler, and replacing it cleanly detaches the previous session.
 */
export function registerProteinViewLiveOpenHandler(handler: ProteinViewLiveOpenHandler): () => void {
	activeHandler = handler;
	return () => {
		if (activeHandler === handler) activeHandler = undefined;
	};
}

/**
 * Ask the active interactive OMP session to mount a live viewer.
 *
 * Returns false in print/RPC/headless modes, where no extension handler is
 * attached and the caller should keep the exact FullHD snapshot as fallback.
 */
export async function requestProteinViewLiveOpen(request: ProteinViewLiveOpenRequest): Promise<boolean> {
	const handler = activeHandler;
	if (handler === undefined) return false;
	await handler(request);
	return true;
}

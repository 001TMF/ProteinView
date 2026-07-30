#!/usr/bin/env python3
"""Fetch or copy one structure and render a persistent ProteinView FullHD PNG."""

from __future__ import annotations

import argparse
import hashlib
import http.client
import json
import os
from pathlib import Path
import re
import shutil
import stat
import struct
import subprocess
import sys
import tempfile
import time
import zlib

WIDTH = 1920
HEIGHT = 1080
MAX_STRUCTURE_BYTES = 64 * 1024 * 1024
MAX_PNG_BYTES = 16 * 1024 * 1024
MAX_CACHE_BYTES = 256 * 1024 * 1024
MAX_CACHE_FILES = 128
MAX_RESIDUE_COLORS = 256
FETCH_TIMEOUT_SECONDS = 30
RENDER_TIMEOUT_SECONDS = 60
PDB_ID = re.compile(r"^[1-9][A-Za-z0-9]{3}$")
SIGNED_I32 = re.compile(r"^[+-]?[0-9]+")
RGB_HEX = re.compile(r"^[0-9A-Fa-f]{6}$")
ALLOWED_EXTENSIONS = {".pdb", ".cif", ".mmcif"}
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
CACHE_FILENAME = re.compile(
    r"^(?:"
    r"[1-9][A-Z0-9]{3}\.cif|"
    r"local-[0-9a-f]{64}\.(?:pdb|cif|mmcif)|"
    r"(?:[1-9][A-Z0-9]{3}|local-[0-9a-f]{64})-"
    r"(?:"
    r"(?:auto|cartoon|backbone|wireframe)-"
    r"(?:structure|element|chain|plddt|bfactor|rainbow)|"
    r"view-[0-9a-f]{16}"
    r")-1920x1080\.png"
    r")$"
)
LIVE_SIDECAR_SUFFIX = ".proteinview.json"
LIVE_SIDECAR_SCHEMA = "proteinview-live-v1"


class RenderError(RuntimeError):
    """A concise user-actionable render failure."""


def normalize_pdb_id(value: str) -> str:
    if not PDB_ID.fullmatch(value):
        raise RenderError(
            f"invalid PDB ID {value!r}; expected four alphanumeric characters "
            "beginning with 1-9"
        )
    return value.upper()


def normalize_chain_id(value: str) -> str:
    return normalize_ascii_identifier(value, "chain ID")


def normalize_ascii_identifier(value: str, label: str) -> str:
    candidate = value.strip()
    try:
        encoded = candidate.encode("ascii")
    except UnicodeEncodeError as error:
        raise argparse.ArgumentTypeError(
            f"{label} must contain 1-32 printable non-whitespace ASCII bytes"
        ) from error
    if (
        not candidate
        or len(encoded) > 32
        or any(byte < 0x21 or byte > 0x7E for byte in encoded)
    ):
        raise argparse.ArgumentTypeError(
            f"{label} must contain 1-32 printable non-whitespace ASCII bytes"
        )
    return candidate


def normalize_residue_color(value: str) -> dict[str, object]:
    try:
        selector, color = value.rsplit("=", 1)
        chain, residue = selector.rsplit(":", 1)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "residue color must use CHAIN:RES[ICODE]=RRGGBB, "
            "for example A:42=FF0000"
        ) from error

    chain = normalize_ascii_identifier(chain, "chain ID")
    match = SIGNED_I32.match(residue)
    if match is None:
        raise argparse.ArgumentTypeError(
            "residue number must begin with a signed 32-bit integer"
        )
    try:
        residue_number = int(match.group(0))
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "residue number must be a signed 32-bit integer"
        ) from error
    if residue_number < -(2**31) or residue_number > 2**31 - 1:
        raise argparse.ArgumentTypeError(
            "residue number must be a signed 32-bit integer"
        )

    suffix = residue[match.end() :]
    if not suffix:
        insertion_code = None
    elif suffix.startswith("[") or suffix.endswith("]"):
        if not (suffix.startswith("[") and suffix.endswith("]") and len(suffix) > 2):
            raise argparse.ArgumentTypeError(
                "bracketed insertion code must use RES[ICODE], for example 42[A]"
            )
        insertion_code = normalize_ascii_identifier(
            suffix[1:-1], "insertion code"
        ).upper()
    else:
        insertion_code = normalize_ascii_identifier(suffix, "insertion code").upper()

    if RGB_HEX.fullmatch(color) is None:
        raise argparse.ArgumentTypeError(
            "residue color RGB must be exactly six hexadecimal digits, "
            "for example FF0000"
        )
    return {
        "chain": chain,
        "residue_number": residue_number,
        "insertion_code": insertion_code,
        "color": color.upper(),
    }


def residue_color_argument(residue: dict[str, object]) -> str:
    insertion = residue["insertion_code"]
    insertion_suffix = f"[{insertion}]" if insertion is not None else ""
    return (
        f"{residue['chain']}:{residue['residue_number']}"
        f"{insertion_suffix}={residue['color']}"
    )


def cache_directory() -> Path:
    configured = os.environ.get("PROTEINVIEW_CACHE_DIR")
    if configured:
        cache = (
            Path(configured).expanduser()
            / "proteinview"
            / "fullhd-v1"
        )
    else:
        user_id = getattr(os, "getuid", lambda: 0)()
        cache = Path(tempfile.gettempdir()) / f"proteinview-codex-{user_id}" / "fullhd-v1"
    if cache.is_symlink():
        raise RenderError(f"ProteinView cache path cannot be a symlink: {cache}")
    cache.mkdir(mode=0o700, parents=True, exist_ok=True)
    if hasattr(os, "getuid") and cache.stat().st_uid != os.getuid():
        raise RenderError(f"ProteinView cache is not owned by the current user: {cache}")
    try:
        cache.chmod(0o700)
    except OSError:
        pass
    return cache


def looks_like_mmcif(data: bytes) -> bool:
    sample = data[:8192].lstrip(b"\xef\xbb\xbf \t\r\n")
    while sample.startswith(b"#"):
        _, separator, sample = sample.partition(b"\n")
        if not separator:
            return False
        sample = sample.lstrip(b" \t\r\n")
    return sample.lower().startswith(b"data_")


def valid_cached_structure(path: Path) -> bool:
    if path.is_symlink():
        return False
    descriptor = None
    try:
        descriptor = os.open(
            path,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        )
        metadata = os.fstat(descriptor)
        if (
            not stat.S_ISREG(metadata.st_mode)
            or metadata.st_size <= 0
            or metadata.st_size > MAX_STRUCTURE_BYTES
        ):
            return False
        return looks_like_mmcif(os.read(descriptor, 8192))
    except OSError:
        return False
    finally:
        if descriptor is not None:
            os.close(descriptor)


def resolve_binary() -> str:
    configured = os.environ.get("PROTEINVIEW_BIN")
    if configured:
        binary = Path(configured).expanduser().resolve(strict=True)
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise RenderError(f"PROTEINVIEW_BIN is not executable: {binary}")
        return str(binary)
    discovered = shutil.which("proteinview")
    if discovered:
        binary = Path(discovered).resolve(strict=True)
        if binary.is_file() and os.access(binary, os.X_OK):
            return str(binary)
    raise RenderError(
        "ProteinView was not found. Install a build with --snapshot support or "
        "set PROTEINVIEW_BIN to that executable."
    )


def write_atomic(destination: Path, chunks) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            total = 0
            for chunk in chunks:
                total += len(chunk)
                if total > MAX_STRUCTURE_BYTES:
                    raise RenderError(
                        f"structure exceeds the {MAX_STRUCTURE_BYTES}-byte safety limit"
                    )
                output.write(chunk)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def fetch_structure(pdb_id: str, cache: Path) -> Path:
    destination = cache / f"{pdb_id}.cif"
    if valid_cached_structure(destination):
        return destination

    deadline = time.monotonic() + FETCH_TIMEOUT_SECONDS
    connection = http.client.HTTPSConnection(
        "files.rcsb.org",
        timeout=FETCH_TIMEOUT_SECONDS,
    )
    try:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise RenderError(
                f"RCSB fetch for PDB {pdb_id} exceeded {FETCH_TIMEOUT_SECONDS} seconds"
            )
        connection.timeout = remaining
        connection.request(
            "GET",
            f"/download/{pdb_id}.cif",
            headers={
                "Accept": "chemical/x-mmcif,text/plain;q=0.9",
                "User-Agent": "ProteinView-Codex/0.1",
            },
        )
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise RenderError(
                f"RCSB fetch for PDB {pdb_id} exceeded {FETCH_TIMEOUT_SECONDS} seconds"
            )
        if connection.sock is not None:
            connection.sock.settimeout(remaining)
        response = connection.getresponse()
        if response.status != 200:
            raise RenderError(
                f"failed to fetch PDB {pdb_id} from RCSB: HTTP {response.status}"
            )
        try:
            content_type = response.headers.get("Content-Type", "").lower()
            if "text/html" in content_type:
                raise RenderError(f"RCSB returned HTML instead of PDB {pdb_id}")
            declared = response.headers.get("Content-Length")
            if declared is not None and int(declared) > MAX_STRUCTURE_BYTES:
                raise RenderError(
                    f"PDB {pdb_id} exceeds the {MAX_STRUCTURE_BYTES}-byte safety limit"
                )

            downloaded = bytearray()
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise RenderError(
                        f"RCSB fetch for PDB {pdb_id} exceeded "
                        f"{FETCH_TIMEOUT_SECONDS} seconds"
                    )
                if connection.sock is not None:
                    connection.sock.settimeout(remaining)
                chunk = response.read(64 * 1024)
                if not chunk:
                    break
                downloaded.extend(chunk)
                if len(downloaded) > MAX_STRUCTURE_BYTES:
                    raise RenderError(
                        f"PDB {pdb_id} exceeds the {MAX_STRUCTURE_BYTES}-byte safety limit"
                    )
            if not looks_like_mmcif(downloaded):
                raise RenderError(f"RCSB returned invalid mmCIF data for PDB {pdb_id}")
            write_atomic(destination, (downloaded,))
        finally:
            response.close()
    except (http.client.HTTPException, OSError, ValueError) as error:
        raise RenderError(f"failed to fetch PDB {pdb_id} from RCSB: {error}") from error
    finally:
        connection.close()
    return destination


def copy_workspace_structure(source: str, cache: Path) -> tuple[Path, str]:
    workspace = Path.cwd().resolve()
    path = Path(source).expanduser().resolve(strict=True)
    try:
        path.relative_to(workspace)
    except ValueError as error:
        raise RenderError(f"local structure must be inside the workspace: {path}") from error
    if path.suffix.lower() not in ALLOWED_EXTENSIONS:
        raise RenderError("local structure must be a regular .pdb, .cif, or .mmcif file")

    digest = hashlib.sha256()
    input_descriptor = None
    output_descriptor = None
    temporary = None
    total = 0
    try:
        input_descriptor = os.open(
            path,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        )
        input_stat = os.fstat(input_descriptor)
        if (
            not stat.S_ISREG(input_stat.st_mode)
            or input_stat.st_size == 0
            or input_stat.st_size > MAX_STRUCTURE_BYTES
        ):
            raise RenderError(
                f"local structure must be a regular 1-{MAX_STRUCTURE_BYTES} byte file: "
                f"{path}"
            )

        output_descriptor, temporary_name = tempfile.mkstemp(
            prefix=".local-structure.", suffix=path.suffix.lower(), dir=cache
        )
        temporary = Path(temporary_name)
        with os.fdopen(input_descriptor, "rb") as input_file, os.fdopen(
            output_descriptor, "wb"
        ) as output:
            input_descriptor = None
            output_descriptor = None
            while chunk := input_file.read(64 * 1024):
                total += len(chunk)
                if total > MAX_STRUCTURE_BYTES:
                    raise RenderError(
                        f"local structure exceeds {MAX_STRUCTURE_BYTES} bytes"
                    )
                digest.update(chunk)
                output.write(chunk)
            output.flush()
            os.fsync(output.fileno())
        key = digest.hexdigest()
        destination = cache / f"local-{key}{path.suffix.lower()}"
        os.replace(temporary, destination)
    finally:
        if input_descriptor is not None:
            os.close(input_descriptor)
        if output_descriptor is not None:
            os.close(output_descriptor)
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return destination, path.name


def validate_png(path: Path) -> None:
    descriptor = os.open(
        path,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise RenderError("ProteinView output is not a regular PNG file")
        if metadata.st_size < 57 or metadata.st_size > MAX_PNG_BYTES:
            raise RenderError(
                f"ProteinView produced an invalid PNG size: {metadata.st_size} bytes"
            )
        with os.fdopen(descriptor, "rb") as image_file:
            descriptor = None
            encoded = image_file.read(MAX_PNG_BYTES + 1)
    finally:
        if descriptor is not None:
            os.close(descriptor)

    if len(encoded) > MAX_PNG_BYTES or not encoded.startswith(PNG_SIGNATURE):
        raise RenderError("ProteinView output is not a valid PNG")

    offset = len(PNG_SIGNATURE)
    width = height = None
    saw_idat = False
    saw_iend = False
    chunk_index = 0
    while offset + 12 <= len(encoded):
        length = struct.unpack(">I", encoded[offset : offset + 4])[0]
        chunk_end = offset + 12 + length
        if chunk_end > len(encoded):
            raise RenderError("ProteinView output contains a truncated PNG chunk")
        chunk_type = encoded[offset + 4 : offset + 8]
        chunk_data = encoded[offset + 8 : offset + 8 + length]
        expected_crc = struct.unpack(">I", encoded[offset + 8 + length : chunk_end])[0]
        actual_crc = zlib.crc32(chunk_type)
        actual_crc = zlib.crc32(chunk_data, actual_crc) & 0xFFFFFFFF
        if actual_crc != expected_crc:
            raise RenderError("ProteinView output contains an invalid PNG checksum")
        if chunk_index == 0:
            if chunk_type != b"IHDR" or length != 13:
                raise RenderError("ProteinView output is missing a valid PNG IHDR")
            width, height = struct.unpack(">II", chunk_data[:8])
        elif chunk_type == b"IDAT":
            saw_idat = True
        elif chunk_type == b"IEND":
            if length != 0 or chunk_end != len(encoded):
                raise RenderError("ProteinView output contains an invalid PNG ending")
            saw_iend = True
            offset = chunk_end
            break
        offset = chunk_end
        chunk_index += 1

    if not saw_idat or not saw_iend or offset != len(encoded):
        raise RenderError("ProteinView output is missing complete PNG image data")
    if (width, height) != (WIDTH, HEIGHT):
        raise RenderError(
            f"ProteinView output is {width}x{height}; expected {WIDTH}x{HEIGHT}"
        )


def render(
    binary: str,
    structure: Path,
    destination: Path,
    mode: str | None,
    color: str | None,
    interface_chain: str | None = None,
    show_interactions: bool = False,
    show_ligands: bool = True,
    residue_colors: list[dict[str, object]] | tuple[dict[str, object], ...] = (),
) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.stem}.", suffix=".png", dir=destination.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    temporary.unlink()
    command = [
        binary,
        str(structure),
        "--snapshot",
        str(temporary),
        "--snapshot-width",
        str(WIDTH),
        "--snapshot-height",
        str(HEIGHT),
    ]
    if mode is not None:
        command.extend(("--mode", mode))
    if color is not None:
        command.extend(("--color", color))
    if interface_chain is not None:
        command.extend(("--snapshot-interface-chain", interface_chain))
    if show_interactions:
        command.append("--snapshot-interactions")
    if not show_ligands:
        command.append("--snapshot-hide-ligands")
    for residue in residue_colors:
        command.extend(("--residue-color", residue_color_argument(residue)))
    try:
        completed = subprocess.run(
            command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=RENDER_TIMEOUT_SECONDS,
            check=False,
        )
        if completed.returncode != 0:
            detail = completed.stderr.strip()[-2000:] or completed.stdout.strip()[-2000:]
            if "snapshot" in detail.lower() and "unexpected" in detail.lower():
                detail += (
                    "\nInstall a newer ProteinView build with --snapshot support or "
                    "set PROTEINVIEW_BIN."
                )
            raise RenderError(
                f"ProteinView exited with status {completed.returncode}: {detail}"
            )
        validate_png(temporary)
        os.replace(temporary, destination)
    except subprocess.TimeoutExpired as error:
        raise RenderError(
            f"ProteinView exceeded the {RENDER_TIMEOUT_SECONDS}-second render timeout"
        ) from error
    finally:
        temporary.unlink(missing_ok=True)


def live_sidecar_path(image: Path) -> Path:
    return image.with_name(image.name + LIVE_SIDECAR_SUFFIX)


def write_live_sidecar(
    image: Path,
    structure: Path,
    source: str,
    mode: str,
    color: str | None,
    interface_chain: str | None,
    show_interactions: bool,
    show_ligands: bool,
    residue_colors: list[dict[str, object]] | tuple[dict[str, object], ...] = (),
) -> Path:
    sidecar = live_sidecar_path(image)
    payload = {
        "schema": LIVE_SIDECAR_SCHEMA,
        "structure_path": str(structure.resolve(strict=True)),
        "source": source,
        "mode": mode,
        "color": None if interface_chain else color,
        "interface_chain": interface_chain,
        "show_interactions": show_interactions,
        "show_ligands": show_ligands,
        "residue_colors": list(residue_colors),
    }
    encoded = (json.dumps(payload, sort_keys=True) + "\n").encode("utf-8")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{sidecar.name}.", suffix=".tmp", dir=sidecar.parent
    )
    temporary = Path(temporary_name)
    try:
        os.fchmod(descriptor, 0o600)
        with os.fdopen(descriptor, "wb") as output:
            descriptor = None
            output.write(encoded)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, sidecar)
    finally:
        if descriptor is not None:
            os.close(descriptor)
        temporary.unlink(missing_ok=True)
    return sidecar


def is_managed_cache_file(path: Path) -> bool:
    if CACHE_FILENAME.fullmatch(path.name):
        return True
    if not path.name.endswith(LIVE_SIDECAR_SUFFIX):
        return False
    image_name = path.name[: -len(LIVE_SIDECAR_SUFFIX)]
    return CACHE_FILENAME.fullmatch(image_name) is not None


def prune_cache(cache: Path, keep: set[Path]) -> None:
    files = [
        path
        for path in cache.iterdir()
        if (
            path.is_file()
            and not path.is_symlink()
            and path not in keep
            and is_managed_cache_file(path)
        )
    ]
    files.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    retained_bytes = sum(path.stat().st_size for path in keep if path.is_file())
    for index, path in enumerate(files):
        size = path.stat().st_size
        if index + len(keep) < MAX_CACHE_FILES and retained_bytes + size <= MAX_CACHE_BYTES:
            retained_bytes += size
            continue
        path.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render one ProteinView FullHD 1920x1080 PNG"
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--pdb", help="four-character RCSB PDB identifier")
    source.add_argument("--file", help="workspace-local PDB or mmCIF path")
    parser.add_argument(
        "--mode",
        choices=("cartoon", "backbone", "wireframe"),
        help="visualization mode; omit to let ProteinView auto-select for large structures",
    )
    parser.add_argument(
        "--color",
        choices=("structure", "element", "chain", "plddt", "bfactor", "rainbow"),
        help="protein color scheme; defaults to structure",
    )
    parser.add_argument(
        "--interface-chain",
        type=normalize_chain_id,
        help="focus chain for interface-residue coloring",
    )
    parser.add_argument(
        "--interactions",
        action="store_true",
        help="overlay classified interface interactions",
    )
    parser.add_argument(
        "--hide-ligands",
        action="store_true",
        help="hide ligands and ions",
    )
    parser.add_argument(
        "--residue-color",
        action="append",
        default=[],
        type=normalize_residue_color,
        metavar="CHAIN:RES[ICODE]=RRGGBB",
        help=(
            "exact polymer-residue color override; repeat for multiple residues "
            f"(maximum {MAX_RESIDUE_COLORS})"
        ),
    )
    args = parser.parse_args()
    if args.interactions and not args.interface_chain:
        parser.error("--interactions requires --interface-chain")
    if args.interface_chain and args.color:
        parser.error("--interface-chain uses the interface palette and cannot use --color")
    if len(args.residue_color) > MAX_RESIDUE_COLORS:
        parser.error(f"at most {MAX_RESIDUE_COLORS} --residue-color values are allowed")
    seen_residues = set()
    for residue in args.residue_color:
        key = (
            residue["chain"],
            residue["residue_number"],
            residue["insertion_code"],
        )
        if key in seen_residues:
            parser.error(
                "duplicate exact residue color target "
                f"{residue_color_argument(residue).rsplit('=', 1)[0]}"
            )
        seen_residues.add(key)
    return args


def main() -> int:
    try:
        args = parse_args()
        cache = cache_directory()
        if args.pdb:
            pdb_id = normalize_pdb_id(args.pdb)
            structure = fetch_structure(pdb_id, cache)
            source_label = f"PDB {pdb_id}"
            output_key = pdb_id
        else:
            structure, source_label = copy_workspace_structure(args.file, cache)
            output_key = structure.stem

        mode_label = args.mode or "auto"
        color_label = "interface" if args.interface_chain else (args.color or "structure")
        presentation = {
            "color": color_label,
            "interface_chain": args.interface_chain,
            "mode": mode_label,
            "residue_colors": args.residue_color,
            "show_interactions": args.interactions,
            "show_ligands": not args.hide_ligands,
        }
        view_key = hashlib.sha256(
            json.dumps(presentation, sort_keys=True).encode("utf-8")
        ).hexdigest()[:16]
        output = cache / f"{output_key}-view-{view_key}-{WIDTH}x{HEIGHT}.png"
        binary = resolve_binary()
        try:
            validate_png(output)
        except (FileNotFoundError, RenderError):
            render(
                binary,
                structure,
                output,
                args.mode,
                args.color,
                args.interface_chain,
                args.interactions,
                not args.hide_ligands,
                args.residue_color,
            )
        sidecar = write_live_sidecar(
            output,
            structure,
            source_label,
            mode_label,
            None if args.interface_chain else (args.color or "structure"),
            args.interface_chain,
            args.interactions,
            not args.hide_ligands,
            args.residue_color,
        )
        prune_cache(cache, {structure, output, sidecar})
        print(
            json.dumps(
                {
                    "path": str(output.resolve()),
                    "source": source_label,
                    "mode": mode_label,
                    "color": color_label,
                    "width": WIDTH,
                    "height": HEIGHT,
                    "interface_chain": args.interface_chain,
                    "interactions": args.interactions,
                    "ligands": not args.hide_ligands,
                    "residue_colors": args.residue_color,
                    "live_sidecar": str(sidecar.resolve()),
                    "renderer": "ProteinView FullHD",
                },
                sort_keys=True,
            )
        )
        return 0
    except (OSError, RenderError) as error:
        print(f"proteinview-render: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

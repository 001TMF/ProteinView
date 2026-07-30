import importlib.util
import io
import json
import os
from pathlib import Path
import struct
import tempfile
import unittest
from unittest import mock
import zlib


SCRIPT = (
    Path(__file__).parents[1]
    / "proteinview"
    / "skills"
    / "render-protein"
    / "scripts"
    / "render_protein.py"
)
SPEC = importlib.util.spec_from_file_location("render_protein", SCRIPT)
RENDER_PROTEIN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RENDER_PROTEIN)


def png_chunk(chunk_type, data):
    checksum = zlib.crc32(chunk_type)
    checksum = zlib.crc32(data, checksum) & 0xFFFFFFFF
    return struct.pack(">I", len(data)) + chunk_type + data + struct.pack(">I", checksum)


def valid_png(width=1920, height=1080):
    ihdr = struct.pack(">IIBBBBB", width, height, 1, 0, 0, 0, 0)
    row = b"\0" + (b"\0" * ((width + 7) // 8))
    pixels = row * height
    return (
        RENDER_PROTEIN.PNG_SIGNATURE
        + png_chunk(b"IHDR", ihdr)
        + png_chunk(b"IDAT", zlib.compress(pixels))
        + png_chunk(b"IEND", b"")
    )


class RenderProteinTests(unittest.TestCase):
    def test_normalizes_strict_pdb_identifier(self):
        self.assertEqual(RENDER_PROTEIN.normalize_pdb_id("1ubq"), "1UBQ")
        for invalid in ("ABCD", "../x", "1UBQ.cif", "1 UQ", "0ABC"):
            with self.assertRaises(RENDER_PROTEIN.RenderError):
                RENDER_PROTEIN.normalize_pdb_id(invalid)

    def test_normalizes_exact_residue_colors_to_canonical_form(self):
        plain = RENDER_PROTEIN.normalize_residue_color("AA:-002=ff00a7")
        inserted = RENDER_PROTEIN.normalize_residue_color("AA:42[a]=00ff7f")

        self.assertEqual(
            plain,
            {
                "chain": "AA",
                "residue_number": -2,
                "insertion_code": None,
                "color": "FF00A7",
            },
        )
        self.assertEqual(
            inserted,
            {
                "chain": "AA",
                "residue_number": 42,
                "insertion_code": "A",
                "color": "00FF7F",
            },
        )
        self.assertEqual(
            RENDER_PROTEIN.residue_color_argument(inserted),
            "AA:42[A]=00FF7F",
        )

        for invalid in (
            "A42=FF0000",
            "A:x=FF0000",
            "A:2147483648=FF0000",
            "A:42[]=FF0000",
            "A:42=#FF0000",
            "Å:42=FF0000",
        ):
            with self.assertRaises(RENDER_PROTEIN.argparse.ArgumentTypeError):
                RENDER_PROTEIN.normalize_residue_color(invalid)

    def test_rejects_duplicate_exact_residue_colors(self):
        stderr = io.StringIO()
        with (
            mock.patch.object(
                RENDER_PROTEIN.sys,
                "argv",
                [
                    "render_protein.py",
                    "--pdb",
                    "1UBQ",
                    "--residue-color",
                    "A:42[a]=ff0000",
                    "--residue-color",
                    "A:42[A]=00FF00",
                ],
            ),
            mock.patch.object(RENDER_PROTEIN.sys, "stderr", stderr),
            self.assertRaises(SystemExit),
        ):
            RENDER_PROTEIN.parse_args()
        self.assertIn("duplicate exact residue color target A:42[A]", stderr.getvalue())

    def test_validates_exact_fullhd_dimensions(self):
        with tempfile.TemporaryDirectory() as directory:
            image = Path(directory) / "render.png"
            image.write_bytes(valid_png())
            RENDER_PROTEIN.validate_png(image)
            image.write_bytes(valid_png(width=1280, height=720))
            with self.assertRaises(RENDER_PROTEIN.RenderError):
                RENDER_PROTEIN.validate_png(image)
            image.write_bytes(valid_png()[:33])
            with self.assertRaises(RENDER_PROTEIN.RenderError):
                RENDER_PROTEIN.validate_png(image)

    def test_configured_cache_uses_dedicated_child_and_preserves_other_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            unrelated = root / "notes.txt"
            unrelated.write_text("keep", encoding="utf-8")
            with mock.patch.dict(
                os.environ, {"PROTEINVIEW_CACHE_DIR": str(root)}, clear=False
            ):
                cache = RENDER_PROTEIN.cache_directory()
            self.assertEqual(cache, root / "proteinview" / "fullhd-v1")
            RENDER_PROTEIN.prune_cache(cache, set())
            self.assertEqual(unrelated.read_text(encoding="utf-8"), "keep")

    def test_rejects_html_as_cached_mmcif(self):
        with tempfile.TemporaryDirectory() as directory:
            cached = Path(directory) / "1UBQ.cif"
            cached.write_text("<html>temporary error</html>", encoding="utf-8")
            self.assertFalse(RENDER_PROTEIN.valid_cached_structure(cached))

    def test_fetch_enforces_an_overall_deadline(self):
        class Response:
            status = 200
            headers = {}
            read_called = False
            closed = False

            def read(self, _size):
                self.read_called = True
                return b"data_1UBQ\n"

            def close(self):
                self.closed = True

        class Connection:
            def __init__(self, response):
                self.response = response
                self.sock = mock.Mock()
                self.closed = False
                self.getresponse_called = False

            def request(self, *_args, **_kwargs):
                return None

            def getresponse(self):
                self.getresponse_called = True
                return self.response

            def close(self):
                self.closed = True

        with tempfile.TemporaryDirectory() as directory:
            response = Response()
            connection = Connection(response)
            with (
                mock.patch.object(
                    RENDER_PROTEIN.http.client,
                    "HTTPSConnection",
                    return_value=connection,
                ),
                mock.patch.object(
                    RENDER_PROTEIN.time,
                    "monotonic",
                    side_effect=(100.0, 100.0, 131.0),
                ),
            ):
                with self.assertRaisesRegex(
                    RENDER_PROTEIN.RenderError, "exceeded 30 seconds"
                ):
                    RENDER_PROTEIN.fetch_structure("1UBQ", Path(directory))
            self.assertFalse(response.read_called)
            self.assertFalse(response.closed)
            self.assertFalse(connection.getresponse_called)
            self.assertTrue(connection.closed)

    def test_copies_only_regular_workspace_structures(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / "cache"
            cache.mkdir()
            structure = root / "input.pdb"
            structure.write_bytes(b"ATOM\n")

            with mock.patch.object(RENDER_PROTEIN.Path, "cwd", return_value=root):
                copied, label = RENDER_PROTEIN.copy_workspace_structure(
                    str(structure), cache
                )
                self.assertEqual(label, "input.pdb")
                self.assertEqual(copied.read_bytes(), b"ATOM\n")

                directory_with_extension = root / "not-a-file.pdb"
                directory_with_extension.mkdir()
                with self.assertRaises(RENDER_PROTEIN.RenderError):
                    RENDER_PROTEIN.copy_workspace_structure(
                        str(directory_with_extension), cache
                    )

    def test_renderer_uses_fixed_argv_and_no_stdin(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            structure = root / "1UBQ.cif"
            structure.write_text("data_1UBQ\n", encoding="utf-8")
            destination = root / "1UBQ.png"

            def run(command, **options):
                output = Path(command[command.index("--snapshot") + 1])
                output.write_bytes(valid_png())
                self.assertEqual(command[0], "/fixed/proteinview")
                self.assertEqual(command[1], str(structure))
                self.assertNotIn("--mode", command)
                self.assertNotIn("--color", command)
                self.assertIs(options["stdin"], RENDER_PROTEIN.subprocess.DEVNULL)
                self.assertFalse(options["check"])
                return mock.Mock(returncode=0, stderr="", stdout="")

            with mock.patch.object(RENDER_PROTEIN.subprocess, "run", side_effect=run):
                RENDER_PROTEIN.render(
                    "/fixed/proteinview",
                    structure,
                    destination,
                    None,
                    None,
                )
            self.assertEqual(destination.read_bytes(), valid_png())

    def test_renderer_passes_repeatable_exact_residue_colors(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            structure = root / "1UBQ.cif"
            structure.write_text("data_1UBQ\n", encoding="utf-8")
            destination = root / "1UBQ.png"
            residues = [
                RENDER_PROTEIN.normalize_residue_color("A:1=ff0000"),
                RENDER_PROTEIN.normalize_residue_color("A:42[b]=00ff7f"),
            ]

            def run(command, **_options):
                output = Path(command[command.index("--snapshot") + 1])
                output.write_bytes(valid_png())
                self.assertEqual(
                    [
                        command[index + 1]
                        for index, value in enumerate(command)
                        if value == "--residue-color"
                    ],
                    ["A:1=FF0000", "A:42[B]=00FF7F"],
                )
                return mock.Mock(returncode=0, stderr="", stdout="")

            with mock.patch.object(RENDER_PROTEIN.subprocess, "run", side_effect=run):
                RENDER_PROTEIN.render(
                    "/fixed/proteinview",
                    structure,
                    destination,
                    None,
                    None,
                    residue_colors=residues,
                )

    def test_writes_private_live_viewer_sidecar(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            image = root / "1UBQ-view-0123456789abcdef-1920x1080.png"
            image.write_bytes(valid_png())
            structure = root / "1UBQ.cif"
            structure.write_text("data_1UBQ\n", encoding="utf-8")

            sidecar = RENDER_PROTEIN.write_live_sidecar(
                image,
                structure,
                "PDB 1UBQ",
                "cartoon",
                "chain",
                None,
                False,
                True,
            )

            self.assertEqual(
                sidecar,
                Path(str(image) + RENDER_PROTEIN.LIVE_SIDECAR_SUFFIX),
            )
            self.assertEqual(
                json.loads(sidecar.read_text(encoding="utf-8")),
                {
                    "schema": "proteinview-live-v1",
                    "structure_path": str(structure.resolve()),
                    "source": "PDB 1UBQ",
                    "mode": "cartoon",
                    "color": "chain",
                    "interface_chain": None,
                    "residue_colors": [],
                    "show_interactions": False,
                    "show_ligands": True,
                },
            )
            self.assertEqual(sidecar.stat().st_mode & 0o077, 0)
            self.assertTrue(RENDER_PROTEIN.is_managed_cache_file(sidecar))

    def test_interface_sidecar_uses_null_regular_color(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            image = root / "1UBQ-view-0123456789abcdef-1920x1080.png"
            structure = root / "1UBQ.cif"
            structure.write_text("data_1UBQ\n", encoding="utf-8")

            sidecar = RENDER_PROTEIN.write_live_sidecar(
                image,
                structure,
                "PDB 1UBQ",
                "auto",
                "structure",
                "A",
                True,
                False,
            )

            payload = json.loads(sidecar.read_text(encoding="utf-8"))
            self.assertIsNone(payload["color"])
            self.assertEqual(payload["interface_chain"], "A")
            self.assertEqual(payload["residue_colors"], [])
            self.assertTrue(payload["show_interactions"])
            self.assertFalse(payload["show_ligands"])

    def test_main_emits_live_sidecar_for_view_image_result(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / "cache"
            cache.mkdir()
            structure = cache / "1UBQ.cif"
            structure.write_text("data_1UBQ\n", encoding="utf-8")
            binary = root / "proteinview"
            binary.write_text("#!/bin/sh\n", encoding="utf-8")
            binary.chmod(0o700)

            def render(_binary, _structure, destination, *_args):
                destination.write_bytes(valid_png())

            stdout = io.StringIO()
            with (
                mock.patch.object(
                    RENDER_PROTEIN.sys,
                    "argv",
                    ["render_protein.py", "--pdb", "1UBQ", "--color", "chain"],
                ),
                mock.patch.object(
                    RENDER_PROTEIN, "cache_directory", return_value=cache
                ),
                mock.patch.object(
                    RENDER_PROTEIN, "fetch_structure", return_value=structure
                ),
                mock.patch.object(
                    RENDER_PROTEIN, "resolve_binary", return_value=str(binary)
                ),
                mock.patch.object(RENDER_PROTEIN, "render", side_effect=render),
                mock.patch.object(RENDER_PROTEIN.sys, "stdout", stdout),
            ):
                self.assertEqual(RENDER_PROTEIN.main(), 0)

            result = json.loads(stdout.getvalue())
            image = Path(result["path"])
            sidecar = Path(result["live_sidecar"])
            self.assertEqual(sidecar, RENDER_PROTEIN.live_sidecar_path(image))
            self.assertTrue(image.is_file())
            self.assertEqual(
                json.loads(sidecar.read_text(encoding="utf-8")),
                {
                    "schema": "proteinview-live-v1",
                    "structure_path": str(structure.resolve()),
                    "source": "PDB 1UBQ",
                    "mode": "auto",
                    "color": "chain",
                    "interface_chain": None,
                    "residue_colors": [],
                    "show_interactions": False,
                    "show_ligands": True,
                },
            )

    def test_main_keys_and_reports_exact_residue_presentation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / "cache"
            cache.mkdir()
            structure = cache / "1UBQ.cif"
            structure.write_text("data_1UBQ\n", encoding="utf-8")
            binary = root / "proteinview"
            binary.write_text("#!/bin/sh\n", encoding="utf-8")
            binary.chmod(0o700)

            def render(_binary, _structure, destination, *_args):
                destination.write_bytes(valid_png())

            def run_with_color(color):
                stdout = io.StringIO()
                with (
                    mock.patch.object(
                        RENDER_PROTEIN.sys,
                        "argv",
                        [
                            "render_protein.py",
                            "--pdb",
                            "1UBQ",
                            "--residue-color",
                            color,
                        ],
                    ),
                    mock.patch.object(
                        RENDER_PROTEIN, "cache_directory", return_value=cache
                    ),
                    mock.patch.object(
                        RENDER_PROTEIN, "fetch_structure", return_value=structure
                    ),
                    mock.patch.object(
                        RENDER_PROTEIN, "resolve_binary", return_value=str(binary)
                    ),
                    mock.patch.object(RENDER_PROTEIN, "render", side_effect=render),
                    mock.patch.object(RENDER_PROTEIN.sys, "stdout", stdout),
                ):
                    self.assertEqual(RENDER_PROTEIN.main(), 0)
                return json.loads(stdout.getvalue())

            red = run_with_color("A:1=ff0000")
            green = run_with_color("A:1=00ff00")

            self.assertNotEqual(red["path"], green["path"])
            self.assertEqual(
                red["residue_colors"],
                [
                    {
                        "chain": "A",
                        "residue_number": 1,
                        "insertion_code": None,
                        "color": "FF0000",
                    }
                ],
            )
            red_sidecar = json.loads(
                Path(red["live_sidecar"]).read_text(encoding="utf-8")
            )
            self.assertEqual(red_sidecar["residue_colors"], red["residue_colors"])

    def test_failed_view_variant_preserves_shared_cached_structure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / "cache"
            cache.mkdir()
            structure = cache / "1UBQ.cif"
            structure.write_text("data_1UBQ\n", encoding="utf-8")
            stderr = io.StringIO()

            with (
                mock.patch.object(
                    RENDER_PROTEIN.sys,
                    "argv",
                    [
                        "render_protein.py",
                        "--pdb",
                        "1UBQ",
                        "--residue-color",
                        "A:9999=FF0000",
                    ],
                ),
                mock.patch.object(
                    RENDER_PROTEIN, "cache_directory", return_value=cache
                ),
                mock.patch.object(
                    RENDER_PROTEIN, "fetch_structure", return_value=structure
                ),
                mock.patch.object(
                    RENDER_PROTEIN, "resolve_binary", return_value="/fixed/proteinview"
                ),
                mock.patch.object(
                    RENDER_PROTEIN,
                    "render",
                    side_effect=RENDER_PROTEIN.RenderError("unknown residue"),
                ),
                mock.patch.object(RENDER_PROTEIN.sys, "stderr", stderr),
            ):
                self.assertEqual(RENDER_PROTEIN.main(), 1)

            self.assertTrue(structure.is_file())
            self.assertIn("unknown residue", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()

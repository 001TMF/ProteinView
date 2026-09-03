//! User-configurable color palette.
//!
//! Every fixed color ProteinView draws comes from a [`Palette`].  The built-in
//! defaults reproduce the previously hardcoded colors exactly, so a user with no
//! config file sees no change.  A TOML file may override any subset of them:
//! anything omitted keeps its default.
//!
//! The palette is resolved once at startup and read-only thereafter, so it lives
//! in a process-wide [`OnceLock`] rather than being threaded through every
//! renderer.  Call [`init`] once from `main`; everything else reads [`palette`].
//!
//! Procedural schemes (Rainbow's HSV sweep) and the depth fog are not covered
//! here yet.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use anyhow::{Context, Result};
use serde::Deserialize;

/// An RGB triple, deserialized from a hex string such as `"FF0080"` or `"#ff0080"`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rgb(pub [u8; 3]);

impl Rgb {
    pub const fn new(r: u8, g: u8, b: u8) -> Self {
        Self([r, g, b])
    }
}

impl<'de> Deserialize<'de> for Rgb {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        use serde::de::Error;
        let raw = String::deserialize(deserializer)?;
        let hex = raw.strip_prefix('#').unwrap_or(&raw);
        if hex.len() != 6 || !hex.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err(D::Error::custom(format!(
                "expected six hexadecimal digits such as \"FF0080\", got {raw:?}"
            )));
        }
        let byte = |i: usize| u8::from_str_radix(&hex[i..i + 2], 16).map_err(D::Error::custom);
        Ok(Rgb([byte(0)?, byte(2)?, byte(4)?]))
    }
}

// ---------------------------------------------------------------------------
// Sections
// ---------------------------------------------------------------------------

/// Secondary-structure colors for the Structure scheme.
#[derive(Debug, Clone, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct StructurePalette {
    pub helix: Rgb,
    pub sheet: Rgb,
    pub turn: Rgb,
    pub coil: Rgb,
}

impl Default for StructurePalette {
    fn default() -> Self {
        Self {
            helix: Rgb::new(255, 0, 128),
            sheet: Rgb::new(255, 200, 0),
            turn: Rgb::new(96, 128, 255),
            coil: Rgb::new(0, 204, 0),
        }
    }
}

/// Per-base colors for nucleic acid residues.
#[derive(Debug, Clone, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct NucleotidePalette {
    pub adenine: Rgb,
    pub uracil: Rgb,
    pub thymine: Rgb,
    pub guanine: Rgb,
    pub cytosine: Rgb,
    pub inosine: Rgb,
}

impl Default for NucleotidePalette {
    fn default() -> Self {
        Self {
            adenine: Rgb::new(220, 60, 60),
            uracil: Rgb::new(60, 60, 220),
            thymine: Rgb::new(60, 60, 220),
            guanine: Rgb::new(60, 180, 60),
            cytosine: Rgb::new(220, 200, 40),
            inosine: Rgb::new(150, 100, 180),
        }
    }
}

/// CPK-style element colors.  `symbols` is keyed by uppercase element symbol.
#[derive(Debug, Clone)]
pub struct ElementPalette {
    pub fallback: Rgb,
    pub symbols: HashMap<String, Rgb>,
}

impl ElementPalette {
    /// Color for an element symbol, case-insensitively.
    pub fn get(&self, symbol: &str) -> Rgb {
        self.symbols
            .get(&symbol.trim().to_ascii_uppercase())
            .copied()
            .unwrap_or(self.fallback)
    }
}

impl Default for ElementPalette {
    fn default() -> Self {
        const CPK: &[(&str, [u8; 3])] = &[
            ("C", [144, 144, 144]),
            ("N", [48, 80, 248]),
            ("O", [255, 13, 13]),
            ("S", [255, 255, 48]),
            ("H", [255, 255, 255]),
            ("P", [255, 128, 0]),
            ("FE", [224, 102, 51]),
            ("MG", [0, 180, 0]),
            ("ZN", [125, 128, 176]),
            ("CA", [61, 255, 0]),
            ("MN", [156, 122, 199]),
            ("CO", [240, 144, 160]),
            ("CU", [200, 128, 51]),
            ("NI", [80, 208, 80]),
            ("CL", [31, 240, 31]),
            ("BR", [166, 41, 41]),
        ];
        Self {
            fallback: Rgb::new(200, 200, 200),
            symbols: CPK
                .iter()
                .map(|(s, c)| ((*s).to_string(), Rgb(*c)))
                .collect(),
        }
    }
}

/// AlphaFold pLDDT confidence bands.
#[derive(Debug, Clone, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PlddtPalette {
    /// pLDDT >= 90
    pub very_high: Rgb,
    /// pLDDT >= 70
    pub high: Rgb,
    /// pLDDT >= 50
    pub low: Rgb,
    /// pLDDT < 50
    pub very_low: Rgb,
}

impl Default for PlddtPalette {
    fn default() -> Self {
        Self {
            very_high: Rgb::new(0, 83, 214),
            high: Rgb::new(101, 203, 243),
            low: Rgb::new(255, 219, 19),
            very_low: Rgb::new(255, 125, 69),
        }
    }
}

/// Endpoints of the B-factor gradient.
#[derive(Debug, Clone, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct BFactorPalette {
    /// Color at the cold end of the range.
    pub low: Rgb,
    /// Color at the hot end of the range.
    pub high: Rgb,
}

impl Default for BFactorPalette {
    fn default() -> Self {
        Self {
            low: Rgb::new(0, 0, 255),
            high: Rgb::new(255, 0, 0),
        }
    }
}

/// Interface-highlighting colors.
#[derive(Debug, Clone, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct InterfacePalette {
    /// Focus chain, in contact with a partner.
    pub focus_contact: Rgb,
    /// Focus chain, away from the interface.
    pub focus_other: Rgb,
    /// Partner chain, in contact with the focus chain.
    pub partner_contact: Rgb,
    /// Partner chain, away from the interface.
    pub partner_other: Rgb,
    /// Ligands, drawn bright so they stand out against the interface coloring.
    pub ligand: Rgb,
}

impl Default for InterfacePalette {
    fn default() -> Self {
        Self {
            focus_contact: Rgb::new(0, 255, 100),
            focus_other: Rgb::new(40, 100, 60),
            partner_contact: Rgb::new(255, 165, 0),
            partner_other: Rgb::new(100, 80, 60),
            ligand: Rgb::new(255, 255, 255),
        }
    }
}

/// Small-molecule colors.
#[derive(Debug, Clone, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct LigandPalette {
    pub ligand: Rgb,
    pub ion: Rgb,
    /// Ligands under the Rainbow scheme, which has no per-residue value for them.
    pub rainbow: Rgb,
}

impl Default for LigandPalette {
    fn default() -> Self {
        Self {
            ligand: Rgb::new(255, 0, 255),
            ion: Rgb::new(0, 255, 255),
            rainbow: Rgb::new(255, 0, 255),
        }
    }
}

// ---------------------------------------------------------------------------
// Palette
// ---------------------------------------------------------------------------

/// The resolved palette: built-in defaults with any configured overrides applied.
#[derive(Debug, Clone)]
pub struct Palette {
    pub structure: StructurePalette,
    pub nucleotide: NucleotidePalette,
    /// Cycled by chain in order; always at least one entry.
    pub chains: Vec<Rgb>,
    pub element: ElementPalette,
    pub plddt: PlddtPalette,
    pub bfactor: BFactorPalette,
    pub interface: InterfacePalette,
    pub ligand: LigandPalette,
}

impl Default for Palette {
    fn default() -> Self {
        Self {
            structure: StructurePalette::default(),
            nucleotide: NucleotidePalette::default(),
            chains: vec![
                Rgb::new(0, 180, 255),
                Rgb::new(255, 100, 0),
                Rgb::new(0, 220, 100),
                Rgb::new(255, 50, 150),
                Rgb::new(180, 100, 255),
                Rgb::new(255, 220, 0),
                Rgb::new(0, 200, 200),
                Rgb::new(255, 150, 150),
            ],
            element: ElementPalette::default(),
            plddt: PlddtPalette::default(),
            bfactor: BFactorPalette::default(),
            interface: InterfacePalette::default(),
            ligand: LigandPalette::default(),
        }
    }
}

impl Palette {
    /// Color for the chain whose id starts with `id`, cycling through `chains`.
    pub fn chain(&self, id: &str) -> Rgb {
        let idx = id.bytes().next().unwrap_or(b'A') as usize % self.chains.len();
        self.chains[idx]
    }
}

// ---------------------------------------------------------------------------
// Config file
// ---------------------------------------------------------------------------

/// Element section as written in the file.  `symbols` is *merged* onto the
/// built-in CPK table rather than replacing it, so overriding carbon does not
/// silently drop every other element.
#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct ElementFile {
    /// Color for elements not listed in `symbols`.
    fallback: Option<Rgb>,
    symbols: HashMap<String, Rgb>,
}

/// Chain section.  A section rather than a bare top-level `chains` key, because
/// in TOML a bare key written after any `[table]` header silently binds to that
/// table instead of the document root.
#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct ChainFile {
    /// Replaces the default cycle outright; order is significant.
    colors: Option<Vec<Rgb>>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct PaletteFile {
    structure: StructurePalette,
    nucleotide: NucleotidePalette,
    chain: ChainFile,
    element: ElementFile,
    plddt: PlddtPalette,
    bfactor: BFactorPalette,
    interface: InterfacePalette,
    ligand: LigandPalette,
}

impl PaletteFile {
    fn resolve(self) -> Result<Palette> {
        let mut palette = Palette {
            structure: self.structure,
            nucleotide: self.nucleotide,
            chains: Palette::default().chains,
            element: ElementPalette::default(),
            plddt: self.plddt,
            bfactor: self.bfactor,
            interface: self.interface,
            ligand: self.ligand,
        };

        if let Some(colors) = self.chain.colors {
            if colors.is_empty() {
                anyhow::bail!("`chain.colors` must list at least one color");
            }
            palette.chains = colors;
        }
        if let Some(fallback) = self.element.fallback {
            palette.element.fallback = fallback;
        }
        for (symbol, color) in self.element.symbols {
            palette
                .element
                .symbols
                .insert(symbol.trim().to_ascii_uppercase(), color);
        }

        Ok(palette)
    }
}

/// Parse a palette from TOML text, filling anything omitted from the defaults.
pub fn parse(text: &str) -> Result<Palette> {
    toml::from_str::<PaletteFile>(text)?.resolve()
}

/// Read and parse a palette file.
pub fn load(path: &Path) -> Result<Palette> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("cannot read palette file {}", path.display()))?;
    parse(&text).with_context(|| format!("invalid palette file {}", path.display()))
}

/// The default palette file location: `$XDG_CONFIG_HOME/proteinview/palette.toml`,
/// falling back to `~/.config/proteinview/palette.toml`.
pub fn default_config_path() -> Option<PathBuf> {
    let base = std::env::var_os("XDG_CONFIG_HOME")
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|home| PathBuf::from(home).join(".config")))?;
    Some(base.join("proteinview").join("palette.toml"))
}

static PALETTE: OnceLock<Palette> = OnceLock::new();

/// The active palette.  Defaults are used if [`init`] was never called.
pub fn palette() -> &'static Palette {
    PALETTE.get_or_init(Palette::default)
}

/// Resolve the palette once, from `explicit` if given, else from the default
/// config path if it exists, else from the built-in defaults.
///
/// An explicit path that cannot be read is an error; so is a malformed file in
/// either location, since silently ignoring a palette the user wrote is worse
/// than refusing to start.  Returns the path that was loaded, if any.
pub fn init(explicit: Option<&Path>) -> Result<Option<PathBuf>> {
    let chosen = match explicit {
        Some(path) => Some(path.to_path_buf()),
        None => default_config_path().filter(|p| p.is_file()),
    };

    let (resolved, loaded) = match &chosen {
        Some(path) => (load(path)?, Some(path.clone())),
        None => (Palette::default(), None),
    };

    let _ = PALETTE.set(resolved);
    Ok(loaded)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_config_is_the_default_palette() {
        let p = parse("").unwrap();
        let d = Palette::default();
        assert_eq!(p.structure.helix, d.structure.helix);
        assert_eq!(p.chains, d.chains);
        assert_eq!(p.element.get("ZN"), d.element.get("ZN"));
        assert_eq!(p.plddt.very_high, d.plddt.very_high);
    }

    #[test]
    fn overrides_apply_and_the_rest_stays_default() {
        let p = parse(
            r##"
            [structure]
            helix = "#112233"
            "##,
        )
        .unwrap();
        assert_eq!(p.structure.helix, Rgb::new(0x11, 0x22, 0x33));
        // Untouched fields keep their defaults.
        assert_eq!(p.structure.sheet, StructurePalette::default().sheet);
        assert_eq!(p.structure.coil, StructurePalette::default().coil);
    }

    #[test]
    fn element_symbols_merge_rather_than_replace() {
        // Overriding carbon must not drop the rest of the CPK table.
        let p = parse(
            r#"
            [element.symbols]
            C = "010203"
            "#,
        )
        .unwrap();
        assert_eq!(p.element.get("C"), Rgb::new(1, 2, 3));
        assert_eq!(p.element.get("ZN"), Rgb::new(125, 128, 176));
        assert_eq!(p.element.get("FE"), Rgb::new(224, 102, 51));
    }

    #[test]
    fn element_lookup_is_case_insensitive_both_ways() {
        let p = parse(
            r#"
            [element.symbols]
            se = "0A0B0C"
            "#,
        )
        .unwrap();
        assert_eq!(p.element.get("SE"), Rgb::new(10, 11, 12));
        assert_eq!(p.element.get("Se"), Rgb::new(10, 11, 12));
        assert_eq!(p.element.get(" fe "), Rgb::new(224, 102, 51));
    }

    #[test]
    fn chains_are_replaced_wholesale_and_cycle() {
        let p = parse("[chain]\ncolors = [\"FF0000\", \"00FF00\"]").unwrap();
        assert_eq!(p.chains.len(), 2);
        // b'A' = 65, so 65 % 2 = 1 picks the second entry.
        assert_eq!(p.chain("A"), Rgb::new(0, 255, 0));
        assert_eq!(p.chain("B"), Rgb::new(255, 0, 0));
        assert_eq!(p.chain("A"), p.chain("C"));
    }

    #[test]
    fn empty_chain_list_is_rejected() {
        let err = parse("[chain]\ncolors = []").unwrap_err().to_string();
        assert!(err.contains("at least one"), "unhelpful error: {err}");
    }

    #[test]
    fn chain_colors_survive_being_written_after_other_sections() {
        // A bare top-level `chains = [...]` would bind to whichever [table]
        // preceded it. Keeping it in its own section makes order irrelevant.
        let p = parse(
            "[structure]\nhelix = \"FF0000\"\n\n[chain]\ncolors = [\"00FF00\"]",
        )
        .unwrap();
        assert_eq!(p.chains, vec![Rgb::new(0, 255, 0)]);
        assert_eq!(p.structure.helix, Rgb::new(255, 0, 0));
    }

    #[test]
    fn hex_accepts_optional_hash_and_any_case() {
        let p = parse(
            r##"
            [structure]
            helix = "#aabbcc"
            sheet = "AABBCC"
            "##,
        )
        .unwrap();
        assert_eq!(p.structure.helix, Rgb::new(0xAA, 0xBB, 0xCC));
        assert_eq!(p.structure.sheet, p.structure.helix);
    }

    #[test]
    fn malformed_color_is_rejected_with_the_offending_value() {
        let err = parse(
            r#"
            [structure]
            helix = "nope"
            "#,
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("six hexadecimal"), "unhelpful error: {err}");
        assert!(err.contains("nope"), "error should quote the value: {err}");
    }

    #[test]
    fn unknown_keys_are_rejected_so_typos_are_not_silent() {
        let err = parse(
            r#"
            [structure]
            helics = "FF0000"
            "#,
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("helics"), "error should name the key: {err}");
    }
}

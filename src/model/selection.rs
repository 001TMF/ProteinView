use std::collections::{HashMap, HashSet};
use std::fmt;
use std::str::FromStr;

use anyhow::{Result, bail};

use crate::model::protein::{Chain, Protein, Residue};

const MAX_IDENTIFIER_BYTES: usize = 32;
pub const MAX_RESIDUE_COLORS: usize = 256;

/// A user-supplied request to color one exact polymer residue.
///
/// Residues are identified by chain ID, sequence number, and insertion code.
/// An omitted insertion code means the blank insertion code, not a wildcard.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResidueColorSpec {
    pub chain_id: String,
    pub residue_number: i32,
    pub insertion_code: Option<String>,
    pub color: [u8; 3],
}

impl ResidueColorSpec {
    pub fn new(
        chain_id: impl Into<String>,
        residue_number: i32,
        insertion_code: Option<String>,
        color: [u8; 3],
    ) -> Result<Self> {
        let chain_id = normalize_identifier(chain_id.into(), "chain ID")?;
        let insertion_code = insertion_code
            .map(|value| normalize_identifier(value, "insertion code"))
            .transpose()?
            .map(|value| value.to_ascii_uppercase());
        Ok(Self {
            chain_id,
            residue_number,
            insertion_code,
            color,
        })
    }
}

impl FromStr for ResidueColorSpec {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        let (selector, color) = value.rsplit_once('=').ok_or_else(|| {
            anyhow::anyhow!(
                "residue color must use CHAIN:RES[ICODE]=RRGGBB, for example A:42=FF0000"
            )
        })?;
        let (chain_id, residue) = selector.rsplit_once(':').ok_or_else(|| {
            anyhow::anyhow!(
                "residue color must use CHAIN:RES[ICODE]=RRGGBB, for example A:42=FF0000"
            )
        })?;
        let (residue_number, insertion_code) = parse_residue_id(residue)?;
        Self::new(
            chain_id.to_string(),
            residue_number,
            insertion_code,
            parse_rgb(color)?,
        )
    }
}

/// A validated residue-color entry, normalized against the parsed structure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedResidueColor {
    pub chain_id: String,
    pub residue_number: i32,
    pub insertion_code: Option<String>,
    pub residue_name: String,
    pub color: [u8; 3],
}

/// Validated exact-residue color overrides with allocation-free render lookup.
#[derive(Debug, Clone, Default)]
pub struct ResidueColorOverrides {
    entries: Vec<ResolvedResidueColor>,
    by_chain: HashMap<String, HashMap<i32, HashMap<String, [u8; 3]>>>,
}

impl ResidueColorOverrides {
    pub fn entries(&self) -> &[ResolvedResidueColor] {
        &self.entries
    }

    pub fn color_for(&self, residue: &Residue, chain: &Chain) -> Option<[u8; 3]> {
        self.by_chain
            .get(chain.id.as_str())?
            .get(&residue.seq_num)?
            .get(residue.insertion_code.as_deref().unwrap_or(""))
            .copied()
    }
}

/// Resolve all requested residues before rendering.
///
/// Unknown or duplicate selectors are rejected as a unit, so callers never
/// apply only part of a batch.
pub fn resolve_residue_colors(
    protein: &Protein,
    specs: &[ResidueColorSpec],
) -> Result<ResidueColorOverrides> {
    if specs.len() > MAX_RESIDUE_COLORS {
        bail!("at most {MAX_RESIDUE_COLORS} exact residue colors are allowed");
    }
    let mut entries = Vec::with_capacity(specs.len());
    let mut by_chain: HashMap<String, HashMap<i32, HashMap<String, [u8; 3]>>> = HashMap::new();
    let mut seen = HashSet::with_capacity(specs.len());

    for spec in specs {
        let chain = protein
            .chains
            .iter()
            .find(|chain| chain.id == spec.chain_id)
            .ok_or_else(|| {
                let available = protein
                    .chains
                    .iter()
                    .map(|chain| chain.id.as_str())
                    .collect::<Vec<_>>()
                    .join(", ");
                anyhow::anyhow!(
                    "residue color chain {:?} was not found; available chains: {}",
                    spec.chain_id,
                    available
                )
            })?;
        let residue = chain
            .residues
            .iter()
            .find(|residue| {
                residue.seq_num == spec.residue_number
                    && residue.insertion_code == spec.insertion_code
            })
            .ok_or_else(|| {
                let requested = format_residue_id(spec.residue_number, spec.insertion_code.as_deref());
                let alternatives = chain
                    .residues
                    .iter()
                    .filter(|residue| residue.seq_num == spec.residue_number)
                    .map(|residue| {
                        format_residue_id(residue.seq_num, residue.insertion_code.as_deref())
                    })
                    .collect::<Vec<_>>();
                if alternatives.is_empty() {
                    anyhow::anyhow!(
                        "residue color target {}:{} was not found",
                        spec.chain_id,
                        requested
                    )
                } else {
                    anyhow::anyhow!(
                        "residue color target {}:{} was not found; matching sequence number uses: {}",
                        spec.chain_id,
                        requested,
                        alternatives.join(", ")
                    )
                }
            })?;

        let key = (
            chain.id.clone(),
            residue.seq_num,
            residue.insertion_code.clone(),
        );
        if !seen.insert(key) {
            bail!(
                "duplicate residue color target {}:{}",
                chain.id,
                format_residue_id(residue.seq_num, residue.insertion_code.as_deref())
            );
        }

        let insertion_key = residue.insertion_code.clone().unwrap_or_default();
        by_chain
            .entry(chain.id.clone())
            .or_default()
            .entry(residue.seq_num)
            .or_default()
            .insert(insertion_key, spec.color);
        entries.push(ResolvedResidueColor {
            chain_id: chain.id.clone(),
            residue_number: residue.seq_num,
            insertion_code: residue.insertion_code.clone(),
            residue_name: residue.name.clone(),
            color: spec.color,
        });
    }

    Ok(ResidueColorOverrides { entries, by_chain })
}

pub fn parse_rgb(value: &str) -> Result<[u8; 3]> {
    if value.len() != 6
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'A'..=b'F').contains(&byte))
    {
        bail!(
            "residue color RGB must be exactly six uppercase hexadecimal digits, for example FF0000"
        );
    }
    Ok([
        u8::from_str_radix(&value[0..2], 16)?,
        u8::from_str_radix(&value[2..4], 16)?,
        u8::from_str_radix(&value[4..6], 16)?,
    ])
}

pub fn format_rgb(color: [u8; 3]) -> String {
    format!("{:02X}{:02X}{:02X}", color[0], color[1], color[2])
}

fn parse_residue_id(value: &str) -> Result<(i32, Option<String>)> {
    if value.is_empty() {
        bail!("residue number cannot be empty");
    }
    let bytes = value.as_bytes();
    let mut end = usize::from(matches!(bytes.first(), Some(b'-') | Some(b'+')));
    let digit_start = end;
    while bytes.get(end).is_some_and(u8::is_ascii_digit) {
        end += 1;
    }
    if end == digit_start {
        bail!("residue number must begin with a signed 32-bit integer");
    }
    let residue_number = value[..end]
        .parse::<i32>()
        .map_err(|_| anyhow::anyhow!("residue number must be a signed 32-bit integer"))?;
    let suffix = &value[end..];
    let insertion_code = if suffix.is_empty() {
        None
    } else if suffix.starts_with('[') || suffix.ends_with(']') {
        if !(suffix.starts_with('[') && suffix.ends_with(']') && suffix.len() > 2) {
            bail!("bracketed insertion code must use RES[ICODE], for example 42[A]");
        }
        Some(suffix[1..suffix.len() - 1].to_string())
    } else {
        Some(suffix.to_string())
    };
    Ok((residue_number, insertion_code))
}

fn normalize_identifier(value: String, label: &str) -> Result<String> {
    let candidate = value.trim();
    if candidate.is_empty()
        || candidate.len() > MAX_IDENTIFIER_BYTES
        || !candidate.bytes().all(|byte| byte.is_ascii_graphic())
    {
        bail!("{label} must contain 1-{MAX_IDENTIFIER_BYTES} printable non-whitespace ASCII bytes");
    }
    Ok(candidate.to_string())
}

pub(crate) fn format_residue_id(number: i32, insertion_code: Option<&str>) -> String {
    match insertion_code {
        Some(code) => format!("{number}[{code}]"),
        None => number.to_string(),
    }
}

impl fmt::Display for ResolvedResidueColor {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{}:{}={}",
            self.chain_id,
            format_residue_id(self.residue_number, self.insertion_code.as_deref()),
            format_rgb(self.color)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::protein::{MoleculeType, SecondaryStructure};

    fn residue(number: i32, insertion_code: Option<&str>) -> Residue {
        Residue {
            name: "ALA".to_string(),
            seq_num: number,
            insertion_code: insertion_code.map(str::to_string),
            atoms: Vec::new(),
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    fn protein() -> Protein {
        Protein {
            name: "selection".to_string(),
            chains: vec![Chain {
                id: "AA".to_string(),
                residues: vec![residue(-2, None), residue(42, None), residue(42, Some("A"))],
                molecule_type: MoleculeType::Protein,
            }],
            ligands: Vec::new(),
        }
    }

    #[test]
    fn parses_exact_residue_and_rgb() {
        let spec: ResidueColorSpec = "AA:-2=0A7FFF".parse().unwrap();
        assert_eq!(spec.chain_id, "AA");
        assert_eq!(spec.residue_number, -2);
        assert_eq!(spec.insertion_code, None);
        assert_eq!(spec.color, [0x0a, 0x7f, 0xff]);

        let inserted: ResidueColorSpec = "AA:42[A]=FF0000".parse().unwrap();
        assert_eq!(inserted.insertion_code.as_deref(), Some("A"));
    }

    #[test]
    fn rejects_malformed_selector_and_rgb() {
        assert!("A42=FF0000".parse::<ResidueColorSpec>().is_err());
        assert!("A:x=FF0000".parse::<ResidueColorSpec>().is_err());
        assert!("A:42=#FF0000".parse::<ResidueColorSpec>().is_err());
        assert!("A:42=ff0000".parse::<ResidueColorSpec>().is_err());
        assert!("A:42=FFFF".parse::<ResidueColorSpec>().is_err());
        assert!("A:42[]=FF0000".parse::<ResidueColorSpec>().is_err());
        assert!("A:42[A=FF0000".parse::<ResidueColorSpec>().is_err());
    }

    #[test]
    fn resolves_blank_and_inserted_residues_independently() {
        let specs = [
            "AA:42=FF0000".parse().unwrap(),
            "AA:42A=00FF00".parse().unwrap(),
        ];
        let protein = protein();
        let overrides = resolve_residue_colors(&protein, &specs).unwrap();
        assert_eq!(
            overrides.color_for(&protein.chains[0].residues[1], &protein.chains[0]),
            Some([255, 0, 0])
        );
        assert_eq!(
            overrides.color_for(&protein.chains[0].residues[2], &protein.chains[0]),
            Some([0, 255, 0])
        );
    }

    #[test]
    fn formats_insertion_codes_without_selector_ambiguity() {
        let resolved = ResolvedResidueColor {
            chain_id: "AA".to_string(),
            residue_number: 4,
            insertion_code: Some("2".to_string()),
            residue_name: "ALA".to_string(),
            color: [0xff, 0, 0],
        };

        assert_eq!(resolved.to_string(), "AA:4[2]=FF0000");
    }

    #[test]
    fn rejects_unknown_and_duplicate_targets_atomically() {
        let protein = protein();
        let missing = ["AA:7=FF0000".parse().unwrap()];
        assert!(resolve_residue_colors(&protein, &missing).is_err());

        let duplicate = [
            "AA:42A=FF0000".parse().unwrap(),
            "AA:42[A]=00FF00".parse().unwrap(),
        ];
        assert!(resolve_residue_colors(&protein, &duplicate).is_err());
    }

    #[test]
    fn enforces_the_bounded_exact_residue_batch() {
        let residues = (0..=MAX_RESIDUE_COLORS)
            .map(|number| residue(number as i32, None))
            .collect::<Vec<_>>();
        let protein = Protein {
            name: "bounded-selection".to_string(),
            chains: vec![Chain {
                id: "AA".to_string(),
                residues,
                molecule_type: MoleculeType::Protein,
            }],
            ligands: Vec::new(),
        };
        let specs = (0..=MAX_RESIDUE_COLORS)
            .map(|number| ResidueColorSpec::new("AA", number as i32, None, [255, 0, 0]).unwrap())
            .collect::<Vec<_>>();

        assert!(resolve_residue_colors(&protein, &specs[..MAX_RESIDUE_COLORS]).is_ok());
        assert!(resolve_residue_colors(&protein, &specs).is_err());
    }
}

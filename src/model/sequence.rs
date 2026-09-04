//! One-letter sequence codes and the wrapped layout of the sequence panel.
//!
//! The layout is computed once per panel width and then used by both the
//! renderer and the cursor navigation, so what the user sees and what the
//! arrow keys move through can never disagree.

use crate::model::protein::{MoleculeType, Protein};

/// Residues per group, separated by a single space in the panel.
pub const GROUP: usize = 10;

/// One row of the sequence panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqRow {
    /// Chain title row: `> Lb  Protein  138 res  1-138`.
    Header(usize),
    /// A wrapped run of residues from one chain.
    Residues {
        chain: usize,
        /// Index of the first residue of the run within the chain.
        start: usize,
        /// Number of residues on this row (`<= wrap`).
        len: usize,
    },
}

/// The full scrollable row list for one panel width.
#[derive(Debug, Clone, Default)]
pub struct SequenceLayout {
    /// Residues per full row.
    pub wrap: usize,
    /// Panel width this layout was built for.
    pub width: u16,
    pub rows: Vec<SeqRow>,
    /// Row index of each chain's header, so a cursor position maps to a row
    /// by arithmetic instead of a search.
    chain_header_row: Vec<usize>,
    /// Residue count per chain, cached for clamping.
    chain_len: Vec<usize>,
}

impl SequenceLayout {
    /// Build the layout for `protein` at a given residues-per-row.
    ///
    /// Empty chains still get a header row: a chain that parsed to zero
    /// polymer residues is information, not something to hide.
    pub fn build(protein: &Protein, wrap: usize, width: u16) -> Self {
        let wrap = wrap.max(1);
        let mut rows = Vec::new();
        let mut chain_header_row = Vec::with_capacity(protein.chains.len());
        let mut chain_len = Vec::with_capacity(protein.chains.len());

        for (chain_index, chain) in protein.chains.iter().enumerate() {
            chain_header_row.push(rows.len());
            chain_len.push(chain.residues.len());
            rows.push(SeqRow::Header(chain_index));
            let mut start = 0;
            while start < chain.residues.len() {
                let len = wrap.min(chain.residues.len() - start);
                rows.push(SeqRow::Residues {
                    chain: chain_index,
                    start,
                    len,
                });
                start += len;
            }
        }

        Self {
            wrap,
            width,
            rows,
            chain_header_row,
            chain_len,
        }
    }

    pub fn row_count(&self) -> usize {
        self.rows.len()
    }

    /// Row index holding `residue` of `chain`, and the column within it.
    pub fn locate(&self, chain: usize, residue: usize) -> Option<(usize, usize)> {
        let header = *self.chain_header_row.get(chain)?;
        if residue >= *self.chain_len.get(chain)? {
            return None;
        }
        Some((header + 1 + residue / self.wrap, residue % self.wrap))
    }

    /// Row index of a chain's header row.
    pub fn header_row(&self, chain: usize) -> Option<usize> {
        self.chain_header_row.get(chain).copied()
    }

    /// The residue at `column` of `row`, clamped to the row's last residue.
    ///
    /// Header rows have no residue, so vertical navigation skips them.
    pub fn residue_at(&self, row: usize, column: usize) -> Option<(usize, usize)> {
        match self.rows.get(row)? {
            SeqRow::Header(_) => None,
            SeqRow::Residues { chain, start, len } => {
                Some((*chain, start + column.min(len.saturating_sub(1))))
            }
        }
    }
}

/// Screen column offset of a residue column, accounting for group spacing.
#[inline]
pub fn column_offset(column: usize) -> usize {
    column + column / GROUP
}

/// Residues that fit in `avail` columns when every group of [`GROUP`] is
/// followed by a space.
pub fn wrap_for_width(avail: usize) -> usize {
    // A group costs GROUP characters plus one separator, except the last one
    // which needs no trailing space.
    let groups = (avail + 1) / (GROUP + 1);
    (groups * GROUP).max(GROUP)
}

/// One-letter code for a residue, given its chain's polymer type.
///
/// Modified residues map to their parent letter (`PSU` -> `U`, `MSE` -> `M`)
/// so a modified base never breaks the reading frame of the sequence; the
/// panel's cursor line always shows the true three-letter name.  Anything
/// unrecognized is `X`.
pub fn one_letter(name: &str, molecule_type: MoleculeType) -> char {
    let name = name.trim();
    match molecule_type {
        MoleculeType::RNA | MoleculeType::DNA => nucleotide_letter(name),
        MoleculeType::Protein | MoleculeType::SmallMolecule => amino_acid_letter(name)
            // A chain classified as protein can still carry stray nucleotides
            // (hybrid or mis-classified chains); fall back before giving up.
            .or_else(|| {
                let letter = nucleotide_letter(name);
                (letter != 'X').then_some(letter)
            })
            .unwrap_or('X'),
    }
}

fn amino_acid_letter(name: &str) -> Option<char> {
    Some(match name {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" | "CYX" | "CYM" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" | "LYN" => 'K',
        "MET" | "MSE" | "FME" => 'M',
        "PHE" => 'F',
        "PRO" | "HYP" => 'P',
        "SER" | "SEP" => 'S',
        "THR" | "TPO" => 'T',
        "TRP" => 'W',
        "TYR" | "PTR" => 'Y',
        "VAL" => 'V',
        "SEC" => 'U',
        "PYL" => 'O',
        "ASX" => 'B',
        "GLX" => 'Z',
        "UNK" => 'X',
        _ => return None,
    })
}

fn nucleotide_letter(name: &str) -> char {
    match name {
        "A" | "DA" | "AMP" | "ADE" | "1MA" | "6MA" | "MA6" | "2MA" => 'A',
        "C" | "DC" | "CMP" | "CYT" | "5MC" | "OMC" | "4OC" | "3MC" => 'C',
        "G" | "DG" | "GMP" | "GUA" | "2MG" | "7MG" | "M2G" | "OMG" | "1MG" | "YG" => 'G',
        "U" | "UMP" | "URA" | "URI" | "PSU" | "4SU" | "H2U" | "5MU" | "UR3" | "OMU" | "3MU" => 'U',
        "T" | "DT" | "THY" | "5MT" => 'T',
        "I" | "DI" => 'I',
        "N" | "DN" => 'N',
        _ => 'X',
    }
}

/// Short label for a chain's polymer type, for the panel header.
pub fn molecule_label(molecule_type: MoleculeType) -> &'static str {
    match molecule_type {
        MoleculeType::Protein => "Protein",
        MoleculeType::RNA => "RNA",
        MoleculeType::DNA => "DNA",
        MoleculeType::SmallMolecule => "Other",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::protein::{Chain, Residue, SecondaryStructure};

    fn chain(id: &str, count: usize, molecule_type: MoleculeType) -> Chain {
        Chain {
            id: id.to_string(),
            molecule_type,
            residues: (0..count)
                .map(|i| Residue {
                    name: "ALA".to_string(),
                    seq_num: i as i32 + 1,
                    insertion_code: None,
                    atoms: Vec::new(),
                    secondary_structure: SecondaryStructure::Coil,
                })
                .collect(),
        }
    }

    fn protein(counts: &[usize]) -> Protein {
        Protein {
            name: "seq".to_string(),
            chains: counts
                .iter()
                .enumerate()
                .map(|(i, n)| chain(&format!("C{i}"), *n, MoleculeType::Protein))
                .collect(),
            ligands: Vec::new(),
        }
    }

    #[test]
    fn locate_and_residue_at_are_inverses() {
        let protein = protein(&[25, 3]);
        let layout = SequenceLayout::build(&protein, 10, 80);
        for chain in 0..2 {
            let len = protein.chains[chain].residues.len();
            for residue in 0..len {
                let (row, column) = layout.locate(chain, residue).unwrap();
                assert_eq!(layout.residue_at(row, column), Some((chain, residue)));
            }
        }
    }

    #[test]
    fn empty_chains_still_get_a_header() {
        let protein = protein(&[0, 4]);
        let layout = SequenceLayout::build(&protein, 10, 80);
        assert_eq!(layout.rows[0], SeqRow::Header(0));
        assert_eq!(layout.rows[1], SeqRow::Header(1));
        assert_eq!(layout.locate(0, 0), None);
    }

    #[test]
    fn short_rows_clamp_the_column() {
        let protein = protein(&[13]);
        let layout = SequenceLayout::build(&protein, 10, 80);
        // Second residue row holds 3 residues; column 7 clamps to the last.
        assert_eq!(layout.residue_at(2, 7), Some((0, 12)));
    }

    #[test]
    fn wrap_accounts_for_group_separators() {
        // Three groups need 10+1+10+1+10 = 32 columns; one fewer fits only two.
        assert_eq!(wrap_for_width(32), 30);
        assert_eq!(wrap_for_width(31), 20);
        // Never returns zero, however narrow the panel is.
        assert_eq!(wrap_for_width(0), GROUP);
    }

    #[test]
    fn modified_residues_keep_the_reading_frame() {
        assert_eq!(one_letter("PSU", MoleculeType::RNA), 'U');
        assert_eq!(one_letter("MSE", MoleculeType::Protein), 'M');
        assert_eq!(one_letter("T1C", MoleculeType::Protein), 'X');
        // A nucleotide inside a chain classified as protein still reads.
        assert_eq!(one_letter("G", MoleculeType::Protein), 'G');
    }
}

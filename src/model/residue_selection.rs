//! The set of residues the user has picked in the sequence panel.
//!
//! The renderer asks "is this residue selected?" once per residue per frame,
//! so membership is a flat `[chain][residue]` bitmap rather than a hash set:
//! lookup is an index, and the whole thing costs one byte per residue.

use crate::model::protein::Protein;
use crate::model::selection::format_residue_id;

/// How a selection should appear in the 3D view.
///
/// Bundled into one borrow so the renderers take a single extra argument, and
/// so a caller with no selection (the snapshot and panel-server paths) simply
/// passes `None`.
#[derive(Debug, Clone, Copy)]
pub struct SelectionView<'a> {
    pub selection: &'a ResidueSelection,
    /// Draw every atom of a picked residue, rather than one marker sphere.
    pub ball_and_stick: bool,
}

/// A per-residue selection over one parsed structure.
#[derive(Debug, Clone, Default)]
pub struct ResidueSelection {
    flags: Vec<Vec<bool>>,
    count: usize,
}

impl ResidueSelection {
    /// An empty selection shaped to `protein`.
    pub fn new(protein: &Protein) -> Self {
        Self {
            flags: protein
                .chains
                .iter()
                .map(|chain| vec![false; chain.residues.len()])
                .collect(),
            count: 0,
        }
    }

    #[inline]
    pub fn contains(&self, chain: usize, residue: usize) -> bool {
        self.flags
            .get(chain)
            .and_then(|c| c.get(residue))
            .copied()
            .unwrap_or(false)
    }

    /// Whether any residue of `chain` is selected.
    pub fn chain_has_any(&self, chain: usize) -> bool {
        self.flags
            .get(chain)
            .is_some_and(|c| c.iter().any(|selected| *selected))
    }

    pub fn set(&mut self, chain: usize, residue: usize, selected: bool) {
        let Some(slot) = self.flags.get_mut(chain).and_then(|c| c.get_mut(residue)) else {
            return;
        };
        if *slot == selected {
            return;
        }
        *slot = selected;
        if selected {
            self.count += 1;
        } else {
            self.count -= 1;
        }
    }

    pub fn toggle(&mut self, chain: usize, residue: usize) {
        let selected = self.contains(chain, residue);
        self.set(chain, residue, !selected);
    }

    /// Set an inclusive residue range within one chain.  The endpoints may be
    /// given in either order.
    pub fn set_range(&mut self, chain: usize, from: usize, to: usize, selected: bool) {
        let (lo, hi) = if from <= to { (from, to) } else { (to, from) };
        for residue in lo..=hi {
            self.set(chain, residue, selected);
        }
    }

    /// Select or deselect an entire chain.
    pub fn set_chain(&mut self, chain: usize, selected: bool) {
        let Some(len) = self.flags.get(chain).map(Vec::len) else {
            return;
        };
        if len > 0 {
            self.set_range(chain, 0, len - 1, selected);
        }
    }

    pub fn clear(&mut self) {
        for chain in &mut self.flags {
            chain.fill(false);
        }
        self.count = 0;
    }

    pub fn count(&self) -> usize {
        self.count
    }

    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Number of chains contributing at least one residue.
    pub fn chain_count(&self) -> usize {
        (0..self.flags.len())
            .filter(|chain| self.chain_has_any(*chain))
            .count()
    }

    /// Mean position of every atom in the selection, for centring the camera.
    pub fn centroid(&self, protein: &Protein) -> Option<[f64; 3]> {
        let mut sum = [0.0f64; 3];
        let mut n = 0usize;
        for (chain_index, chain) in protein.chains.iter().enumerate() {
            for (residue_index, residue) in chain.residues.iter().enumerate() {
                if !self.contains(chain_index, residue_index) {
                    continue;
                }
                for atom in &residue.atoms {
                    sum[0] += atom.x;
                    sum[1] += atom.y;
                    sum[2] += atom.z;
                    n += 1;
                }
            }
        }
        (n > 0).then(|| [sum[0] / n as f64, sum[1] / n as f64, sum[2] / n as f64])
    }

    /// Compact human-readable form: `Lb:12-18,40 LC:7`.
    ///
    /// Truncated after `max_chains` chains so it always fits on one line.
    pub fn describe(&self, protein: &Protein, max_chains: usize) -> String {
        let mut parts: Vec<String> = Vec::new();
        let mut skipped = 0usize;

        for (chain_index, chain) in protein.chains.iter().enumerate() {
            let ranges = self.chain_ranges(chain_index, chain.residues.len());
            if ranges.is_empty() {
                continue;
            }
            if parts.len() == max_chains {
                skipped += 1;
                continue;
            }
            let spans: Vec<String> = ranges
                .iter()
                .map(|(from, to)| {
                    let first = &chain.residues[*from];
                    let last = &chain.residues[*to];
                    let first_id =
                        format_residue_id(first.seq_num, first.insertion_code.as_deref());
                    if from == to {
                        first_id
                    } else {
                        let last_id =
                            format_residue_id(last.seq_num, last.insertion_code.as_deref());
                        format!("{first_id}-{last_id}")
                    }
                })
                .collect();
            parts.push(format!("{}:{}", chain.id, spans.join(",")));
        }

        if skipped > 0 {
            parts.push(format!("+{skipped} more"));
        }
        parts.join("  ")
    }

    /// Contiguous selected runs within one chain, as residue index pairs.
    fn chain_ranges(&self, chain: usize, len: usize) -> Vec<(usize, usize)> {
        let mut ranges = Vec::new();
        let mut start: Option<usize> = None;
        for residue in 0..len {
            match (self.contains(chain, residue), start) {
                (true, None) => start = Some(residue),
                (false, Some(from)) => {
                    ranges.push((from, residue - 1));
                    start = None;
                }
                _ => {}
            }
        }
        if let Some(from) = start {
            ranges.push((from, len - 1));
        }
        ranges
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::protein::{Atom, Chain, MoleculeType, Residue, SecondaryStructure};

    fn protein() -> Protein {
        let chain = |id: &str, count: usize, x0: f64| Chain {
            id: id.to_string(),
            molecule_type: MoleculeType::Protein,
            residues: (0..count)
                .map(|i| Residue {
                    name: "ALA".to_string(),
                    seq_num: i as i32 + 1,
                    insertion_code: None,
                    atoms: vec![Atom {
                        name: "CA".to_string(),
                        element: "C".to_string(),
                        x: x0 + i as f64,
                        y: 0.0,
                        z: 0.0,
                        b_factor: 0.0,
                        is_backbone: true,
                        is_hetero: false,
                    }],
                    secondary_structure: SecondaryStructure::Coil,
                })
                .collect(),
        };
        Protein {
            name: "sel".to_string(),
            chains: vec![chain("A", 10, 0.0), chain("B", 5, 100.0)],
            ligands: Vec::new(),
        }
    }

    #[test]
    fn count_tracks_set_toggle_and_clear() {
        let protein = protein();
        let mut selection = ResidueSelection::new(&protein);
        assert!(selection.is_empty());

        selection.set_range(0, 2, 5, true);
        assert_eq!(selection.count(), 4);
        // Re-selecting an already selected residue must not double count.
        selection.set(0, 3, true);
        assert_eq!(selection.count(), 4);

        selection.toggle(0, 3);
        assert_eq!(selection.count(), 3);
        assert!(!selection.contains(0, 3));

        selection.set_chain(1, true);
        assert_eq!(selection.count(), 8);
        assert_eq!(selection.chain_count(), 2);

        selection.clear();
        assert!(selection.is_empty());
        assert_eq!(selection.chain_count(), 0);
    }

    #[test]
    fn out_of_range_indices_are_ignored() {
        let protein = protein();
        let mut selection = ResidueSelection::new(&protein);
        selection.set(9, 0, true);
        selection.set(0, 999, true);
        assert!(selection.is_empty());
        assert!(!selection.contains(9, 0));
    }

    #[test]
    fn describe_collapses_runs_and_reports_residue_numbers() {
        let protein = protein();
        let mut selection = ResidueSelection::new(&protein);
        selection.set_range(0, 1, 3, true);
        selection.set(0, 7, true);
        selection.set(1, 0, true);
        assert_eq!(selection.describe(&protein, 8), "A:2-4,8  B:1");
        assert_eq!(selection.describe(&protein, 1), "A:2-4,8  +1 more");
    }

    #[test]
    fn centroid_averages_selected_atoms_only() {
        let protein = protein();
        let mut selection = ResidueSelection::new(&protein);
        assert!(selection.centroid(&protein).is_none());
        selection.set_range(0, 0, 2, true);
        let centroid = selection.centroid(&protein).unwrap();
        assert!((centroid[0] - 1.0).abs() < 1e-9);
    }
}

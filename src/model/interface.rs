use std::collections::HashSet;

use crate::model::protein::Protein;

/// A contact between two residues on different chains.
#[derive(Debug, Clone)]
pub struct Contact {
    /// Index into `protein.chains` for the first residue's chain.
    pub chain_a: usize,
    /// Index into `chain.residues` for the first residue.
    pub residue_a: usize,
    /// Index into `protein.chains` for the second residue's chain.
    pub chain_b: usize,
    /// Index into `chain.residues` for the second residue.
    pub residue_b: usize,
    /// Minimum heavy-atom distance between the two residues in Angstroms.
    pub min_distance: f64,
}

/// Full interface analysis result.
#[derive(Debug, Clone)]
pub struct InterfaceAnalysis {
    /// All inter-chain residue-residue contacts.
    pub contacts: Vec<Contact>,
    /// Set of (chain_idx, residue_idx) pairs that lie at the interface.
    pub interface_residues: HashSet<(usize, usize)>,
    /// Per-chain count of interface residues (indexed by chain position).
    pub chain_interface_counts: Vec<usize>,
    /// Total number of unique interface residues across all chains.
    pub total_interface_residues: usize,
}

/// Squared Euclidean distance between two atoms.
#[inline]
fn dist_sq(ax: f64, ay: f64, az: f64, bx: f64, by: f64, bz: f64) -> f64 {
    let dx = ax - bx;
    let dy = ay - by;
    let dz = az - bz;
    dx * dx + dy * dy + dz * dz
}

/// Analyze the interface between all chain pairs in a protein.
///
/// A contact exists when any heavy atom (non-hydrogen) of residue A is within
/// `cutoff` Angstroms of any heavy atom of residue B, where A and B belong to
/// different chains.
///
/// The default cutoff used in most structural biology tools is 4.5 A.
pub fn analyze_interface(protein: &Protein, cutoff: f64) -> InterfaceAnalysis {
    let cutoff_sq = cutoff * cutoff;
    let num_chains = protein.chains.len();

    let mut contacts: Vec<Contact> = Vec::new();
    let mut interface_residues: HashSet<(usize, usize)> = HashSet::new();

    // Compare every pair of chains (i < j).
    for i in 0..num_chains {
        for j in (i + 1)..num_chains {
            let chain_i = &protein.chains[i];
            let chain_j = &protein.chains[j];

            for (ri, res_i) in chain_i.residues.iter().enumerate() {
                for (rj, res_j) in chain_j.residues.iter().enumerate() {
                    let mut min_d_sq = f64::MAX;
                    let mut found_contact = false;

                    // Compare all heavy-atom pairs between the two residues.
                    for atom_a in &res_i.atoms {
                        if atom_a.element == "H" {
                            continue;
                        }
                        for atom_b in &res_j.atoms {
                            if atom_b.element == "H" {
                                continue;
                            }
                            let d_sq =
                                dist_sq(atom_a.x, atom_a.y, atom_a.z, atom_b.x, atom_b.y, atom_b.z);
                            if d_sq < min_d_sq {
                                min_d_sq = d_sq;
                            }
                            if d_sq <= cutoff_sq {
                                found_contact = true;
                            }
                        }
                    }

                    if found_contact {
                        contacts.push(Contact {
                            chain_a: i,
                            residue_a: ri,
                            chain_b: j,
                            residue_b: rj,
                            min_distance: min_d_sq.sqrt(),
                        });
                        interface_residues.insert((i, ri));
                        interface_residues.insert((j, rj));
                    }
                }
            }
        }
    }

    // Sort contacts by minimum distance (closest first).
    contacts.sort_by(|a, b| a.min_distance.partial_cmp(&b.min_distance).unwrap_or(std::cmp::Ordering::Equal));

    // Count interface residues per chain.
    let mut chain_interface_counts = vec![0usize; num_chains];
    for &(chain_idx, _) in &interface_residues {
        chain_interface_counts[chain_idx] += 1;
    }

    let total_interface_residues = interface_residues.len();

    InterfaceAnalysis {
        contacts,
        interface_residues,
        chain_interface_counts,
        total_interface_residues,
    }
}

impl InterfaceAnalysis {
    /// Convert interface residues to (chain_id, seq_num) pairs using the protein.
    pub fn interface_residues_by_id_with_protein(&self, protein: &Protein) -> HashSet<(String, i32)> {
        let mut set = HashSet::new();
        for &(chain_idx, res_idx) in &self.interface_residues {
            if let Some(chain) = protein.chains.get(chain_idx) {
                if let Some(residue) = chain.residues.get(res_idx) {
                    set.insert((chain.id.clone(), residue.seq_num));
                }
            }
        }
        set
    }

    /// Return all contacts between a specific pair of chains.
    pub fn contacts_between(&self, chain_a: usize, chain_b: usize) -> Vec<&Contact> {
        self.contacts
            .iter()
            .filter(|c| {
                (c.chain_a == chain_a && c.chain_b == chain_b)
                    || (c.chain_a == chain_b && c.chain_b == chain_a)
            })
            .collect()
    }

    /// Produce human-readable summary lines suitable for TUI display.
    ///
    /// Format:
    /// ```text
    /// Interface: 24 residues (Chain A: 12, Chain B: 12)
    /// Chain A-B: 18 contacts, min dist 2.8A
    /// Top contacts: A:ARG45-B:ASP102 (2.8A), A:TYR32-B:GLU156 (3.1A), ...
    /// ```
    pub fn summary(&self, protein: &Protein) -> Vec<String> {
        let mut lines: Vec<String> = Vec::new();

        if self.contacts.is_empty() {
            lines.push("Interface: no inter-chain contacts detected".to_string());
            return lines;
        }

        // Line 1 -- overall residue counts.
        let per_chain: Vec<String> = protein
            .chains
            .iter()
            .enumerate()
            .filter(|(idx, _)| self.chain_interface_counts.get(*idx).copied().unwrap_or(0) > 0)
            .map(|(idx, chain)| {
                format!(
                    "Chain {}: {}",
                    chain.id,
                    self.chain_interface_counts[idx]
                )
            })
            .collect();

        lines.push(format!(
            "Interface: {} residues ({})",
            self.total_interface_residues,
            per_chain.join(", ")
        ));

        // Lines 2..N -- per chain-pair statistics.
        let num_chains = protein.chains.len();
        for i in 0..num_chains {
            for j in (i + 1)..num_chains {
                let pair_contacts = self.contacts_between(i, j);
                if pair_contacts.is_empty() {
                    continue;
                }
                let min_dist = pair_contacts
                    .iter()
                    .map(|c| c.min_distance)
                    .fold(f64::MAX, f64::min);

                lines.push(format!(
                    "Chain {}-{}: {} contacts, min dist {:.1}\u{00C5}",
                    protein.chains[i].id,
                    protein.chains[j].id,
                    pair_contacts.len(),
                    min_dist,
                ));

                // Top 5 closest contacts for this pair.
                let top_n = 5.min(pair_contacts.len());
                // pair_contacts are already sorted by distance (inherited from
                // the globally-sorted contacts vec) but a local sort is cheap
                // and guarantees correctness for the filtered subset.
                let mut sorted: Vec<&&Contact> = pair_contacts.iter().collect();
                sorted.sort_by(|a, b| {
                    a.min_distance
                        .partial_cmp(&b.min_distance)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });

                let labels: Vec<String> = sorted[..top_n]
                    .iter()
                    .map(|c| {
                        let res_a = &protein.chains[c.chain_a].residues[c.residue_a];
                        let res_b = &protein.chains[c.chain_b].residues[c.residue_b];
                        format!(
                            "{}:{}{}-{}:{}{} ({:.1}\u{00C5})",
                            protein.chains[c.chain_a].id,
                            res_a.name,
                            res_a.seq_num,
                            protein.chains[c.chain_b].id,
                            res_b.name,
                            res_b.seq_num,
                            c.min_distance,
                        )
                    })
                    .collect();

                lines.push(format!("Top contacts: {}", labels.join(", ")));
            }
        }

        lines
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::protein::{Atom, Chain, MoleculeType, Protein, Residue, SecondaryStructure};

    /// Helper: make a single-atom residue at the given position.
    fn make_residue(name: &str, seq_num: i32, x: f64, y: f64, z: f64) -> Residue {
        Residue {
            name: name.to_string(),
            seq_num,
            atoms: vec![Atom {
                name: "CA".to_string(),
                element: "C".to_string(),
                x,
                y,
                z,
                b_factor: 0.0,
                is_backbone: true,
            }],
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    fn two_chain_protein() -> Protein {
        Protein {
            name: "test".to_string(),
            chains: vec![
                Chain {
                    id: "A".to_string(),
                    residues: vec![
                        make_residue("ALA", 1, 0.0, 0.0, 0.0),
                        make_residue("GLY", 2, 10.0, 0.0, 0.0), // far away
                    ],
                    molecule_type: MoleculeType::Protein,
                },
                Chain {
                    id: "B".to_string(),
                    residues: vec![
                        make_residue("ASP", 1, 3.0, 0.0, 0.0), // within 4.5 of A:ALA1
                        make_residue("LEU", 2, 20.0, 0.0, 0.0), // far away
                    ],
                    molecule_type: MoleculeType::Protein,
                },
            ],
        }
    }

    #[test]
    fn test_contact_detected() {
        let protein = two_chain_protein();
        let analysis = analyze_interface(&protein, 4.5);

        assert_eq!(analysis.contacts.len(), 1);
        assert_eq!(analysis.total_interface_residues, 2);

        let c = &analysis.contacts[0];
        assert_eq!(c.chain_a, 0);
        assert_eq!(c.residue_a, 0);
        assert_eq!(c.chain_b, 1);
        assert_eq!(c.residue_b, 0);
        assert!((c.min_distance - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_is_interface_residue() {
        let protein = two_chain_protein();
        let analysis = analyze_interface(&protein, 4.5);

        assert!(analysis.interface_residues.contains(&(0, 0)));
        assert!(analysis.interface_residues.contains(&(1, 0)));
        assert!(!analysis.interface_residues.contains(&(0, 1)));
        assert!(!analysis.interface_residues.contains(&(1, 1)));
    }

    #[test]
    fn test_contacts_between() {
        let protein = two_chain_protein();
        let analysis = analyze_interface(&protein, 4.5);

        let ab = analysis.contacts_between(0, 1);
        assert_eq!(ab.len(), 1);

        // Reversed order should also work.
        let ba = analysis.contacts_between(1, 0);
        assert_eq!(ba.len(), 1);
    }

    #[test]
    fn test_no_contacts_below_cutoff() {
        let protein = two_chain_protein();
        let analysis = analyze_interface(&protein, 2.0);

        assert!(analysis.contacts.is_empty());
        assert_eq!(analysis.total_interface_residues, 0);
    }

    #[test]
    fn test_hydrogen_atoms_skipped() {
        let protein = Protein {
            name: "htest".to_string(),
            chains: vec![
                Chain {
                    id: "A".to_string(),
                    residues: vec![Residue {
                        name: "ALA".to_string(),
                        seq_num: 1,
                        atoms: vec![Atom {
                            name: "H".to_string(),
                            element: "H".to_string(),
                            x: 0.0,
                            y: 0.0,
                            z: 0.0,
                            b_factor: 0.0,
                            is_backbone: false,
                        }],
                        secondary_structure: SecondaryStructure::Coil,
                    }],
                    molecule_type: MoleculeType::Protein,
                },
                Chain {
                    id: "B".to_string(),
                    residues: vec![Residue {
                        name: "ASP".to_string(),
                        seq_num: 1,
                        atoms: vec![Atom {
                            name: "H".to_string(),
                            element: "H".to_string(),
                            x: 1.0,
                            y: 0.0,
                            z: 0.0,
                            b_factor: 0.0,
                            is_backbone: false,
                        }],
                        secondary_structure: SecondaryStructure::Coil,
                    }],
                    molecule_type: MoleculeType::Protein,
                },
            ],
        };

        let analysis = analyze_interface(&protein, 4.5);
        assert!(analysis.contacts.is_empty(), "hydrogen-only atoms must not produce contacts");
    }

    #[test]
    fn test_summary_format() {
        let protein = two_chain_protein();
        let analysis = analyze_interface(&protein, 4.5);
        let lines = analysis.summary(&protein);

        assert!(!lines.is_empty());
        assert!(lines[0].starts_with("Interface:"));
        assert!(lines[0].contains("2 residues"));
    }

    #[test]
    fn test_chain_interface_counts() {
        let protein = two_chain_protein();
        let analysis = analyze_interface(&protein, 4.5);

        assert_eq!(analysis.chain_interface_counts[0], 1);
        assert_eq!(analysis.chain_interface_counts[1], 1);
    }

    #[test]
    fn test_empty_protein() {
        let protein = Protein {
            name: "empty".to_string(),
            chains: vec![],
        };
        let analysis = analyze_interface(&protein, 4.5);

        assert!(analysis.contacts.is_empty());
        assert_eq!(analysis.total_interface_residues, 0);
        let lines = analysis.summary(&protein);
        assert_eq!(lines.len(), 1);
        assert!(lines[0].contains("no inter-chain contacts"));
    }
}

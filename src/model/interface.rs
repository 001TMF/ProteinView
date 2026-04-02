use std::collections::HashSet;

use crate::model::protein::Protein;

/// Classification of inter-residue interactions at a chain interface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InteractionType {
    HydrogenBond,
    SaltBridge,
    HydrophobicContact,
    Other,
}

/// A classified interaction between two atoms across a chain interface.
#[derive(Debug, Clone)]
pub struct Interaction {
    pub interaction_type: InteractionType,
    pub atom_a: [f64; 3],
    pub atom_b: [f64; 3],
    /// Distance in angstroms between the two atoms. Stored for future
    /// display in tooltips or filtering by distance threshold.
    #[allow(dead_code)]
    pub distance: f64,
}

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

/// A contact between a ligand and a polymer residue.
#[derive(Debug, Clone)]
pub struct LigandContact {
    pub ligand_idx: usize,
    pub chain_idx: usize,
    pub residue_idx: usize,
    pub min_distance: f64,
}

/// Binding pocket analysis for all ligands.
#[derive(Debug, Clone)]
pub struct BindingPocketAnalysis {
    pub contacts: Vec<LigandContact>,
    /// Per-ligand: set of (chain_idx, residue_idx) forming the binding pocket.
    pub pockets: Vec<HashSet<(usize, usize)>>,
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
    pub binding_pockets: Option<BindingPocketAnalysis>,
    /// Classified atom-atom interactions across chain interfaces.
    pub interactions: Vec<Interaction>,
}

/// Squared Euclidean distance between two atoms.
#[inline]
fn dist_sq(ax: f64, ay: f64, az: f64, bx: f64, by: f64, bz: f64) -> f64 {
    let dx = ax - bx;
    let dy = ay - by;
    let dz = az - bz;
    dx * dx + dy * dy + dz * dz
}

/// Returns true if the atom is part of the positively-charged group
/// (ARG guanidinium NE/NH1/NH2, LYS ammonium NZ).
fn is_charged_positive(res_name: &str, atom_name: &str) -> bool {
    match res_name {
        "ARG" => atom_name.starts_with("NH") || atom_name == "NE",
        "LYS" => atom_name == "NZ",
        _ => false,
    }
}

/// Returns true if the atom on a negatively-charged residue side-chain
/// carries formal negative charge (ASP OD*, GLU OE*).
fn is_charged_negative(res_name: &str, atom_name: &str) -> bool {
    match res_name {
        "ASP" => atom_name.starts_with("OD"),
        "GLU" => atom_name.starts_with("OE"),
        _ => false,
    }
}

/// Returns true if the residue is typically hydrophobic.
fn is_hydrophobic_residue(name: &str) -> bool {
    matches!(
        name,
        "ALA" | "VAL" | "LEU" | "ILE" | "PHE" | "TRP" | "MET" | "PRO"
    )
}

/// Classify inter-chain interactions from the given contacts.
///
/// For each [`Contact`] the function finds the closest heavy-atom pair between
/// the two residues and classifies the interaction by distance and chemistry.
/// Classification uses a priority order — the first matching rule wins:
///
///   1. **Salt bridge** (distance <= 4.0 A): a positively-charged sidechain
///      atom (ARG NE/NH1/NH2, LYS NZ) paired with a negatively-charged
///      sidechain atom (ASP OD*, GLU OE*), checked symmetrically.
///   2. **H-bond** (distance <= 3.5 A): both atoms drawn from {N, O, S},
///      with at least one being N or O.
///   3. **Hydrophobic contact** (distance <= 4.5 A): C–C pair where both
///      residues are hydrophobic (ALA, VAL, LEU, ILE, PHE, TRP, MET, PRO).
///   4. **Other**: anything not matching the above.
///
/// **Limitation:** classification is based on the single closest heavy-atom
/// pair per contact, which may miss secondary interactions (e.g. a salt bridge
/// at 3.8 A when a closer C–C pair exists at 3.5 A).
fn classify_interactions(protein: &Protein, contacts: &[Contact]) -> Vec<Interaction> {
    let mut interactions = Vec::with_capacity(contacts.len());

    for contact in contacts {
        let Some(chain_a) = protein.chains.get(contact.chain_a) else {
            continue;
        };
        let Some(chain_b) = protein.chains.get(contact.chain_b) else {
            continue;
        };
        let Some(res_a) = chain_a.residues.get(contact.residue_a) else {
            continue;
        };
        let Some(res_b) = chain_b.residues.get(contact.residue_b) else {
            continue;
        };

        // Find the closest heavy-atom pair.
        let mut best_dist_sq = f64::MAX;
        let mut best_a: Option<&crate::model::protein::Atom> = None;
        let mut best_b: Option<&crate::model::protein::Atom> = None;

        for atom_a in &res_a.atoms {
            if atom_a.element.trim() == "H" {
                continue;
            }
            for atom_b in &res_b.atoms {
                if atom_b.element.trim() == "H" {
                    continue;
                }
                let d_sq = dist_sq(atom_a.x, atom_a.y, atom_a.z, atom_b.x, atom_b.y, atom_b.z);
                if d_sq < best_dist_sq {
                    best_dist_sq = d_sq;
                    best_a = Some(atom_a);
                    best_b = Some(atom_b);
                }
            }
        }

        let (atom_a, atom_b) = match (best_a, best_b) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };

        let distance = best_dist_sq.sqrt();

        // Classify
        let interaction_type = if distance <= 4.0
            && ((is_charged_positive(&res_a.name, &atom_a.name)
                && is_charged_negative(&res_b.name, &atom_b.name))
                || (is_charged_negative(&res_a.name, &atom_a.name)
                    && is_charged_positive(&res_b.name, &atom_b.name)))
        {
            InteractionType::SaltBridge
        } else if distance <= 3.5 {
            let el_a = atom_a.element.as_str();
            let el_b = atom_b.element.as_str();
            let donor_acceptor = matches!(el_a, "N" | "O") && matches!(el_b, "N" | "O" | "S");
            let acceptor_donor = matches!(el_b, "N" | "O") && matches!(el_a, "N" | "O" | "S");
            if donor_acceptor || acceptor_donor {
                InteractionType::HydrogenBond
            } else if el_a == "C"
                && el_b == "C"
                && is_hydrophobic_residue(&res_a.name)
                && is_hydrophobic_residue(&res_b.name)
            {
                InteractionType::HydrophobicContact
            } else {
                InteractionType::Other
            }
        } else if distance <= 4.5
            && atom_a.element == "C"
            && atom_b.element == "C"
            && is_hydrophobic_residue(&res_a.name)
            && is_hydrophobic_residue(&res_b.name)
        {
            InteractionType::HydrophobicContact
        } else {
            InteractionType::Other
        };

        interactions.push(Interaction {
            interaction_type,
            atom_a: [atom_a.x, atom_a.y, atom_a.z],
            atom_b: [atom_b.x, atom_b.y, atom_b.z],
            distance,
        });
    }

    interactions
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
                        if atom_a.element.trim() == "H" {
                            continue;
                        }
                        for atom_b in &res_j.atoms {
                            if atom_b.element.trim() == "H" {
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
                        let min_distance = min_d_sq.sqrt();
                        contacts.push(Contact {
                            chain_a: i,
                            residue_a: ri,
                            chain_b: j,
                            residue_b: rj,
                            min_distance,
                        });
                        interface_residues.insert((i, ri));
                        interface_residues.insert((j, rj));
                    }
                }
            }
        }
    }

    // Sort contacts by minimum distance (closest first).
    contacts.sort_by(|a, b| {
        a.min_distance
            .partial_cmp(&b.min_distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Count interface residues per chain.
    let mut chain_interface_counts = vec![0usize; num_chains];
    for &(chain_idx, _) in &interface_residues {
        chain_interface_counts[chain_idx] += 1;
    }

    let total_interface_residues = interface_residues.len();

    let interactions = classify_interactions(protein, &contacts);

    InterfaceAnalysis {
        contacts,
        interface_residues,
        chain_interface_counts,
        total_interface_residues,
        binding_pockets: None,
        interactions,
    }
}

/// Analyze ligand binding pockets: find polymer residues within cutoff of each ligand.
pub fn analyze_binding_pockets(protein: &Protein, cutoff: f64) -> BindingPocketAnalysis {
    let cutoff_sq = cutoff * cutoff;
    let mut contacts: Vec<LigandContact> = Vec::new();
    let mut pockets: Vec<HashSet<(usize, usize)>> = vec![HashSet::new(); protein.ligands.len()];

    for (li, ligand) in protein.ligands.iter().enumerate() {
        for (ci, chain) in protein.chains.iter().enumerate() {
            for (ri, residue) in chain.residues.iter().enumerate() {
                let mut min_d_sq = f64::MAX;
                let mut found = false;

                for latom in &ligand.atoms {
                    if latom.element.trim() == "H" {
                        continue;
                    }
                    for ratom in &residue.atoms {
                        if ratom.element.trim() == "H" {
                            continue;
                        }
                        let d = dist_sq(latom.x, latom.y, latom.z, ratom.x, ratom.y, ratom.z);
                        if d < min_d_sq {
                            min_d_sq = d;
                        }
                        if d <= cutoff_sq {
                            found = true;
                        }
                    }
                }

                if found {
                    contacts.push(LigandContact {
                        ligand_idx: li,
                        chain_idx: ci,
                        residue_idx: ri,
                        min_distance: min_d_sq.sqrt(),
                    });
                    pockets[li].insert((ci, ri));
                }
            }
        }
    }

    contacts.sort_by(|a, b| {
        a.min_distance
            .partial_cmp(&b.min_distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    BindingPocketAnalysis { contacts, pockets }
}

impl InterfaceAnalysis {
    /// Return counts of each interaction type: (hbonds, salt_bridges, hydrophobic, other).
    pub fn interaction_counts(&self) -> (usize, usize, usize, usize) {
        let mut hbonds = 0;
        let mut salt_bridges = 0;
        let mut hydrophobic = 0;
        let mut other = 0;
        for interaction in &self.interactions {
            match interaction.interaction_type {
                InteractionType::HydrogenBond => hbonds += 1,
                InteractionType::SaltBridge => salt_bridges += 1,
                InteractionType::HydrophobicContact => hydrophobic += 1,
                InteractionType::Other => other += 1,
            }
        }
        (hbonds, salt_bridges, hydrophobic, other)
    }

    /// Convert interface residues to (chain_id, seq_num) pairs using the protein.
    pub fn interface_residues_by_id_with_protein(
        &self,
        protein: &Protein,
    ) -> HashSet<(String, i32)> {
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
            .map(|(idx, chain)| format!("Chain {}: {}", chain.id, self.chain_interface_counts[idx]))
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

        // Ligand binding pocket summary.
        if let Some(ref bp) = self.binding_pockets {
            if !bp.contacts.is_empty() {
                lines.push(String::new());
                lines.push("Ligand Contacts".to_string());

                for (li, ligand) in protein.ligands.iter().enumerate() {
                    let pocket_size = bp.pockets.get(li).map(|p| p.len()).unwrap_or(0);
                    if pocket_size == 0 {
                        continue;
                    }

                    let min_dist = bp
                        .contacts
                        .iter()
                        .filter(|c| c.ligand_idx == li)
                        .map(|c| c.min_distance)
                        .fold(f64::MAX, f64::min);

                    let type_label = match ligand.ligand_type {
                        crate::model::protein::LigandType::Ion => "Ion",
                        crate::model::protein::LigandType::Ligand => "Ligand",
                    };

                    lines.push(format!(
                        "{} {} ({}:{}): {} res, min {:.1}\u{00C5}",
                        type_label,
                        ligand.name,
                        ligand.chain_id,
                        ligand.seq_num,
                        pocket_size,
                        min_dist,
                    ));

                    // Top 3 closest residues
                    let mut pocket_contacts: Vec<_> =
                        bp.contacts.iter().filter(|c| c.ligand_idx == li).collect();
                    pocket_contacts.sort_by(|a, b| {
                        a.min_distance
                            .partial_cmp(&b.min_distance)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    let top_n = 3.min(pocket_contacts.len());
                    let labels: Vec<String> = pocket_contacts[..top_n]
                        .iter()
                        .map(|c| {
                            let res = &protein.chains[c.chain_idx].residues[c.residue_idx];
                            format!(
                                "{}:{}{} ({:.1}\u{00C5})",
                                protein.chains[c.chain_idx].id,
                                res.name,
                                res.seq_num,
                                c.min_distance
                            )
                        })
                        .collect();
                    if !labels.is_empty() {
                        lines.push(format!("  {}", labels.join(", ")));
                    }
                }
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
                is_hetero: false,
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
                        make_residue("ASP", 1, 3.0, 0.0, 0.0),  // within 4.5 of A:ALA1
                        make_residue("LEU", 2, 20.0, 0.0, 0.0), // far away
                    ],
                    molecule_type: MoleculeType::Protein,
                },
            ],
            ligands: Vec::new(),
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
                            is_hetero: false,
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
                            is_hetero: false,
                        }],
                        secondary_structure: SecondaryStructure::Coil,
                    }],
                    molecule_type: MoleculeType::Protein,
                },
            ],
            ligands: Vec::new(),
        };

        let analysis = analyze_interface(&protein, 4.5);
        assert!(
            analysis.contacts.is_empty(),
            "hydrogen-only atoms must not produce contacts"
        );
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
            ligands: Vec::new(),
        };
        let analysis = analyze_interface(&protein, 4.5);

        assert!(analysis.contacts.is_empty());
        assert_eq!(analysis.total_interface_residues, 0);
        let lines = analysis.summary(&protein);
        assert_eq!(lines.len(), 1);
        assert!(lines[0].contains("no inter-chain contacts"));
    }

    #[test]
    fn test_binding_pocket_detection() {
        use crate::model::protein::{Ligand, LigandType};

        let mut protein = two_chain_protein();
        // Add a ligand near chain A residue 0 (at 0,0,0)
        protein.ligands.push(Ligand {
            name: "HEM".to_string(),
            chain_id: "A".to_string(),
            seq_num: 100,
            atoms: vec![Atom {
                name: "FE".to_string(),
                element: "Fe".to_string(),
                x: 2.0,
                y: 0.0,
                z: 0.0,
                b_factor: 0.0,
                is_backbone: false,
                is_hetero: true,
            }],
            ligand_type: LigandType::Ligand,
        });

        let bp = analyze_binding_pockets(&protein, 4.5);
        assert_eq!(bp.pockets.len(), 1);
        // Should detect chain A residue 0 (at 0,0,0) as in pocket (dist = 2.0)
        assert!(bp.pockets[0].contains(&(0, 0)));
        // Should also detect chain B residue 0 (at 3,0,0) as in pocket (dist = 1.0)
        assert!(bp.pockets[0].contains(&(1, 0)));
        // Chain A residue 1 (at 10,0,0) should NOT be in pocket
        assert!(!bp.pockets[0].contains(&(0, 1)));
    }

    #[test]
    fn test_binding_pocket_empty_ligands() {
        let protein = two_chain_protein(); // has no ligands
        let bp = analyze_binding_pockets(&protein, 4.5);
        assert!(bp.contacts.is_empty());
        assert!(bp.pockets.is_empty());
    }

    #[test]
    fn test_binding_pocket_ion() {
        use crate::model::protein::{Ligand, LigandType};

        let mut protein = two_chain_protein();
        protein.ligands.push(Ligand {
            name: "ZN".to_string(),
            chain_id: "A".to_string(),
            seq_num: 200,
            atoms: vec![Atom {
                name: "ZN".to_string(),
                element: "Zn".to_string(),
                x: 1.0,
                y: 0.0,
                z: 0.0,
                b_factor: 0.0,
                is_backbone: false,
                is_hetero: true,
            }],
            ligand_type: LigandType::Ion,
        });

        let bp = analyze_binding_pockets(&protein, 4.5);
        assert!(!bp.contacts.is_empty());
        assert!(bp.pockets[0].contains(&(0, 0))); // chain A res 0 at origin, dist 1.0
    }

    // ---------------------------------------------------------------
    // Helper for classify_interactions tests
    // ---------------------------------------------------------------

    /// Build an `Atom` with the given name, element, and coordinates.
    fn make_atom(name: &str, element: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom {
            name: name.to_string(),
            element: element.to_string(),
            x,
            y,
            z,
            b_factor: 0.0,
            is_backbone: name == "CA" || name == "C" || name == "N" || name == "O",
            is_hetero: false,
        }
    }

    /// Build a two-chain protein where each chain has a single residue with
    /// the given atoms. Returns the protein and a pre-built contact vector
    /// pointing at residue 0 on each chain, ready for `classify_interactions`.
    fn interaction_protein(
        res_a_name: &str,
        atoms_a: Vec<Atom>,
        res_b_name: &str,
        atoms_b: Vec<Atom>,
    ) -> (Protein, Vec<Contact>) {
        let protein = Protein {
            name: "interaction_test".to_string(),
            chains: vec![
                Chain {
                    id: "A".to_string(),
                    residues: vec![Residue {
                        name: res_a_name.to_string(),
                        seq_num: 1,
                        atoms: atoms_a,
                        secondary_structure: SecondaryStructure::Coil,
                    }],
                    molecule_type: MoleculeType::Protein,
                },
                Chain {
                    id: "B".to_string(),
                    residues: vec![Residue {
                        name: res_b_name.to_string(),
                        seq_num: 1,
                        atoms: atoms_b,
                        secondary_structure: SecondaryStructure::Coil,
                    }],
                    molecule_type: MoleculeType::Protein,
                },
            ],
            ligands: Vec::new(),
        };
        let contacts = vec![Contact {
            chain_a: 0,
            residue_a: 0,
            chain_b: 1,
            residue_b: 0,
            min_distance: 0.0, // placeholder; classify_interactions recomputes
        }];
        (protein, contacts)
    }

    // ---------------------------------------------------------------
    // classify_interactions tests
    // ---------------------------------------------------------------

    #[test]
    fn test_classify_salt_bridge() {
        // ARG NH2 at origin, ASP OD1 at 3.5 A along x => distance 3.5 <= 4.0
        let (protein, contacts) = interaction_protein(
            "ARG",
            vec![make_atom("NH2", "N", 0.0, 0.0, 0.0)],
            "ASP",
            vec![make_atom("OD1", "O", 3.5, 0.0, 0.0)],
        );
        let interactions = classify_interactions(&protein, &contacts);
        assert_eq!(interactions.len(), 1);
        assert_eq!(
            interactions[0].interaction_type,
            InteractionType::SaltBridge
        );
        assert!((interactions[0].distance - 3.5).abs() < 1e-9);
    }

    #[test]
    fn test_classify_hydrogen_bond() {
        // SER N at origin, THR O at 3.0 A => H-bond (N/O donor-acceptor, <= 3.5)
        let (protein, contacts) = interaction_protein(
            "SER",
            vec![make_atom("N", "N", 0.0, 0.0, 0.0)],
            "THR",
            vec![make_atom("O", "O", 3.0, 0.0, 0.0)],
        );
        let interactions = classify_interactions(&protein, &contacts);
        assert_eq!(interactions.len(), 1);
        assert_eq!(
            interactions[0].interaction_type,
            InteractionType::HydrogenBond
        );
    }

    #[test]
    fn test_classify_hydrophobic_contact() {
        // ALA C at origin, VAL C at 4.0 A => hydrophobic (C-C on hydrophobic residues, <= 4.5)
        let (protein, contacts) = interaction_protein(
            "ALA",
            vec![make_atom("CB", "C", 0.0, 0.0, 0.0)],
            "VAL",
            vec![make_atom("CB", "C", 4.0, 0.0, 0.0)],
        );
        let interactions = classify_interactions(&protein, &contacts);
        assert_eq!(interactions.len(), 1);
        assert_eq!(
            interactions[0].interaction_type,
            InteractionType::HydrophobicContact
        );
    }

    #[test]
    fn test_classify_other_fallback() {
        // ASP C at origin, SER C at 3.0 A => C-C but neither residue is hydrophobic
        // and elements are C (not N/O), so falls through to Other.
        let (protein, contacts) = interaction_protein(
            "ASP",
            vec![make_atom("CB", "C", 0.0, 0.0, 0.0)],
            "SER",
            vec![make_atom("CB", "C", 3.0, 0.0, 0.0)],
        );
        let interactions = classify_interactions(&protein, &contacts);
        assert_eq!(interactions.len(), 1);
        assert_eq!(interactions[0].interaction_type, InteractionType::Other);
    }

    #[test]
    fn test_classify_salt_bridge_priority_over_hbond() {
        // ARG NH2 (N) and ASP OD1 (O) at 3.0 A satisfies BOTH salt bridge (<= 4.0)
        // AND H-bond (<= 3.5 with N/O). Salt bridge should win because it is
        // checked first in the classification cascade.
        let (protein, contacts) = interaction_protein(
            "ARG",
            vec![make_atom("NH2", "N", 0.0, 0.0, 0.0)],
            "ASP",
            vec![make_atom("OD1", "O", 3.0, 0.0, 0.0)],
        );
        let interactions = classify_interactions(&protein, &contacts);
        assert_eq!(interactions.len(), 1);
        assert_eq!(
            interactions[0].interaction_type,
            InteractionType::SaltBridge
        );
    }

    // ---------------------------------------------------------------
    // Helper predicate tests
    // ---------------------------------------------------------------

    #[test]
    fn test_is_charged_positive() {
        assert!(is_charged_positive("ARG", "NH1"));
        assert!(is_charged_positive("ARG", "NH2"));
        assert!(is_charged_positive("ARG", "NE"));
        assert!(is_charged_positive("LYS", "NZ"));
        // Non-charged atoms on charged residues should be false.
        assert!(!is_charged_positive("LYS", "CA"));
        assert!(!is_charged_positive("ARG", "CA"));
        // Unrelated residues.
        assert!(!is_charged_positive("ALA", "N"));
        assert!(!is_charged_positive("SER", "OG"));
    }

    #[test]
    fn test_is_charged_negative() {
        assert!(is_charged_negative("ASP", "OD1"));
        assert!(is_charged_negative("ASP", "OD2"));
        assert!(is_charged_negative("GLU", "OE1"));
        assert!(is_charged_negative("GLU", "OE2"));
        // Non-charged atoms on charged residues.
        assert!(!is_charged_negative("ASP", "CA"));
        assert!(!is_charged_negative("GLU", "CA"));
        // Unrelated residues.
        assert!(!is_charged_negative("ALA", "O"));
    }

    #[test]
    fn test_is_hydrophobic_residue() {
        // All hydrophobic residues.
        for name in &["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "MET", "PRO"] {
            assert!(
                is_hydrophobic_residue(name),
                "{} should be hydrophobic",
                name
            );
        }
        // Polar / charged residues.
        for name in &[
            "ASP", "GLU", "LYS", "ARG", "SER", "THR", "ASN", "GLN", "HIS", "GLY", "CYS", "TYR",
        ] {
            assert!(
                !is_hydrophobic_residue(name),
                "{} should NOT be hydrophobic",
                name
            );
        }
    }

    // ---------------------------------------------------------------
    // interaction_counts() test
    // ---------------------------------------------------------------

    #[test]
    fn test_interaction_counts() {
        let analysis = InterfaceAnalysis {
            contacts: Vec::new(),
            interface_residues: HashSet::new(),
            chain_interface_counts: Vec::new(),
            total_interface_residues: 0,
            binding_pockets: None,
            interactions: vec![
                Interaction {
                    interaction_type: InteractionType::HydrogenBond,
                    atom_a: [0.0; 3],
                    atom_b: [1.0, 0.0, 0.0],
                    distance: 1.0,
                },
                Interaction {
                    interaction_type: InteractionType::HydrogenBond,
                    atom_a: [0.0; 3],
                    atom_b: [2.0, 0.0, 0.0],
                    distance: 2.0,
                },
                Interaction {
                    interaction_type: InteractionType::SaltBridge,
                    atom_a: [0.0; 3],
                    atom_b: [3.0, 0.0, 0.0],
                    distance: 3.0,
                },
                Interaction {
                    interaction_type: InteractionType::HydrophobicContact,
                    atom_a: [0.0; 3],
                    atom_b: [4.0, 0.0, 0.0],
                    distance: 4.0,
                },
                Interaction {
                    interaction_type: InteractionType::HydrophobicContact,
                    atom_a: [0.0; 3],
                    atom_b: [4.0, 0.0, 0.0],
                    distance: 4.0,
                },
                Interaction {
                    interaction_type: InteractionType::HydrophobicContact,
                    atom_a: [0.0; 3],
                    atom_b: [4.0, 0.0, 0.0],
                    distance: 4.0,
                },
                Interaction {
                    interaction_type: InteractionType::Other,
                    atom_a: [0.0; 3],
                    atom_b: [5.0, 0.0, 0.0],
                    distance: 5.0,
                },
            ],
        };

        let (hbonds, salt_bridges, hydrophobic, other) = analysis.interaction_counts();
        assert_eq!(hbonds, 2);
        assert_eq!(salt_bridges, 1);
        assert_eq!(hydrophobic, 3);
        assert_eq!(other, 1);
    }
}

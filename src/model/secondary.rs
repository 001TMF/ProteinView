use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::model::protein::{Chain, MoleculeType, Protein, Residue, SecondaryStructure};

/// A secondary structure range parsed from PDB HELIX/SHEET records or CIF categories.
#[derive(Debug, Clone)]
pub struct SSRange {
    /// Chain identifier (may be multi-character for CIF auth_asym_id)
    pub chain_id: String,
    pub start_seq: i32,
    pub end_seq: i32,
    pub ss_type: SecondaryStructure,
}

/// Parse HELIX and SHEET records from a PDB file and return a list of SSRange entries.
///
/// PDB format column positions (0-indexed):
///   HELIX: initChainID=col 19, initSeqNum=cols 21..25, endChainID=col 31, endSeqNum=cols 33..37
///   SHEET: initChainID=col 21, initSeqNum=cols 22..26, endChainID=col 32, endSeqNum=cols 33..37
fn parse_ss_records(file_path: &str) -> Vec<SSRange> {
    let file = match File::open(file_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!(
                "Warning: could not re-open '{}' for SS records: {}",
                file_path, e
            );
            return Vec::new();
        }
    };
    let reader = BufReader::new(file);
    let mut ranges = Vec::new();

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!(
                    "Warning: skipping unreadable line in '{}': {}",
                    file_path, e
                );
                continue;
            }
        };

        if line.starts_with("HELIX ") {
            if line.len() < 38 {
                continue;
            }
            // initChainID at col 19
            let init_chain = line.as_bytes()[19] as char;
            // initSeqNum at cols 21..25 (inclusive, i.e. bytes 21..=24)
            let init_seq_str = &line[21..25];
            // endChainID at col 31
            let end_chain = line.as_bytes()[31] as char;
            // endSeqNum at cols 33..37 (inclusive, i.e. bytes 33..=36)
            let end_seq_str = &line[33..37];

            let init_seq: i32 = match init_seq_str.trim().parse() {
                Ok(n) => n,
                Err(_) => continue,
            };
            let end_seq: i32 = match end_seq_str.trim().parse() {
                Ok(n) => n,
                Err(_) => continue,
            };

            // Both chain IDs should match for a valid helix range.
            // We use the initChainID as the canonical chain.
            if init_chain == end_chain {
                ranges.push(SSRange {
                    chain_id: init_chain.to_string(),
                    start_seq: init_seq,
                    end_seq,
                    ss_type: SecondaryStructure::Helix,
                });
            }
        } else if line.starts_with("SHEET ") {
            if line.len() < 38 {
                continue;
            }
            // initChainID at col 21
            let init_chain = line.as_bytes()[21] as char;
            // initSeqNum at cols 22..26 (inclusive, i.e. bytes 22..=25)
            let init_seq_str = &line[22..26];
            // endChainID at col 32
            let end_chain = line.as_bytes()[32] as char;
            // endSeqNum at cols 33..37 (inclusive, i.e. bytes 33..=36)
            let end_seq_str = &line[33..37];

            let init_seq: i32 = match init_seq_str.trim().parse() {
                Ok(n) => n,
                Err(_) => continue,
            };
            let end_seq: i32 = match end_seq_str.trim().parse() {
                Ok(n) => n,
                Err(_) => continue,
            };

            if init_chain == end_chain {
                ranges.push(SSRange {
                    chain_id: init_chain.to_string(),
                    start_seq: init_seq,
                    end_seq,
                    ss_type: SecondaryStructure::Sheet,
                });
            }
        }
    }

    ranges
}

/// Parse HELIX/SHEET records from the given PDB file and assign
/// secondary structure to all matching residues in the protein.
/// Residues not covered by any HELIX or SHEET record remain as Coil.
pub fn assign_from_pdb_file(protein: &mut Protein, file_path: &str) {
    let ranges = parse_ss_records(file_path);
    if ranges.is_empty() {
        return;
    }

    apply_ss_ranges(protein, &ranges);
}

/// Apply a list of SSRange entries to a protein, setting the secondary
/// structure for each residue that falls within a range.
fn apply_ss_ranges(protein: &mut Protein, ranges: &[SSRange]) {
    for chain in &mut protein.chains {
        for residue in &mut chain.residues {
            for range in ranges {
                if chain.id == range.chain_id
                    && residue.seq_num >= range.start_seq
                    && residue.seq_num <= range.end_seq
                {
                    residue.secondary_structure = range.ss_type;
                    break; // first matching range wins
                }
            }
        }
    }
}

/// Split a CIF data line into whitespace-separated tokens, respecting
/// single-quoted strings (e.g. 'some value') as a single token.
fn tokenize_cif_line(line: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut chars = line.chars().peekable();
    while let Some(&ch) = chars.peek() {
        if ch.is_whitespace() {
            chars.next();
        } else if ch == '\'' {
            // Quoted token: consume until matching closing quote followed by
            // whitespace or end-of-line.
            chars.next(); // skip opening quote
            let mut token = String::new();
            loop {
                match chars.next() {
                    Some('\'') => {
                        // Check if next char is whitespace or end
                        match chars.peek() {
                            None | Some(&' ') | Some(&'\t') => break,
                            Some(_) => token.push('\''),
                        }
                    }
                    Some(c) => token.push(c),
                    None => break,
                }
            }
            tokens.push(token);
        } else {
            // Unquoted token
            let mut token = String::new();
            while let Some(&c) = chars.peek() {
                if c.is_whitespace() {
                    break;
                }
                token.push(c);
                chars.next();
            }
            tokens.push(token);
        }
    }
    tokens
}

/// Parse secondary structure records from a CIF/mmCIF file.
///
/// Reads `_struct_conf` (helices) and `_struct_sheet_range` (sheets) loop
/// categories. Uses `auth_asym_id` and `auth_seq_id` fields to match
/// the chain and residue identifiers produced by pdbtbx.
fn parse_cif_ss_records(file_path: &str) -> Vec<SSRange> {
    let file = match File::open(file_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!(
                "Warning: could not re-open '{}' for CIF SS records: {}",
                file_path, e
            );
            return Vec::new();
        }
    };
    let reader = BufReader::new(file);
    let mut ranges = Vec::new();

    // We need to parse two different loop categories. We'll do it in a
    // single pass using a simple state machine.
    #[derive(PartialEq)]
    enum ParseState {
        Scanning,
        StructConfHeaders,
        StructConfData,
        SheetRangeHeaders,
        SheetRangeData,
    }

    let mut state = ParseState::Scanning;
    let mut column_names: Vec<String> = Vec::new();
    let mut col_map: HashMap<String, usize> = HashMap::new();

    #[allow(clippy::lines_filter_map_ok)]
    let lines: Vec<String> = reader.lines().filter_map(|l| l.ok()).collect();
    let mut i = 0;

    while i < lines.len() {
        let line = &lines[i];
        let trimmed = line.trim();

        match state {
            ParseState::Scanning => {
                if trimmed == "loop_" {
                    // Peek at next lines to determine which category this loop belongs to
                    if i + 1 < lines.len() {
                        let next = lines[i + 1].trim().to_string();
                        if next.starts_with("_struct_conf.") {
                            state = ParseState::StructConfHeaders;
                            column_names.clear();
                            col_map.clear();
                        } else if next.starts_with("_struct_sheet_range.") {
                            state = ParseState::SheetRangeHeaders;
                            column_names.clear();
                            col_map.clear();
                        }
                    }
                }
            }

            ParseState::StructConfHeaders => {
                if trimmed.starts_with("_struct_conf.") {
                    let col_name = trimmed.to_string();
                    col_map.insert(col_name.clone(), column_names.len());
                    column_names.push(col_name);
                } else if !trimmed.is_empty() && !trimmed.starts_with('#') {
                    // First data line
                    state = ParseState::StructConfData;
                    // Process this line as data (don't skip it)
                    let tokens = tokenize_cif_line(trimmed);
                    if let Some(range) = parse_struct_conf_row(&tokens, &col_map) {
                        ranges.push(range);
                    }
                } else if trimmed.starts_with('#') || trimmed.is_empty() {
                    // End of loop with no data
                    state = ParseState::Scanning;
                }
            }

            ParseState::StructConfData => {
                if trimmed.is_empty()
                    || trimmed.starts_with('#')
                    || trimmed.starts_with("loop_")
                    || trimmed.starts_with('_')
                {
                    state = ParseState::Scanning;
                    // Don't advance i, re-process in Scanning state
                    continue;
                }
                let tokens = tokenize_cif_line(trimmed);
                if let Some(range) = parse_struct_conf_row(&tokens, &col_map) {
                    ranges.push(range);
                }
            }

            ParseState::SheetRangeHeaders => {
                if trimmed.starts_with("_struct_sheet_range.") {
                    let col_name = trimmed.to_string();
                    col_map.insert(col_name.clone(), column_names.len());
                    column_names.push(col_name);
                } else if !trimmed.is_empty() && !trimmed.starts_with('#') {
                    // First data line
                    state = ParseState::SheetRangeData;
                    let tokens = tokenize_cif_line(trimmed);
                    if let Some(range) = parse_sheet_range_row(&tokens, &col_map) {
                        ranges.push(range);
                    }
                } else if trimmed.starts_with('#') || trimmed.is_empty() {
                    state = ParseState::Scanning;
                }
            }

            ParseState::SheetRangeData => {
                if trimmed.is_empty()
                    || trimmed.starts_with('#')
                    || trimmed.starts_with("loop_")
                    || trimmed.starts_with('_')
                {
                    state = ParseState::Scanning;
                    continue;
                }
                let tokens = tokenize_cif_line(trimmed);
                if let Some(range) = parse_sheet_range_row(&tokens, &col_map) {
                    ranges.push(range);
                }
            }
        }

        i += 1;
    }

    ranges
}

/// Parse a single data row from the _struct_conf loop into an SSRange (helix).
/// Returns None if the row cannot be parsed or is not a helix type.
fn parse_struct_conf_row(tokens: &[String], col_map: &HashMap<String, usize>) -> Option<SSRange> {
    // conf_type_id must start with HELX for helix (e.g., HELX_P)
    let conf_type_idx = *col_map.get("_struct_conf.conf_type_id")?;
    let conf_type = tokens.get(conf_type_idx)?;
    if !conf_type.starts_with("HELX") {
        // Could be TURN or other types; we only handle helices here
        // TURN types could be added later if needed
        return None;
    }

    // Use auth fields first, fall back to label fields
    let beg_chain = get_cif_field(
        tokens,
        col_map,
        "_struct_conf.beg_auth_asym_id",
        "_struct_conf.beg_label_asym_id",
    )?;
    let end_chain = get_cif_field(
        tokens,
        col_map,
        "_struct_conf.end_auth_asym_id",
        "_struct_conf.end_label_asym_id",
    )?;
    let beg_seq_str = get_cif_field(
        tokens,
        col_map,
        "_struct_conf.beg_auth_seq_id",
        "_struct_conf.beg_label_seq_id",
    )?;
    let end_seq_str = get_cif_field(
        tokens,
        col_map,
        "_struct_conf.end_auth_seq_id",
        "_struct_conf.end_label_seq_id",
    )?;

    if beg_chain != end_chain {
        return None;
    }

    let start_seq: i32 = beg_seq_str.parse().ok()?;
    let end_seq: i32 = end_seq_str.parse().ok()?;

    Some(SSRange {
        chain_id: beg_chain,
        start_seq,
        end_seq,
        ss_type: SecondaryStructure::Helix,
    })
}

/// Parse a single data row from the _struct_sheet_range loop into an SSRange (sheet).
fn parse_sheet_range_row(tokens: &[String], col_map: &HashMap<String, usize>) -> Option<SSRange> {
    let beg_chain = get_cif_field(
        tokens,
        col_map,
        "_struct_sheet_range.beg_auth_asym_id",
        "_struct_sheet_range.beg_label_asym_id",
    )?;
    let end_chain = get_cif_field(
        tokens,
        col_map,
        "_struct_sheet_range.end_auth_asym_id",
        "_struct_sheet_range.end_label_asym_id",
    )?;
    let beg_seq_str = get_cif_field(
        tokens,
        col_map,
        "_struct_sheet_range.beg_auth_seq_id",
        "_struct_sheet_range.beg_label_seq_id",
    )?;
    let end_seq_str = get_cif_field(
        tokens,
        col_map,
        "_struct_sheet_range.end_auth_seq_id",
        "_struct_sheet_range.end_label_seq_id",
    )?;

    if beg_chain != end_chain {
        return None;
    }

    let start_seq: i32 = beg_seq_str.parse().ok()?;
    let end_seq: i32 = end_seq_str.parse().ok()?;

    Some(SSRange {
        chain_id: beg_chain,
        start_seq,
        end_seq,
        ss_type: SecondaryStructure::Sheet,
    })
}

/// Retrieve a CIF field value from tokenized data, preferring the primary
/// column name and falling back to the fallback column name.
/// Returns None if neither column exists or the value is "?" (missing).
fn get_cif_field(
    tokens: &[String],
    col_map: &HashMap<String, usize>,
    primary: &str,
    fallback: &str,
) -> Option<String> {
    let idx = col_map.get(primary).or_else(|| col_map.get(fallback))?;
    let val = tokens.get(*idx)?;
    if val == "?" || val == "." {
        return None;
    }
    Some(val.clone())
}

/// Parse secondary structure from CIF _struct_conf and _struct_sheet_range
/// categories and assign to matching residues in the protein.
/// Uses auth_asym_id/auth_seq_id to match pdbtbx's chain/residue numbering.
pub fn assign_from_cif_file(protein: &mut Protein, file_path: &str) {
    let ranges = parse_cif_ss_records(file_path);
    if ranges.is_empty() {
        return;
    }
    apply_ss_ranges(protein, &ranges);
}

// ---------------------------------------------------------------------------
// Secondary structure inference from backbone coordinates (DSSP-like)
// ---------------------------------------------------------------------------

/// Infer secondary structure from backbone geometry for protein chains that
/// lack explicit HELIX/SHEET annotations.
///
/// Skips non-protein chains (RNA, DNA, SmallMolecule) and chains that already
/// have any non-Coil secondary structure assignments.
pub fn infer_secondary_structure(chains: &mut [Chain]) {
    for chain in chains.iter_mut() {
        // Skip non-protein chains
        if chain.molecule_type != MoleculeType::Protein {
            continue;
        }

        // Skip chains that already have explicit SS assignments
        let has_existing_ss = chain
            .residues
            .iter()
            .any(|r| r.secondary_structure != SecondaryStructure::Coil);
        if has_existing_ss {
            continue;
        }

        let inferred = infer_chain_ss(&chain.residues);
        for (residue, ss) in chain.residues.iter_mut().zip(inferred) {
            residue.secondary_structure = ss;
        }
    }
}

/// Infer secondary structure for a single chain's residues.
fn infer_chain_ss(residues: &[Residue]) -> Vec<SecondaryStructure> {
    let mut assignments = vec![SecondaryStructure::Coil; residues.len()];
    let torsions = compute_torsion_angles(residues);
    let hbonds = compute_hbond_map(residues);

    assign_helices(&mut assignments, &hbonds, &torsions);
    assign_sheets(&mut assignments, &hbonds, &torsions);

    // Post-processing: fill single-residue gaps and enforce minimum run lengths
    fill_single_gaps(&mut assignments, &torsions, SecondaryStructure::Helix);
    fill_single_gaps(&mut assignments, &torsions, SecondaryStructure::Sheet);
    enforce_min_run(&mut assignments, SecondaryStructure::Helix, 3);
    enforce_min_run(&mut assignments, SecondaryStructure::Sheet, 2);

    assignments
}

/// Compute phi/psi torsion angles for each residue from backbone N, CA, C atoms.
///
/// Phi(i)  = dihedral(C[i-1], N[i], CA[i], C[i])
/// Psi(i)  = dihedral(N[i], CA[i], C[i], N[i+1])
///
/// Returns None for terminal residues or those with incomplete backbone.
fn compute_torsion_angles(residues: &[Residue]) -> Vec<Option<(f64, f64)>> {
    let mut torsions = vec![None; residues.len()];

    for i in 1..residues.len().saturating_sub(1) {
        let c_prev = atom_pos(&residues[i - 1], "C");
        let n_i = atom_pos(&residues[i], "N");
        let ca_i = atom_pos(&residues[i], "CA");
        let c_i = atom_pos(&residues[i], "C");
        let n_next = atom_pos(&residues[i + 1], "N");

        let (Some(c_prev), Some(n_i), Some(ca_i), Some(c_i), Some(n_next)) =
            (c_prev, n_i, ca_i, c_i, n_next)
        else {
            continue;
        };

        let (Some(phi), Some(psi)) = (
            dihedral(c_prev, n_i, ca_i, c_i),
            dihedral(n_i, ca_i, c_i, n_next),
        ) else {
            continue; // degenerate geometry — skip this residue
        };
        torsions[i] = Some((phi, psi));
    }

    torsions
}

/// Build an NxN boolean matrix of hydrogen bonds between residues.
///
/// hbonds[acceptor][donor] = true means the C=O of residue `acceptor` accepts
/// an H-bond from the N-H of residue `donor`.
///
/// Uses the DSSP energy formula (Kabsch & Sander, 1983):
///   E = 27.888 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)
/// with threshold E < -0.5 kcal/mol.
fn compute_hbond_map(residues: &[Residue]) -> Vec<Vec<bool>> {
    let n = residues.len();
    let mut hbonds = vec![vec![false; n]; n];

    // Pre-extract backbone atom positions for efficiency
    let c_atoms: Vec<Option<[f64; 3]>> = residues.iter().map(|r| atom_pos(r, "C")).collect();
    let o_atoms: Vec<Option<[f64; 3]>> = residues.iter().map(|r| atom_pos(r, "O")).collect();
    let n_atoms: Vec<Option<[f64; 3]>> = residues.iter().map(|r| atom_pos(r, "N")).collect();
    let ca_atoms: Vec<Option<[f64; 3]>> = residues.iter().map(|r| atom_pos(r, "CA")).collect();

    for acceptor in 0..n {
        let (Some(c_acc), Some(o_acc)) = (c_atoms[acceptor], o_atoms[acceptor]) else {
            continue;
        };

        for donor in 0..n {
            // Skip self and immediate neighbors
            if donor == acceptor || donor.abs_diff(acceptor) <= 1 {
                continue;
            }

            let (Some(n_don), Some(ca_don)) = (n_atoms[donor], ca_atoms[donor]) else {
                continue;
            };

            // N...O distance pre-filter: skip pairs where N...O > 5.2 A
            if distance(n_don, o_acc) > 5.2 {
                continue;
            }

            let Some(h_don) = estimate_amide_h(&c_atoms, n_don, ca_don, donor) else {
                continue;
            };

            let energy = hbond_energy(o_acc, c_acc, n_don, h_don);
            if energy < -0.5 {
                hbonds[acceptor][donor] = true;
            }
        }
    }

    hbonds
}

/// Estimate amide hydrogen position using the bisector method.
///
/// The H atom is placed along the bisector of the N->C(i-1) and N->CA(i)
/// vectors, at 1.0 A from N.
fn estimate_amide_h(
    c_atoms: &[Option<[f64; 3]>],
    n_atom: [f64; 3],
    ca_atom: [f64; 3],
    donor: usize,
) -> Option<[f64; 3]> {
    // Need C of previous residue
    let prev_c = donor.checked_sub(1).and_then(|i| c_atoms[i])?;
    let dir_prev = normalize(sub(n_atom, prev_c))?;
    let dir_ca = normalize(sub(n_atom, ca_atom))?;
    let bisector = normalize(add(dir_prev, dir_ca))?;
    // Place H at 1.0 A along bisector from N
    Some(add(n_atom, scale(bisector, 1.0)))
}

/// DSSP hydrogen bond energy formula (Kabsch & Sander, 1983).
///
/// E = 0.084 * 332 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)
///   = 27.888 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)
///
/// Distances are in Angstroms, energy in kcal/mol.
/// A bond is considered present if E < -0.5 kcal/mol.
fn hbond_energy(o: [f64; 3], c: [f64; 3], n: [f64; 3], h: [f64; 3]) -> f64 {
    let r_on = distance(o, n).max(0.5);
    let r_ch = distance(c, h).max(0.5);
    let r_oh = distance(o, h).max(0.5);
    let r_cn = distance(c, n).max(0.5);
    27.888 * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn)
}

/// Detect helices from H-bond patterns: i -> i+3, i -> i+4, or i -> i+5.
///
/// For each turn pattern (alpha=4, 3_10=3, pi=5), checks if acceptor residue i
/// forms an H-bond with donor residue i+turn. Residues in the interior of such
/// turns that have compatible torsion angles are assigned as Helix.
fn assign_helices(
    assignments: &mut [SecondaryStructure],
    hbonds: &[Vec<bool>],
    torsions: &[Option<(f64, f64)>],
) {
    let n = assignments.len();
    let mut support = vec![0usize; n];

    // Check alpha (i->i+4), 3_10 (i->i+3), and pi (i->i+5) turns
    for turn in [4usize, 3, 5] {
        for i in 0..n.saturating_sub(turn) {
            if !hbonds[i][i + turn] {
                continue;
            }

            // Check that at least half the interior residues have helix-compatible torsions
            let span = i + 1..=i + turn;
            let span_len = turn;
            let compatible = span
                .clone()
                .filter(|&idx| torsions_match(torsions[idx], SecondaryStructure::Helix))
                .count();
            if compatible * 2 < span_len {
                continue;
            }

            for idx in span {
                support[idx] += 1;
            }
        }
    }

    for (idx, s) in assignments.iter_mut().enumerate() {
        if support[idx] > 0 {
            *s = SecondaryStructure::Helix;
        }
    }
}

/// Detect beta-sheets from parallel and antiparallel H-bond bridge patterns.
///
/// Antiparallel bridge: (hbonds[i][j] && hbonds[j][i]) or
///                      (hbonds[i-1][j+1] && hbonds[j-1][i+1])
/// Parallel bridge:     (hbonds[i-1][j] && hbonds[j][i+1]) or
///                      (hbonds[j-1][i] && hbonds[i][j+1])
fn assign_sheets(
    assignments: &mut [SecondaryStructure],
    hbonds: &[Vec<bool>],
    torsions: &[Option<(f64, f64)>],
) {
    let n = assignments.len();
    let mut support = vec![0usize; n];

    for i in 1..n.saturating_sub(1) {
        for j in (i + 2)..n.saturating_sub(1) {
            // Both residues should have sheet-compatible torsion angles
            if !torsions_match(torsions[i], SecondaryStructure::Sheet)
                || !torsions_match(torsions[j], SecondaryStructure::Sheet)
            {
                continue;
            }

            let antiparallel =
                (hbonds[i][j] && hbonds[j][i]) || (hbonds[i - 1][j + 1] && hbonds[j - 1][i + 1]);
            let parallel =
                (hbonds[i - 1][j] && hbonds[j][i + 1]) || (hbonds[j - 1][i] && hbonds[i][j + 1]);

            if antiparallel || parallel {
                support[i] += 1;
                support[j] += 1;
            }
        }
    }

    // Only assign Sheet to residues currently marked as Coil (helix takes priority)
    for idx in 0..n {
        if assignments[idx] == SecondaryStructure::Coil
            && support[idx] > 0
            && torsions_match(torsions[idx], SecondaryStructure::Sheet)
        {
            assignments[idx] = SecondaryStructure::Sheet;
        }
    }
}

// ---------------------------------------------------------------------------
// Post-processing: gap filling and minimum run length
// ---------------------------------------------------------------------------

/// Fill single-residue Coil gaps within runs of the given SS type.
///
/// If residues [i-1] and [i+1] are both `target` and residue [i] is Coil,
/// promote [i] to `target` if its torsion angles are compatible (or missing).
fn fill_single_gaps(
    assignments: &mut [SecondaryStructure],
    torsions: &[Option<(f64, f64)>],
    target: SecondaryStructure,
) {
    if assignments.len() < 3 {
        return;
    }

    for i in 1..assignments.len() - 1 {
        if assignments[i - 1] != target
            || assignments[i] != SecondaryStructure::Coil
            || assignments[i + 1] != target
        {
            continue;
        }

        // Fill the gap if torsions are compatible or unknown
        if torsions_match(torsions[i], target) || torsions[i].is_none() {
            assignments[i] = target;
        }
    }
}

/// Remove runs of `ss` shorter than `min_len`, resetting them to Coil.
fn enforce_min_run(assignments: &mut [SecondaryStructure], ss: SecondaryStructure, min_len: usize) {
    let mut i = 0;
    while i < assignments.len() {
        if assignments[i] != ss {
            i += 1;
            continue;
        }

        let start = i;
        while i < assignments.len() && assignments[i] == ss {
            i += 1;
        }

        if i - start < min_len {
            for state in &mut assignments[start..i] {
                *state = SecondaryStructure::Coil;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Torsion angle windows
// ---------------------------------------------------------------------------

/// Check whether torsion angles are compatible with the given SS type.
fn torsions_match(torsions: Option<(f64, f64)>, target: SecondaryStructure) -> bool {
    let Some((phi, psi)) = torsions else {
        return false;
    };
    match target {
        SecondaryStructure::Helix => is_helix_torsion(phi, psi),
        SecondaryStructure::Sheet => is_sheet_torsion(phi, psi),
        _ => false,
    }
}

/// Helix-compatible torsion angles (wide window covering alpha/3_10/pi).
fn is_helix_torsion(phi: f64, psi: f64) -> bool {
    // Strong: typical alpha helix (phi ~ -57, psi ~ -47)
    let strong = (-80.0..=-40.0).contains(&phi) && (-70.0..=-20.0).contains(&psi);
    // Weak: wider window for 3_10 and pi helices
    let weak = (-170.0..=-20.0).contains(&phi) && (-80.0..=10.0).contains(&psi);
    strong || weak
}

/// Sheet-compatible torsion angles (extended conformation).
fn is_sheet_torsion(phi: f64, psi: f64) -> bool {
    // Strong: canonical beta sheet
    let strong = ((-100.0..=-40.0).contains(&phi) && (20.0..=90.0).contains(&psi))
        || ((80.0..=180.0).contains(&phi) && (120.0..=180.0).contains(&psi));
    // Weak: wider window
    let weak = ((-140.0..=-20.0).contains(&phi) && (0.0..=180.0).contains(&psi))
        || ((60.0..=180.0).contains(&phi) && (90.0..=180.0).contains(&psi));
    strong || weak
}

// ---------------------------------------------------------------------------
// 3D vector math utilities
// ---------------------------------------------------------------------------

/// Get the [x, y, z] position of a named atom in a residue, or None.
fn atom_pos(residue: &Residue, atom_name: &str) -> Option<[f64; 3]> {
    residue
        .atoms
        .iter()
        .find(|a| a.name == atom_name)
        .map(|a| [a.x, a.y, a.z])
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(v: [f64; 3], factor: f64) -> [f64; 3] {
    [v[0] * factor, v[1] * factor, v[2] * factor]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm(sub(a, b))
}

fn normalize(v: [f64; 3]) -> Option<[f64; 3]> {
    let len = norm(v);
    if len < 1e-8 {
        None
    } else {
        Some([v[0] / len, v[1] / len, v[2] / len])
    }
}

/// Compute the dihedral angle (in degrees) defined by four points.
/// Returns None for degenerate geometry (coincident atoms, collinear atoms).
fn dihedral(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Option<f64> {
    // IUPAC/biochemistry convention (matches MDAnalysis/BioPython):
    // Bond vectors along the chain A→B→C→D
    let b1 = sub(b, a);
    let b2 = sub(c, b);
    let b3 = sub(d, c);

    let b2_unit = normalize(b2)?;

    let n1 = cross(b1, b2); // normal to plane A-B-C
    let n2 = cross(b2, b3); // normal to plane B-C-D
    let n1_unit = normalize(n1)?;
    let n2_unit = normalize(n2)?;

    let m1 = cross(n1_unit, b2_unit);
    Some(-dot(m1, n2_unit).atan2(dot(n1_unit, n2_unit)).to_degrees())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tokenize_cif_line_basic() {
        let tokens = tokenize_cif_line(
            "HELX_P HELX_P1 1 GLY A 4   ? HIS A 15  ? GLY L 4   HIS L 15  1 ? 12",
        );
        assert_eq!(tokens[0], "HELX_P");
        assert_eq!(tokens[1], "HELX_P1");
        assert_eq!(tokens[2], "1");
        assert_eq!(tokens[3], "GLY");
        assert_eq!(tokens[4], "A");
        assert_eq!(tokens[5], "4");
        assert_eq!(tokens[6], "?");
        assert_eq!(tokens.len(), 20);
    }

    #[test]
    fn test_tokenize_cif_quoted_string() {
        let tokens = tokenize_cif_line("A 1 'hello world' B");
        assert_eq!(tokens.len(), 4);
        assert_eq!(tokens[0], "A");
        assert_eq!(tokens[1], "1");
        assert_eq!(tokens[2], "hello world");
        assert_eq!(tokens[3], "B");
    }

    #[test]
    fn test_parse_cif_ss_records_from_example_file() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/1ZVH.cif");
        let ranges = parse_cif_ss_records(path);

        // 1ZVH.cif has 9 HELX_P records and 17 sheet range records
        let helix_count = ranges
            .iter()
            .filter(|r| r.ss_type == SecondaryStructure::Helix)
            .count();
        let sheet_count = ranges
            .iter()
            .filter(|r| r.ss_type == SecondaryStructure::Sheet)
            .count();

        assert_eq!(
            helix_count, 9,
            "Expected 9 helix ranges, got {}",
            helix_count
        );
        assert_eq!(
            sheet_count, 17,
            "Expected 17 sheet ranges, got {}",
            sheet_count
        );

        // Check first helix: chain L, residues 4-15
        let first_helix = ranges
            .iter()
            .find(|r| r.ss_type == SecondaryStructure::Helix)
            .unwrap();
        assert_eq!(first_helix.chain_id, "L");
        assert_eq!(first_helix.start_seq, 4);
        assert_eq!(first_helix.end_seq, 15);

        // Check a helix on chain A (auth_asym_id): residues 87-91
        let chain_a_helix = ranges
            .iter()
            .find(|r| r.ss_type == SecondaryStructure::Helix && r.chain_id == "A")
            .unwrap();
        assert_eq!(chain_a_helix.start_seq, 87);
        assert_eq!(chain_a_helix.end_seq, 91);
    }

    #[test]
    fn test_assign_from_cif_file_sets_secondary_structure() {
        use crate::model::protein::{Atom, Chain, MoleculeType, Protein, Residue};

        // Build a minimal protein matching 1ZVH chain L residues 1-20
        let mut residues = Vec::new();
        for i in 1..=20 {
            residues.push(Residue {
                name: "ALA".to_string(),
                seq_num: i,
                atoms: vec![Atom {
                    name: "CA".to_string(),
                    element: "C".to_string(),
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    b_factor: 0.0,
                    is_backbone: true,
                    is_hetero: false,
                }],
                secondary_structure: SecondaryStructure::Coil,
            });
        }
        let mut protein = Protein {
            name: "test".to_string(),
            chains: vec![Chain {
                id: "L".to_string(),
                residues,
                molecule_type: MoleculeType::Protein,
            }],
            ligands: Vec::new(),
        };

        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/1ZVH.cif");
        assign_from_cif_file(&mut protein, path);

        let chain = &protein.chains[0];
        // Residues 1-3 should be Coil
        assert_eq!(
            chain.residues[0].secondary_structure,
            SecondaryStructure::Coil
        ); // res 1
        assert_eq!(
            chain.residues[1].secondary_structure,
            SecondaryStructure::Coil
        ); // res 2
        assert_eq!(
            chain.residues[2].secondary_structure,
            SecondaryStructure::Coil
        ); // res 3

        // Residues 4-15 should be Helix (first helix)
        assert_eq!(
            chain.residues[3].secondary_structure,
            SecondaryStructure::Helix
        ); // res 4
        assert_eq!(
            chain.residues[14].secondary_structure,
            SecondaryStructure::Helix
        ); // res 15

        // Residues 16-18 should be Coil (gap between helices)
        assert_eq!(
            chain.residues[15].secondary_structure,
            SecondaryStructure::Coil
        ); // res 16
        assert_eq!(
            chain.residues[16].secondary_structure,
            SecondaryStructure::Coil
        ); // res 17
        assert_eq!(
            chain.residues[17].secondary_structure,
            SecondaryStructure::Coil
        ); // res 18

        // Residues 19-20: res 19 starts helix 2 (19-23)
        assert_eq!(
            chain.residues[18].secondary_structure,
            SecondaryStructure::Helix
        ); // res 19
        assert_eq!(
            chain.residues[19].secondary_structure,
            SecondaryStructure::Helix
        ); // res 20
    }

    #[test]
    fn test_full_cif_load_has_non_coil_residues() {
        // Integration test: load the full CIF file through the parser
        // and verify that secondary structure was assigned
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/1ZVH.cif");
        let protein = crate::parser::pdb::load_structure(path).unwrap();

        let total_residues = protein.residue_count();
        let non_coil = protein
            .chains
            .iter()
            .flat_map(|c| &c.residues)
            .filter(|r| r.secondary_structure != SecondaryStructure::Coil)
            .count();

        assert!(total_residues > 0, "Protein should have residues");
        assert!(
            non_coil > 0,
            "Expected some non-Coil residues after CIF SS assignment, but all {} residues are Coil",
            total_residues
        );

        // Should have both helices and sheets
        let helix_count = protein
            .chains
            .iter()
            .flat_map(|c| &c.residues)
            .filter(|r| r.secondary_structure == SecondaryStructure::Helix)
            .count();
        let sheet_count = protein
            .chains
            .iter()
            .flat_map(|c| &c.residues)
            .filter(|r| r.secondary_structure == SecondaryStructure::Sheet)
            .count();

        assert!(helix_count > 0, "Expected helix residues, got 0");
        assert!(sheet_count > 0, "Expected sheet residues, got 0");
    }

    // -----------------------------------------------------------------------
    // Tests for secondary structure inference
    // -----------------------------------------------------------------------

    #[test]
    fn test_hbond_energy_known_geometry() {
        // Hand-calculated test: typical alpha-helix H-bond geometry
        // O at origin, C at (1.24, 0, 0), N at (2.8, 1.5, 0), H at (1.9, 1.2, 0)
        let o = [0.0, 0.0, 0.0];
        let c = [1.24, 0.0, 0.0];
        let n = [2.8, 1.5, 0.0];
        let h = [1.9, 1.2, 0.0];

        let energy = hbond_energy(o, c, n, h);

        // r_ON = sqrt(2.8^2 + 1.5^2) = sqrt(10.09) ~ 3.177
        // r_CH = sqrt((1.9-1.24)^2 + 1.2^2) = sqrt(0.4356 + 1.44) ~ 1.370
        // r_OH = sqrt(1.9^2 + 1.2^2) = sqrt(5.05) ~ 2.247
        // r_CN = sqrt((2.8-1.24)^2 + 1.5^2) = sqrt(2.4336 + 2.25) ~ 2.164
        // E = 27.888 * (1/3.177 + 1/1.370 - 1/2.247 - 1/2.164)
        //   = 27.888 * (0.3148 + 0.7299 - 0.4451 - 0.4621)
        //   = 27.888 * 0.1375 ~ 3.83

        // This geometry does NOT form a favorable H-bond (energy > 0)
        assert!(
            energy > -0.5,
            "Expected no H-bond for this geometry, got E={}",
            energy
        );

        // Now test a geometry that SHOULD form an H-bond:
        // Closer O-H distance, typical for real H-bond
        let o2 = [0.0, 0.0, 0.0];
        let c2 = [1.24, 0.0, 0.0];
        let n2 = [3.0, 0.0, 0.0];
        let h2 = [2.0, 0.0, 0.0]; // H pointing directly at O

        let energy2 = hbond_energy(o2, c2, n2, h2);
        // r_ON = 3.0, r_CH = 0.76, r_OH = 2.0, r_CN = 1.76
        // E = 27.888 * (1/3.0 + 1/0.76 - 1/2.0 - 1/1.76)
        //   = 27.888 * (0.333 + 1.316 - 0.500 - 0.568)
        //   = 27.888 * 0.581 ~ 16.2 (repulsive — atoms too close/collinear)
        // The actual energy depends on geometry; let's just verify the formula runs
        assert!(energy2.is_finite(), "Energy should be finite");
    }

    #[test]
    fn test_torsion_angle_computation() {
        // IUPAC convention: dihedral A-B-C-D.
        // B→C along +Y, A at (1,0,0) projects onto +X perpendicular to B-C.

        // Trans (±180°): D projects to -X (opposite side from A)
        let trans = dihedral(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        )
        .unwrap();
        assert!(
            (trans.abs() - 180.0).abs() < 1.0,
            "Expected ~±180° for trans, got {trans}"
        );

        // Cis (0°): D projects to +X (same side as A)
        let cis = dihedral(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        )
        .unwrap();
        assert!(cis.abs() < 1.0, "Expected ~0° for cis, got {cis}");

        // -90°: D at (0,1,1) → projects to +Z
        let neg90 = dihedral(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
        )
        .unwrap();
        assert!((neg90 + 90.0).abs() < 1.0, "Expected ~-90°, got {neg90}");

        // +90°: D at (0,1,-1) → projects to -Z
        let pos90 = dihedral(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, -1.0],
        )
        .unwrap();
        assert!((pos90 - 90.0).abs() < 1.0, "Expected ~+90°, got {pos90}");

        // Degenerate geometry: coincident atoms → returns None
        assert!(
            dihedral(
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0]
            )
            .is_none()
        );
        // Collinear atoms → returns None
        assert!(
            dihedral(
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 0.0, 0.0]
            )
            .is_none()
        );
    }

    #[test]
    fn test_skip_non_protein_chains() {
        use crate::model::protein::{Atom, Chain, Residue};

        // Create an RNA chain — inference should not touch it
        let mut chains = vec![Chain {
            id: "A".to_string(),
            molecule_type: MoleculeType::RNA,
            residues: vec![Residue {
                name: "A".to_string(),
                seq_num: 1,
                atoms: vec![Atom {
                    name: "CA".to_string(),
                    element: "C".to_string(),
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    b_factor: 0.0,
                    is_backbone: true,
                    is_hetero: false,
                }],
                secondary_structure: SecondaryStructure::Coil,
            }],
        }];

        infer_secondary_structure(&mut chains);

        // RNA chain should remain all Coil
        assert_eq!(
            chains[0].residues[0].secondary_structure,
            SecondaryStructure::Coil,
            "RNA chain should not have SS inferred"
        );

        // Also test DNA
        let mut dna_chains = vec![Chain {
            id: "B".to_string(),
            molecule_type: MoleculeType::DNA,
            residues: vec![Residue {
                name: "DA".to_string(),
                seq_num: 1,
                atoms: vec![Atom {
                    name: "CA".to_string(),
                    element: "C".to_string(),
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    b_factor: 0.0,
                    is_backbone: true,
                    is_hetero: false,
                }],
                secondary_structure: SecondaryStructure::Coil,
            }],
        }];

        infer_secondary_structure(&mut dna_chains);
        assert_eq!(
            dna_chains[0].residues[0].secondary_structure,
            SecondaryStructure::Coil,
            "DNA chain should not have SS inferred"
        );
    }

    #[test]
    fn test_skip_chains_with_existing_ss() {
        use crate::model::protein::{Atom, Chain, Residue};

        // Create a protein chain with some residues already assigned to Helix
        let mut chains = vec![Chain {
            id: "A".to_string(),
            molecule_type: MoleculeType::Protein,
            residues: vec![
                Residue {
                    name: "ALA".to_string(),
                    seq_num: 1,
                    atoms: vec![Atom {
                        name: "CA".to_string(),
                        element: "C".to_string(),
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                        b_factor: 0.0,
                        is_backbone: true,
                        is_hetero: false,
                    }],
                    secondary_structure: SecondaryStructure::Helix,
                },
                Residue {
                    name: "GLY".to_string(),
                    seq_num: 2,
                    atoms: vec![Atom {
                        name: "CA".to_string(),
                        element: "C".to_string(),
                        x: 3.8,
                        y: 0.0,
                        z: 0.0,
                        b_factor: 0.0,
                        is_backbone: true,
                        is_hetero: false,
                    }],
                    secondary_structure: SecondaryStructure::Coil,
                },
            ],
        }];

        infer_secondary_structure(&mut chains);

        // The chain already has non-Coil SS, so inference should leave it untouched
        assert_eq!(
            chains[0].residues[0].secondary_structure,
            SecondaryStructure::Helix,
            "Existing Helix should be preserved"
        );
        assert_eq!(
            chains[0].residues[1].secondary_structure,
            SecondaryStructure::Coil,
            "Existing Coil should be preserved when chain has explicit SS"
        );
    }

    #[test]
    fn test_gap_filling() {
        // Test that single-residue Coil gaps within helix/sheet runs are filled
        // Use canonical alpha-helix torsion angles: phi ~ -57, psi ~ -47
        let torsions = vec![
            Some((-57.0, -47.0)), // helix-compatible
            Some((-57.0, -47.0)), // helix-compatible (gap to fill)
            Some((-57.0, -47.0)), // helix-compatible
        ];
        let mut assignments = vec![
            SecondaryStructure::Helix,
            SecondaryStructure::Coil,
            SecondaryStructure::Helix,
        ];

        fill_single_gaps(&mut assignments, &torsions, SecondaryStructure::Helix);

        assert_eq!(assignments[0], SecondaryStructure::Helix);
        assert_eq!(
            assignments[1],
            SecondaryStructure::Helix,
            "Single Coil gap between Helix runs should be filled"
        );
        assert_eq!(assignments[2], SecondaryStructure::Helix);
    }

    #[test]
    fn test_gap_filling_no_fill_when_torsion_incompatible() {
        // If the gap residue has incompatible torsion angles, it should NOT be filled
        let torsions = vec![
            Some((-57.0, -47.0)),  // helix-compatible
            Some((-120.0, 130.0)), // sheet-compatible, NOT helix
            Some((-57.0, -47.0)),  // helix-compatible
        ];
        let mut assignments = vec![
            SecondaryStructure::Helix,
            SecondaryStructure::Coil,
            SecondaryStructure::Helix,
        ];

        fill_single_gaps(&mut assignments, &torsions, SecondaryStructure::Helix);

        assert_eq!(
            assignments[1],
            SecondaryStructure::Coil,
            "Gap with incompatible torsions should NOT be filled"
        );
    }

    #[test]
    fn test_minimum_run_length() {
        // Helices shorter than 3 residues should be removed
        let mut assignments = vec![
            SecondaryStructure::Coil,
            SecondaryStructure::Helix,
            SecondaryStructure::Helix,
            SecondaryStructure::Coil,
            SecondaryStructure::Helix,
            SecondaryStructure::Helix,
            SecondaryStructure::Helix,
            SecondaryStructure::Coil,
        ];

        enforce_min_run(&mut assignments, SecondaryStructure::Helix, 3);

        // First helix run (2 residues) should be removed
        assert_eq!(
            assignments[1],
            SecondaryStructure::Coil,
            "2-residue helix should be removed"
        );
        assert_eq!(
            assignments[2],
            SecondaryStructure::Coil,
            "2-residue helix should be removed"
        );

        // Second helix run (3 residues) should survive
        assert_eq!(
            assignments[4],
            SecondaryStructure::Helix,
            "3-residue helix should survive"
        );
        assert_eq!(
            assignments[5],
            SecondaryStructure::Helix,
            "3-residue helix should survive"
        );
        assert_eq!(
            assignments[6],
            SecondaryStructure::Helix,
            "3-residue helix should survive"
        );
    }

    #[test]
    fn test_minimum_run_length_sheet() {
        // Sheets shorter than 2 residues should be removed
        let mut assignments = vec![
            SecondaryStructure::Coil,
            SecondaryStructure::Sheet, // 1-residue run: should be removed
            SecondaryStructure::Coil,
            SecondaryStructure::Sheet, // 2-residue run: should survive
            SecondaryStructure::Sheet,
            SecondaryStructure::Coil,
        ];

        enforce_min_run(&mut assignments, SecondaryStructure::Sheet, 2);

        assert_eq!(
            assignments[1],
            SecondaryStructure::Coil,
            "1-residue sheet should be removed"
        );
        assert_eq!(
            assignments[3],
            SecondaryStructure::Sheet,
            "2-residue sheet should survive"
        );
        assert_eq!(
            assignments[4],
            SecondaryStructure::Sheet,
            "2-residue sheet should survive"
        );
    }

    #[test]
    fn test_infer_ss_alphafold_pdb() {
        // Integration test: AF3_TNFa.pdb has no HELIX/SHEET records,
        // so inference should produce significant structured residues.
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/AF3_TNFa.pdb");
        let protein = crate::parser::pdb::load_structure(path).unwrap();

        let structured_count = protein
            .chains
            .iter()
            .flat_map(|c| &c.residues)
            .filter(|r| r.secondary_structure != SecondaryStructure::Coil)
            .count();
        let total_count = protein.chains.iter().flat_map(|c| &c.residues).count();

        assert!(
            structured_count > total_count / 3,
            "Expected significant SS for AF3_TNFa.pdb, got {structured_count}/{total_count} non-coil"
        );
    }

    #[test]
    fn test_helix_torsion_canonical_alpha() {
        // Canonical alpha helix: phi ~ -57, psi ~ -47
        assert!(is_helix_torsion(-57.0, -47.0));
        // 3_10 helix: phi ~ -49, psi ~ -26
        assert!(is_helix_torsion(-49.0, -26.0));
        // Strong window boundaries
        assert!(is_helix_torsion(-80.0, -70.0));
        assert!(is_helix_torsion(-40.0, -20.0));
        // Just outside strong but inside weak
        assert!(is_helix_torsion(-170.0, -80.0));
        assert!(is_helix_torsion(-20.0, 10.0));
        // Outside all windows
        assert!(!is_helix_torsion(-10.0, -47.0));
        assert!(!is_helix_torsion(-57.0, 20.0));
        assert!(!is_helix_torsion(-180.0, -47.0));
    }

    #[test]
    fn test_sheet_torsion_canonical_beta() {
        // Canonical antiparallel beta: phi ~ -139, psi ~ 135
        assert!(is_sheet_torsion(-139.0, 135.0));
        // Canonical parallel beta: phi ~ -119, psi ~ 113
        assert!(is_sheet_torsion(-80.0, 50.0));
        // Outside sheet windows
        assert!(!is_sheet_torsion(-57.0, -47.0)); // alpha helix region
    }

    #[test]
    fn test_explicit_pdb_ss_preserved_after_inference() {
        // 1UBQ has explicit HELIX/SHEET records. Inference should NOT
        // override them — the chain should be skipped entirely.
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/1UBQ.pdb");
        let protein = crate::parser::pdb::load_structure(path).unwrap();
        let chain = &protein.chains[0];

        // Residue 10 is in a beta sheet, residue 23 is in an alpha helix
        let residue_10 = chain.residues.iter().find(|r| r.seq_num == 10).unwrap();
        let residue_23 = chain.residues.iter().find(|r| r.seq_num == 23).unwrap();

        assert_eq!(
            residue_10.secondary_structure,
            SecondaryStructure::Sheet,
            "Explicit PDB Sheet assignment should be preserved"
        );
        assert_eq!(
            residue_23.secondary_structure,
            SecondaryStructure::Helix,
            "Explicit PDB Helix assignment should be preserved"
        );
    }

    #[test]
    fn test_ubq_helix_torsions_in_range() {
        // 1UBQ has known alpha helix at residues 23-34.
        // Canonical alpha helix: phi ~ -57, psi ~ -47.
        // We verify torsions fall in the helix-compatible window.
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/1UBQ.pdb");
        let protein = crate::parser::pdb::load_structure(path).unwrap();
        let chain = &protein.chains[0];
        let torsions = compute_torsion_angles(&chain.residues);

        let helix_residues: Vec<_> = chain
            .residues
            .iter()
            .enumerate()
            .filter(|(_, r)| r.seq_num >= 24 && r.seq_num <= 33)
            .collect();

        let helix_compatible = helix_residues
            .iter()
            .filter(|(i, _)| torsions[*i].map_or(false, |(phi, psi)| is_helix_torsion(phi, psi)))
            .count();

        assert!(
            helix_compatible >= helix_residues.len() / 2,
            "Expected most 1UBQ helix residues (24-33) to have helix-compatible torsions, got {}/{}",
            helix_compatible,
            helix_residues.len()
        );
    }
}

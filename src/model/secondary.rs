use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::model::protein::{MoleculeType, Protein, SecondaryStructure};

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
        Err(_) => return Vec::new(),
    };
    let reader = BufReader::new(file);
    let mut ranges = Vec::new();

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
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
                    end_seq: end_seq,
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
                    end_seq: end_seq,
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
        Err(_) => return Vec::new(),
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

/// Infer protein secondary structure directly from backbone geometry.
///
/// This is a fallback for structures such as AlphaFold/AlphaFold3 PDBs that
/// often omit HELIX/SHEET records entirely. We use phi/psi torsion angle
/// windows and keep only contiguous runs long enough to look like real
/// helices/sheets in the cartoon renderer.
pub fn infer_protein_secondary_structure(protein: &mut Protein) {
    for chain in &mut protein.chains {
        if chain.molecule_type != MoleculeType::Protein {
            continue;
        }

        let has_existing_ss = chain
            .residues
            .iter()
            .any(|r| r.secondary_structure != SecondaryStructure::Coil);
        if has_existing_ss {
            continue;
        }

        let inferred = infer_chain_secondary_structure(&chain.residues);
        for (residue, ss) in chain.residues.iter_mut().zip(inferred) {
            residue.secondary_structure = ss;
        }
    }
}

fn infer_chain_secondary_structure(
    residues: &[crate::model::protein::Residue],
) -> Vec<SecondaryStructure> {
    let mut assignments = vec![SecondaryStructure::Coil; residues.len()];
    let torsions = compute_torsions(residues);
    let hbonds = compute_hbond_map(residues);

    assign_helices_from_hbonds(&mut assignments, &hbonds, &torsions);
    assign_sheets_from_hbonds(&mut assignments, &hbonds, &torsions);
    fill_single_residue_gaps(&mut assignments, &torsions, SecondaryStructure::Helix);
    fill_single_residue_gaps(&mut assignments, &torsions, SecondaryStructure::Sheet);
    retain_runs(&mut assignments, SecondaryStructure::Helix, 3);
    retain_runs(&mut assignments, SecondaryStructure::Sheet, 2);
    assignments
}

fn compute_torsions(residues: &[crate::model::protein::Residue]) -> Vec<Option<(f64, f64)>> {
    let mut torsions = vec![None; residues.len()];
    for i in 1..residues.len().saturating_sub(1) {
        let c_prev = atom_pos(&residues[i - 1], "C");
        let n = atom_pos(&residues[i], "N");
        let ca = atom_pos(&residues[i], "CA");
        let c = atom_pos(&residues[i], "C");
        let n_next = atom_pos(&residues[i + 1], "N");

        let (Some(c_prev), Some(n), Some(ca), Some(c), Some(n_next)) = (c_prev, n, ca, c, n_next)
        else {
            continue;
        };

        let phi = dihedral(c_prev, n, ca, c);
        let psi = dihedral(n, ca, c, n_next);
        torsions[i] = Some((phi, psi));
    }
    torsions
}

fn compute_hbond_map(residues: &[crate::model::protein::Residue]) -> Vec<Vec<bool>> {
    let n = residues.len();
    let mut hbonds = vec![vec![false; n]; n];

    let c_atoms: Vec<Option<[f64; 3]>> = residues.iter().map(|r| atom_pos(r, "C")).collect();
    let o_atoms: Vec<Option<[f64; 3]>> = residues.iter().map(|r| atom_pos(r, "O")).collect();
    let n_atoms: Vec<Option<[f64; 3]>> = residues.iter().map(|r| atom_pos(r, "N")).collect();
    let ca_atoms: Vec<Option<[f64; 3]>> = residues.iter().map(|r| atom_pos(r, "CA")).collect();

    for acceptor in 0..n {
        let (Some(c), Some(o)) = (c_atoms[acceptor], o_atoms[acceptor]) else {
            continue;
        };

        for donor in 0..n {
            if donor == acceptor || donor.abs_diff(acceptor) <= 1 {
                continue;
            }

            let (Some(n_atom), Some(ca_atom)) = (n_atoms[donor], ca_atoms[donor]) else {
                continue;
            };
            let Some(h_atom) =
                estimate_amide_h(&c_atoms, &n_atoms, &ca_atoms, donor, n_atom, ca_atom)
            else {
                continue;
            };

            let energy = hbond_energy(o, c, n_atom, h_atom);
            if energy < -0.5 {
                hbonds[acceptor][donor] = true;
            }
        }
    }

    hbonds
}

fn estimate_amide_h(
    c_atoms: &[Option<[f64; 3]>],
    _n_atoms: &[Option<[f64; 3]>],
    _ca_atoms: &[Option<[f64; 3]>],
    donor: usize,
    n_atom: [f64; 3],
    ca_atom: [f64; 3],
) -> Option<[f64; 3]> {
    let prev_c = donor.checked_sub(1).and_then(|i| c_atoms[i])?;
    let dir_prev = normalize(sub(n_atom, prev_c))?;
    let dir_ca = normalize(sub(n_atom, ca_atom))?;
    let bisector = normalize(add(dir_prev, dir_ca))?;
    Some(add(n_atom, scale(bisector, 1.0)))
}

fn hbond_energy(o: [f64; 3], c: [f64; 3], n: [f64; 3], h: [f64; 3]) -> f64 {
    let r_on = distance(o, n).max(0.5);
    let r_ch = distance(c, h).max(0.5);
    let r_oh = distance(o, h).max(0.5);
    let r_cn = distance(c, n).max(0.5);
    27.888 * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn)
}

fn assign_helices_from_hbonds(
    assignments: &mut [SecondaryStructure],
    hbonds: &[Vec<bool>],
    torsions: &[Option<(f64, f64)>],
) {
    let mut support = vec![0usize; assignments.len()];

    for turn in [4usize, 3, 5] {
        for i in 0..assignments.len().saturating_sub(turn) {
            if !hbonds[i][i + turn] {
                continue;
            }

            let span = i + 1..=i + turn;
            let compatible = span
                .clone()
                .filter(|&idx| torsions_match_target(torsions[idx], SecondaryStructure::Helix))
                .count();
            let span_len = turn;
            if compatible * 2 < span_len {
                continue;
            }

            for idx in span {
                support[idx] += 1;
            }
        }
    }

    for (idx, state) in assignments.iter_mut().enumerate() {
        if support[idx] > 0 {
            *state = SecondaryStructure::Helix;
        }
    }
}

fn assign_sheets_from_hbonds(
    assignments: &mut [SecondaryStructure],
    hbonds: &[Vec<bool>],
    torsions: &[Option<(f64, f64)>],
) {
    let mut support = vec![0usize; assignments.len()];
    let n = assignments.len();

    for i in 1..n.saturating_sub(1) {
        for j in i + 2..n.saturating_sub(1) {
            if !torsions_match_target(torsions[i], SecondaryStructure::Sheet)
                || !torsions_match_target(torsions[j], SecondaryStructure::Sheet)
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

    for idx in 0..n {
        if assignments[idx] == SecondaryStructure::Coil
            && support[idx] > 0
            && torsions_match_target(torsions[idx], SecondaryStructure::Sheet)
        {
            assignments[idx] = SecondaryStructure::Sheet;
        }
    }
}

fn atom_pos(residue: &crate::model::protein::Residue, atom_name: &str) -> Option<[f64; 3]> {
    residue
        .atoms
        .iter()
        .find(|atom| atom.name == atom_name)
        .map(|atom| [atom.x, atom.y, atom.z])
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

fn dihedral(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let b0 = sub(a, b);
    let b1 = sub(c, b);
    let b2 = sub(d, c);

    let Some(b1_unit) = normalize(b1) else {
        return 0.0;
    };

    let n0 = cross(b0, b1);
    let n1 = cross(b1, b2);
    let Some(n0_unit) = normalize(n0) else {
        return 0.0;
    };
    let Some(n1_unit) = normalize(n1) else {
        return 0.0;
    };

    let m1 = cross(n0_unit, b1_unit);
    dot(m1, n1_unit).atan2(dot(n0_unit, n1_unit)).to_degrees()
}

fn is_strong_helix_torsion(phi: f64, psi: f64) -> bool {
    (-140.0..=-80.0).contains(&phi) && (-170.0..=-100.0).contains(&psi)
}

fn is_weak_helix_torsion(phi: f64, psi: f64) -> bool {
    (-170.0..=-40.0).contains(&phi) && (-180.0..=-60.0).contains(&psi)
}

fn is_strong_sheet_torsion(phi: f64, psi: f64) -> bool {
    ((-100.0..=-40.0).contains(&phi) && (20.0..=90.0).contains(&psi))
        || ((80.0..=180.0).contains(&phi) && (120.0..=180.0).contains(&psi))
}

fn is_weak_sheet_torsion(phi: f64, psi: f64) -> bool {
    ((-140.0..=-20.0).contains(&phi) && (0.0..=180.0).contains(&psi))
        || ((60.0..=180.0).contains(&phi) && (90.0..=180.0).contains(&psi))
}

fn fill_single_residue_gaps(
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

        if torsions_match_target(torsions[i], target) || torsions[i].is_none() {
            assignments[i] = target;
        }
    }
}

fn torsions_match_target(torsions: Option<(f64, f64)>, target: SecondaryStructure) -> bool {
    let Some((phi, psi)) = torsions else {
        return false;
    };
    match target {
        SecondaryStructure::Helix => {
            is_strong_helix_torsion(phi, psi) || is_weak_helix_torsion(phi, psi)
        }
        SecondaryStructure::Sheet => {
            is_strong_sheet_torsion(phi, psi) || is_weak_sheet_torsion(phi, psi)
        }
        _ => false,
    }
}

fn retain_runs(assignments: &mut [SecondaryStructure], ss: SecondaryStructure, min_len: usize) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::pdb::load_structure;

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

    #[test]
    fn test_infer_secondary_structure_for_alphafold_pdb() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/AF3_TNFa.pdb");
        let protein = load_structure(path).unwrap();

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
        let structured_count = helix_count + sheet_count;

        assert!(
            helix_count > 0,
            "expected inferred helices for AF3_TNFa.pdb"
        );
        assert!(
            structured_count > 0,
            "expected inferred secondary structure for AF3_TNFa.pdb"
        );
    }

    #[test]
    fn test_explicit_pdb_secondary_structure_still_used() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/1UBQ.pdb");
        let protein = load_structure(path).unwrap();
        let chain = &protein.chains[0];

        let residue_10 = chain.residues.iter().find(|r| r.seq_num == 10).unwrap();
        let residue_23 = chain.residues.iter().find(|r| r.seq_num == 23).unwrap();

        assert_eq!(residue_10.secondary_structure, SecondaryStructure::Sheet);
        assert_eq!(residue_23.secondary_structure, SecondaryStructure::Helix);
    }
}

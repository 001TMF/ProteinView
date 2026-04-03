/// Covalent-radii-based bond detection.
///
/// Two atoms are considered bonded when their distance is less than
/// `1.3 × (r₁ + r₂)` where r₁, r₂ are the single-bond covalent radii
/// (Cordero et al., Dalton Trans. 2008).  The 1.3 scale factor matches
/// the convention used by Open Babel, ASE, and other molecular toolkits.
const BOND_SCALE: f64 = 1.3;

/// Look up the single-bond covalent radius (Å) for a given element symbol.
/// Falls back to 1.5 Å for unknown elements (reasonable for transition metals).
pub fn covalent_radius(element: &str) -> f64 {
    let e = element.trim();
    match e {
        // Period 1
        "H" | "h" => 0.31,
        "He" | "HE" | "he" => 0.28,
        // Period 2
        "Li" | "LI" | "li" => 1.28,
        "Be" | "BE" | "be" => 0.96,
        "B" | "b" => 0.84,
        "C" | "c" => 0.76,
        "N" | "n" => 0.71,
        "O" | "o" => 0.66,
        "F" | "f" => 0.57,
        "Ne" | "NE" | "ne" => 0.58,
        // Period 3
        "Na" | "NA" | "na" => 1.66,
        "Mg" | "MG" | "mg" => 1.41,
        "Al" | "AL" | "al" => 1.21,
        "Si" | "SI" | "si" => 1.11,
        "P" | "p" => 1.07,
        "S" | "s" => 1.05,
        "Cl" | "CL" | "cl" => 1.02,
        "Ar" | "AR" | "ar" => 1.06,
        // Period 4
        "K" | "k" => 2.03,
        "Ca" | "CA" | "ca" => 1.76,
        "Sc" | "SC" | "sc" => 1.70,
        "Ti" | "TI" | "ti" => 1.60,
        "V" | "v" => 1.53,
        "Cr" | "CR" | "cr" => 1.39,
        "Mn" | "MN" | "mn" => 1.39,
        "Fe" | "FE" | "fe" => 1.32,
        "Co" | "CO" | "co" => 1.26,
        "Ni" | "NI" | "ni" => 1.24,
        "Cu" | "CU" | "cu" => 1.32,
        "Zn" | "ZN" | "zn" => 1.22,
        "Ga" | "GA" | "ga" => 1.22,
        "Ge" | "GE" | "ge" => 1.20,
        "As" | "AS" | "as" => 1.19,
        "Se" | "SE" | "se" => 1.20,
        "Br" | "BR" | "br" => 1.20,
        "Kr" | "KR" | "kr" => 1.16,
        // Period 5
        "Rb" | "RB" | "rb" => 2.20,
        "Sr" | "SR" | "sr" => 1.95,
        "Y" | "y" => 1.90,
        "Zr" | "ZR" | "zr" => 1.75,
        "Nb" | "NB" | "nb" => 1.64,
        "Mo" | "MO" | "mo" => 1.54,
        "Ru" | "RU" | "ru" => 1.46,
        "Rh" | "RH" | "rh" => 1.42,
        "Pd" | "PD" | "pd" => 1.39,
        "Ag" | "AG" | "ag" => 1.45,
        "Cd" | "CD" | "cd" => 1.44,
        "In" | "IN" | "in" => 1.42,
        "Sn" | "SN" | "sn" => 1.39,
        "Sb" | "SB" | "sb" => 1.39,
        "Te" | "TE" | "te" => 1.38,
        "I" | "i" => 1.39,
        "Xe" | "XE" | "xe" => 1.40,
        // Period 6 (common)
        "Cs" | "CS" | "cs" => 2.44,
        "Ba" | "BA" | "ba" => 2.15,
        "W" | "w" => 1.62,
        "Re" | "RE" | "re" => 1.51,
        "Os" | "OS" | "os" => 1.44,
        "Ir" | "IR" | "ir" => 1.41,
        "Pt" | "PT" | "pt" => 1.36,
        "Au" | "AU" | "au" => 1.36,
        "Hg" | "HG" | "hg" => 1.32,
        "Tl" | "TL" | "tl" => 1.45,
        "Pb" | "PB" | "pb" => 1.46,
        "Bi" | "BI" | "bi" => 1.48,
        // Lanthanides (representative)
        "La" | "LA" | "la" => 2.07,
        "Ce" | "CE" | "ce" => 2.04,
        "Eu" | "EU" | "eu" => 1.98,
        "Gd" | "GD" | "gd" => 1.96,
        "Lu" | "LU" | "lu" => 1.87,
        // Actinides (representative)
        "U" | "u" => 1.96,
        // Fallback
        _ => 1.50,
    }
}

/// Check whether two atoms are likely bonded based on covalent radii.
///
/// Returns true when distance < BOND_SCALE × (r₁ + r₂).
#[allow(clippy::too_many_arguments)]
pub fn atoms_bonded(
    elem1: &str,
    x1: f64,
    y1: f64,
    z1: f64,
    elem2: &str,
    x2: f64,
    y2: f64,
    z2: f64,
) -> bool {
    let dx = x2 - x1;
    let dy = y2 - y1;
    let dz = z2 - z1;
    let dist_sq = dx * dx + dy * dy + dz * dz;
    let r_sum = covalent_radius(elem1) + covalent_radius(elem2);
    let cutoff = BOND_SCALE * r_sum;
    dist_sq <= cutoff * cutoff
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_carbon_carbon_single_bond() {
        // Typical C-C single bond: ~1.54 Å
        // Cutoff: 1.3 * (0.76 + 0.76) = 1.976 Å
        assert!(atoms_bonded("C", 0.0, 0.0, 0.0, "C", 1.54, 0.0, 0.0));
    }

    #[test]
    fn test_carbon_carbon_too_far() {
        assert!(!atoms_bonded("C", 0.0, 0.0, 0.0, "C", 2.5, 0.0, 0.0));
    }

    #[test]
    fn test_carbon_hydrogen_bond() {
        // Typical C-H: ~1.09 Å, cutoff: 1.3 * (0.76 + 0.31) = 1.391 Å
        assert!(atoms_bonded("C", 0.0, 0.0, 0.0, "H", 1.09, 0.0, 0.0));
    }

    #[test]
    fn test_hydrogen_hydrogen_not_bonded() {
        // Cutoff: 1.3 * (0.31 + 0.31) = 0.806 Å
        assert!(!atoms_bonded("H", 0.0, 0.0, 0.0, "H", 1.0, 0.0, 0.0));
    }

    #[test]
    fn test_metal_ligand_bond() {
        // Fe-N bond in heme: ~2.0 Å
        // Cutoff: 1.3 * (1.32 + 0.71) = 2.639 Å
        assert!(atoms_bonded("Fe", 0.0, 0.0, 0.0, "N", 2.0, 0.0, 0.0));
    }

    #[test]
    fn test_metal_ligand_not_bonded() {
        assert!(!atoms_bonded("Fe", 0.0, 0.0, 0.0, "N", 3.0, 0.0, 0.0));
    }

    #[test]
    fn test_case_insensitive() {
        assert!(atoms_bonded("c", 0.0, 0.0, 0.0, "C", 1.54, 0.0, 0.0));
        assert!(atoms_bonded("FE", 0.0, 0.0, 0.0, "n", 2.0, 0.0, 0.0));
    }

    #[test]
    fn test_unknown_element_uses_fallback() {
        // Cutoff: 1.3 * (1.5 + 1.5) = 3.9 Å
        assert!(atoms_bonded("Xx", 0.0, 0.0, 0.0, "Yy", 3.5, 0.0, 0.0));
        assert!(!atoms_bonded("Xx", 0.0, 0.0, 0.0, "Yy", 4.5, 0.0, 0.0));
    }

    #[test]
    fn test_old_cutoff_missed_fe_n() {
        // The old fixed 1.9 Å cutoff would miss this Fe-N bond at 2.0 Å
        let old_cutoff_sq = 1.9_f64 * 1.9;
        let dist_sq = 2.0_f64 * 2.0;
        assert!(dist_sq > old_cutoff_sq, "old cutoff would miss this bond");
        assert!(atoms_bonded("Fe", 0.0, 0.0, 0.0, "N", 2.0, 0.0, 0.0));
    }

    #[test]
    fn test_water_bonds() {
        // O-H bond: ~0.96 Å, cutoff: 1.3 * (0.66 + 0.31) = 1.261 Å — bonded
        assert!(atoms_bonded("O", 0.0, 0.0, 0.117, "H", 0.0, 0.757, -0.469));
        // H-H in water: ~1.51 Å, cutoff: 0.806 Å — NOT bonded
        assert!(!atoms_bonded(
            "H", 0.0, 0.757, -0.469, "H", 0.0, -0.757, -0.469
        ));
    }
}

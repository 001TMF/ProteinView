use super::{looks_like_mmcif, normalize_pdb_id};

#[test]
fn normalizes_a_valid_rcsb_identifier() {
    assert_eq!(normalize_pdb_id("1ubq").unwrap(), "1UBQ");
}

#[test]
fn rejects_identifiers_that_could_escape_the_cache_path() {
    for invalid in ["../x", "1/../", "1ABC.cif", "ABCD", "0ABC", "1 A?"] {
        assert!(normalize_pdb_id(invalid).is_err(), "{invalid}");
    }
}

#[test]
fn recognizes_mmcif_and_rejects_html() {
    assert!(looks_like_mmcif(
        b"# comment\n\ndata_1UBQ\n_entry.id 1UBQ\n"
    ));
    assert!(!looks_like_mmcif(b"<html>temporary error</html>"));
    assert!(!looks_like_mmcif(b"\xff\xfe"));
}

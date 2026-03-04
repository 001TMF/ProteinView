use anyhow::Result;

/// Download a PDB file from RCSB by ID (requires "fetch" feature)
pub fn fetch_pdb(pdb_id: &str) -> Result<String> {
    #[cfg(feature = "fetch")]
    {
        let url = format!("https://files.rcsb.org/download/{}.cif", pdb_id.to_uppercase());
        let response = reqwest::blocking::get(&url)?;
        if !response.status().is_success() {
            anyhow::bail!("Failed to fetch PDB {}: HTTP {}", pdb_id, response.status());
        }
        let tmp_path = std::env::temp_dir().join(format!("{}.cif", pdb_id.to_uppercase()));
        std::fs::write(&tmp_path, response.bytes()?)?;
        Ok(tmp_path.to_string_lossy().to_string())
    }
    #[cfg(not(feature = "fetch"))]
    {
        let _ = pdb_id;
        anyhow::bail!("Fetch feature not enabled. Build with: cargo build --features fetch")
    }
}

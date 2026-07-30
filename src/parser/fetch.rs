use anyhow::Result;

#[cfg(feature = "fetch")]
const MAX_DOWNLOAD_BYTES: u64 = 64 * 1024 * 1024;

/// Download a PDB file from RCSB by ID (requires "fetch" feature)
pub fn fetch_pdb(pdb_id: &str) -> Result<String> {
    let normalized_id = normalize_pdb_id(pdb_id)?;

    #[cfg(feature = "fetch")]
    {
        use std::io::{Read, Write};
        use std::time::Duration;

        let url = format!("https://files.rcsb.org/download/{normalized_id}.cif");
        let client = reqwest::blocking::Client::builder()
            .timeout(Duration::from_secs(30))
            .redirect(reqwest::redirect::Policy::none())
            .build()?;
        let response = client.get(&url).send()?;
        if !response.status().is_success() {
            anyhow::bail!(
                "Failed to fetch PDB {normalized_id}: HTTP {}",
                response.status()
            );
        }
        if response
            .content_length()
            .is_some_and(|length| length > MAX_DOWNLOAD_BYTES)
        {
            anyhow::bail!(
                "PDB {normalized_id} exceeds the {MAX_DOWNLOAD_BYTES}-byte download limit"
            );
        }
        if response
            .headers()
            .get(reqwest::header::CONTENT_TYPE)
            .and_then(|value| value.to_str().ok())
            .is_some_and(|value| value.to_ascii_lowercase().contains("text/html"))
        {
            anyhow::bail!("RCSB returned HTML instead of PDB {normalized_id}");
        }

        let mut bytes = Vec::new();
        response
            .take(MAX_DOWNLOAD_BYTES + 1)
            .read_to_end(&mut bytes)?;
        if bytes.len() as u64 > MAX_DOWNLOAD_BYTES {
            anyhow::bail!(
                "PDB {normalized_id} exceeds the {MAX_DOWNLOAD_BYTES}-byte download limit"
            );
        }
        if !looks_like_mmcif(&bytes) {
            anyhow::bail!("RCSB returned invalid mmCIF data for PDB {normalized_id}");
        }

        let mut temporary = tempfile::Builder::new()
            .prefix(&format!("proteinview-{normalized_id}-"))
            .suffix(".cif")
            .tempfile()?;
        temporary.write_all(&bytes)?;
        temporary.flush()?;
        temporary.as_file().sync_all()?;
        let (_file, tmp_path) = temporary.keep()?;
        Ok(tmp_path.to_string_lossy().to_string())
    }
    #[cfg(not(feature = "fetch"))]
    {
        let _ = normalized_id;
        anyhow::bail!("Fetch feature not enabled. Build with: cargo build --features fetch")
    }
}

#[cfg(any(feature = "fetch", test))]
fn looks_like_mmcif(bytes: &[u8]) -> bool {
    let Ok(text) = std::str::from_utf8(bytes) else {
        return false;
    };
    text.trim_start_matches(['\u{feff}', ' ', '\t', '\r', '\n'])
        .lines()
        .map(str::trim)
        .find(|line| !line.is_empty() && !line.starts_with('#'))
        .is_some_and(|line| line.to_ascii_lowercase().starts_with("data_"))
}

fn normalize_pdb_id(pdb_id: &str) -> Result<String> {
    let bytes = pdb_id.as_bytes();
    if bytes.len() != 4
        || !(b'1'..=b'9').contains(&bytes[0])
        || !bytes.iter().all(u8::is_ascii_alphanumeric)
    {
        anyhow::bail!(
            "invalid RCSB PDB ID `{pdb_id}`: expected four ASCII alphanumeric characters beginning with 1-9"
        );
    }
    Ok(pdb_id.to_ascii_uppercase())
}

#[cfg(test)]
#[path = "fetch_tests.rs"]
mod tests;

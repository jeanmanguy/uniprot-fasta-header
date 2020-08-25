use crate::{
    error::UniProtHeaderError,
    parser::{
        chevron, db, entry_name, iso_id, optional_gene_name, organism_id, os_until_ox, pipe, space,
        until_os,
    },
    Database,
};
use nom::{error::ParseError, IResult};

/// UniProtKB isoform header
#[derive(Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct UniProtKBIsoform {
    /// UniProtKB database
    pub database: Database,
    /// UniProt accession number (https://www.uniprot.org/help/accession%5Fnumbers)
    pub identifier: String,
    /// Isoform number
    pub isoform: String,
    /// UniProt entry name (https://www.uniprot.org/help/entry%5Fname)
    pub entry_name: String,
    /// Protein name (https://www.uniprot.org/help/protein%5Fnames)
    pub protein_name: String,
    /// Organism name
    pub organism_name: String,
    /// NCBI iaxonomic identifier (https://www.uniprot.org/help/taxonomic%5Fidentifier)
    pub organism_identifier: String,
    /// Gene name (https://www.uniprot.org/help/gene%5Fname)
    pub gene_name: Option<String>,
}

impl Default for UniProtKBIsoform {
    fn default() -> Self {
        Self {
            database: Database::SwissProt,
            identifier: String::default(),
            isoform: String::default(),
            entry_name: String::default(),
            protein_name: String::default(),
            organism_name: String::default(),
            organism_identifier: String::default(),
            gene_name: None,
        }
    }
}

/// Parse a UniProtKB isoform fasta header
pub fn uniprotkb_iso(string: &[u8]) -> Result<UniProtKBIsoform, UniProtHeaderError> {
    match parse_uniprotkb_iso(string) {
        Ok((_, parsed)) => Ok(parsed),
        Err(err) => match err {
            nom::Err::Incomplete(_i) => Err(UniProtHeaderError::Incomplete),
            nom::Err::Error((rest, kind)) => Err(UniProtHeaderError::from_error_kind(rest, kind)),
            nom::Err::Failure((rest, kind)) => Err(UniProtHeaderError::from_error_kind(rest, kind)),
        },
    }
}

fn parse_uniprotkb_iso(input: &[u8]) -> IResult<&[u8], UniProtKBIsoform> {
    let (input, _) = chevron(input)?;
    let (input, database) = db(input)?;
    let (input, _) = pipe(input)?;
    let (input, (id, iso)) = iso_id(input)?;
    let (input, _) = pipe(input)?;
    let (input, entry) = entry_name(input)?;
    let (input, _) = space(input)?;
    let (input, protein) = until_os(input)?;
    let (input, _) = space(input)?;
    let (input, organism) = os_until_ox(input)?;
    let (input, _) = space(input)?;
    let (input, organism_id) = organism_id(input)?;
    let (input, _) = space(input)?;
    let (input, (gene, _)) = optional_gene_name(input)?; // + optional space

    let gene_name = match gene {
        Some(g) => Some(String::from_utf8_lossy(g).to_string()),
        None => None,
    };

    Ok((
        input,
        UniProtKBIsoform {
            database,
            identifier: String::from_utf8_lossy(id).to_string(),
            isoform: String::from_utf8_lossy(iso).to_string(),
            entry_name: String::from_utf8(entry).unwrap(),
            protein_name: String::from_utf8_lossy(protein).trim().to_string(),
            organism_name: String::from_utf8_lossy(organism).trim().to_string(),
            organism_identifier: String::from_utf8_lossy(organism_id).trim().to_string(),
            gene_name,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_uniprot_1433b_macfa_iso() {
        let entry = UniProtKBIsoform {
            database: Database::SwissProt,
            identifier: "Q4R572".to_string(),
            isoform: "2".to_string(),
            entry_name: "1433B_MACFA".to_string(),
            protein_name: "Isoform Short of 14-3-3 protein beta/alpha".to_string(),
            organism_name: "Macaca fascicularis".to_string(),
            organism_identifier: "9541".to_string(),
            gene_name: Some("YWHAB".to_string()),
        };
        let test_header =
            ">sp|Q4R572-2|1433B_MACFA Isoform Short of 14-3-3 protein beta/alpha OS=Macaca fascicularis OX=9541 GN=YWHAB".as_bytes();
        assert_eq!(uniprotkb_iso(test_header).unwrap(), entry)
    }

    #[test]
    fn test_uniprot_alg2_human_iso() {
        let entry = UniProtKBIsoform {
            database: Database::SwissProt,
            identifier: "Q9H553".to_string(),
            isoform: "2".to_string(),
            entry_name: "ALG2_HUMAN".to_string(),
            protein_name: "Isoform 2 of Alpha-1,3/1,6-mannosyltransferase ALG2".to_string(),
            organism_name: "Homo sapiens".to_string(),
            organism_identifier: "9606".to_string(),
            gene_name: Some("ALG2".to_string()),
        };
        let test_header =
            ">sp|Q9H553-2|ALG2_HUMAN Isoform 2 of Alpha-1,3/1,6-mannosyltransferase ALG2 OS=Homo sapiens OX=9606 GN=ALG2".as_bytes();
        assert_eq!(uniprotkb_iso(test_header).unwrap(), entry)
    }

    #[test]
    fn test_uniprot_agl27_arath_iso() {
        let entry = UniProtKBIsoform {
            database: Database::SwissProt,
            identifier: "Q9AT76".to_string(),
            isoform: "4".to_string(),
            entry_name: "AGL27_ARATH".to_string(),
            protein_name: "Isoform 4 of Agamous-like MADS-box protein AGL27".to_string(),
            organism_name: "Arabidopsis thaliana".to_string(),
            organism_identifier: "3702".to_string(),
            gene_name: Some("AGL27".to_string()),
        };
        let test_header =
            ">sp|Q9AT76-4|AGL27_ARATH Isoform 4 of Agamous-like MADS-box protein AGL27 OS=Arabidopsis thaliana OX=3702 GN=AGL27".as_bytes();
        assert_eq!(uniprotkb_iso(test_header).unwrap(), entry)
    }

    #[test]
    fn test_uniprot_ters_bpspp_iso() {
        let entry = UniProtKBIsoform {
            database: Database::SwissProt,
            identifier: "P54307".to_string(),
            isoform: "2".to_string(),
            entry_name: "TERS_BPSPP".to_string(),
            protein_name: "Isoform G1P* of Terminase small subunit".to_string(),
            organism_name: "Bacillus phage SPP1".to_string(),
            organism_identifier: "10724".to_string(),
            gene_name: Some("1".to_string()),
        };
        let test_header =
            ">sp|P54307-2|TERS_BPSPP Isoform G1P* of Terminase small subunit OS=Bacillus phage SPP1 OX=10724 GN=1".as_bytes();
        assert_eq!(uniprotkb_iso(test_header).unwrap(), entry)
    }
}

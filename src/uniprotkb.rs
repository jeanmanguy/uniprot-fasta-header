use crate::{
    error::UniProtHeaderError,
    parser::{
        chevron, db, entry_name, evidence, optional_gene_name, organism_id, os_until_ox, pipe,
        space, unique_id, until_os, version,
    },
    Database, ProteinExistence,
};
use nom::{error::ParseError, IResult};

/// UniProtKB header
#[derive(Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct UniProtKB {
    /// UniProtKB database
    pub database: Database,
    /// UniProt accession number (https://www.uniprot.org/help/accession%5Fnumbers)
    pub identifier: String,
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
    /// Protein existence (https://www.uniprot.org/help/protein%5Fexistence)
    pub protein_existence: ProteinExistence,
    /// Sequence version (https://www.uniprot.org/help/entry%5Fhistory)
    pub sequence_version: String,
}

impl Default for UniProtKB {
    fn default() -> Self {
        Self {
            database: Database::SwissProt,
            identifier: String::default(),
            entry_name: String::default(),
            protein_name: String::default(),
            organism_name: String::default(),
            organism_identifier: String::default(),
            gene_name: None,
            protein_existence: ProteinExistence::Uncertain,
            sequence_version: String::default(),
        }
    }
}

/// Parse a UniProtKB fasta header
pub fn uniprotkb(string: &[u8]) -> Result<UniProtKB, UniProtHeaderError> {
    match parse_uniprotkb(string) {
        Ok((_, parsed)) => Ok(parsed),
        Err(err) => match err {
            nom::Err::Incomplete(_i) => Err(UniProtHeaderError::Incomplete),
            nom::Err::Error((rest, kind)) => Err(UniProtHeaderError::from_error_kind(rest, kind)),
            nom::Err::Failure((rest, kind)) => Err(UniProtHeaderError::from_error_kind(rest, kind)),
        },
    }
}

fn parse_uniprotkb(input: &[u8]) -> IResult<&[u8], UniProtKB> {
    let (input, _) = chevron(input)?;
    let (input, database) = db(input)?;
    let (input, _) = pipe(input)?;
    let (input, id) = unique_id(input)?;
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
    let (input, evidence) = evidence(input)?;
    let (input, _) = space(input)?;
    let (input, version) = version(input)?;

    let gene_name = match gene {
        Some(g) => Some(String::from_utf8_lossy(g).to_string()),
        None => None,
    };

    Ok((
        input,
        UniProtKB {
            database,
            identifier: String::from_utf8_lossy(id).trim().to_string(),
            entry_name: String::from_utf8(entry).unwrap(),
            protein_name: String::from_utf8_lossy(protein).trim().to_string(),
            organism_name: String::from_utf8_lossy(organism).trim().to_string(),
            organism_identifier: String::from_utf8_lossy(organism_id).trim().to_string(),
            gene_name,
            protein_existence: evidence,
            sequence_version: String::from_utf8_lossy(version).trim().to_string(),
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_uniprot_acn2_acago() {
        let entry = UniProtKB {
            database: Database::SwissProt,
            identifier: "Q8I6R7".to_string(),
            entry_name: "ACN2_ACAGO".to_string(),
            protein_name: "Acanthoscurrin-2 (Fragment)".to_string(),
            organism_name: "Acanthoscurria gomesiana".to_string(),
            organism_identifier: "115339".to_string(),
            gene_name: Some("acantho2".to_string()),
            protein_existence: ProteinExistence::ExperimentalEvidenceProtein,
            sequence_version: "1".to_string(),
        };
        let test_header = ">sp|Q8I6R7|ACN2_ACAGO Acanthoscurrin-2 (Fragment) OS=Acanthoscurria gomesiana OX=115339 GN=acantho2 PE=1 SV=1".as_bytes();
        assert_eq!(uniprotkb(test_header).unwrap(), entry)
    }

    #[test]
    fn test_uniprot_acox_cupnh() {
        let entry = UniProtKB {
            database: Database::SwissProt,
            identifier: "P27748".to_string(),
            entry_name: "ACOX_CUPNH".to_string(),
            protein_name: "Acetoin catabolism protein X".to_string(),
            organism_name: "Cupriavidus necator (strain ATCC 17699 / H16 / DSM 428 / Stanier 337)"
                .to_string(),
            organism_identifier: "381666".to_string(),
            gene_name: Some("acoX".to_string()),
            protein_existence: ProteinExistence::Predicted,
            sequence_version: "2".to_string(),
        };
        let test_header = ">sp|P27748|ACOX_CUPNH Acetoin catabolism protein X OS=Cupriavidus necator (strain ATCC 17699 / H16 / DSM 428 / Stanier 337) OX=381666 GN=acoX PE=4 SV=2".as_bytes();
        assert_eq!(uniprotkb(test_header).unwrap(), entry)
    }

    #[test]
    fn test_uniprot_ha22_mouse() {
        let entry = UniProtKB {
            database: Database::SwissProt,
            identifier: "P04224".to_string(),
            entry_name: "HA22_MOUSE".to_string(),
            protein_name: "H-2 class II histocompatibility antigen, E-K alpha chain".to_string(),
            organism_name: "Mus musculus".to_string(),
            organism_identifier: "10090".to_string(),
            gene_name: None,
            protein_existence: ProteinExistence::ExperimentalEvidenceProtein,
            sequence_version: "1".to_string(),
        };
        let test_header = ">sp|P04224|HA22_MOUSE H-2 class II histocompatibility antigen, E-K alpha chain OS=Mus musculus OX=10090 PE=1 SV=1".as_bytes();
        assert_eq!(uniprotkb(test_header).unwrap(), entry)
    }

    #[test]
    fn test_uniprot_q3sa23_9hiv1() {
        let entry = UniProtKB {
            database: Database::TrEMBL,
            identifier: "Q3SA23".to_string(),
            entry_name: "Q3SA23_9HIV1".to_string(),
            protein_name: "Protein Nef (Fragment)".to_string(),
            organism_name: "Human immunodeficiency virus 1".to_string(),
            organism_identifier: "11676".to_string(),
            gene_name: Some("nef".to_string()),
            protein_existence: ProteinExistence::InferredHomology,
            sequence_version: "1".to_string(),
        };
        let test_header = ">tr|Q3SA23|Q3SA23_9HIV1 Protein Nef (Fragment) OS=Human immunodeficiency virus 1  OX=11676 GN=nef PE=3 SV=1".as_bytes();
        assert_eq!(uniprotkb(test_header).unwrap(), entry)
    }

    #[test]
    fn test_uniprot_cask_bovin() {
        let entry = UniProtKB {
            database: Database::SwissProt,
            identifier: "P02668".to_string(),
            entry_name: "CASK_BOVIN".to_string(),
            protein_name: "Kappa-casein".to_string(),
            organism_name: "Bos taurus".to_string(),
            organism_identifier: "9913".to_string(),
            gene_name: Some("CSN3".to_string()),
            protein_existence: ProteinExistence::ExperimentalEvidenceProtein,
            sequence_version: "1".to_string(),
        };
        let test_header =
            ">sp|P02668|CASK_BOVIN Kappa-casein OS=Bos taurus OX=9913 GN=CSN3 PE=1 SV=1".as_bytes();
        assert_eq!(uniprotkb(test_header).unwrap(), entry)
    }

    #[test]
    fn test_uniprot_ypfu_ecoli() {
        let entry = UniProtKB {
            database: Database::SwissProt,
            identifier: "P18355".to_string(),
            entry_name: "YPFU_ECOLI".to_string(),
            protein_name: "Uncharacterized protein in traD-traI intergenic region".to_string(),
            organism_name: "Escherichia coli (strain K12)".to_string(),
            organism_identifier: "83333".to_string(),
            gene_name: None,
            protein_existence: ProteinExistence::InferredHomology,
            sequence_version: "1".to_string(),
        };
        let test_header =
            ">sp|P18355|YPFU_ECOLI Uncharacterized protein in traD-traI intergenic region OS=Escherichia coli (strain K12) OX=83333 PE=3 SV=1".as_bytes();
        assert_eq!(uniprotkb(test_header).unwrap(), entry)
    }
}

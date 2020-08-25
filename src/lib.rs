//! Parse UniProt Fasta Header
//!
//! ## UniProtKB
//!
//! ### Format
//!
//! `>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion`
//!
//!
//! ### Usage
//!
//! ```rust
//! let header = ">sp|P18355|YPFU_ECOLI Uncharacterized protein in traD-traI intergenic region OS=Escherichia coli (strain K12) OX=83333 PE=3 SV=1".as_bytes();
//!
//! let parsed = uniprot_fasta_header::uniprotkb(header).unwrap();
//!
//! let expected = uniprot_fasta_header::UniProtKB {
//!     database: uniprot_fasta_header::Database::SwissProt,
//!     identifier: "P18355".to_string(),
//!     entry_name: "YPFU_ECOLI".to_string(),
//!     protein_name: "Uncharacterized protein in traD-traI intergenic region".to_string(),
//!     organism_name: "Escherichia coli (strain K12)".to_string(),
//!     organism_identifier: "83333".to_string(),
//!     gene_name: None,
//!     protein_existence: uniprot_fasta_header::ProteinExistence::InferredHomology,
//!     sequence_version: "1".to_string(),
//! };
//!
//! assert_eq!(parsed, expected);
//! ```
//!
//!
//! ## UniProtKB isoform
//!
//! ### Format
//!
//! `>sp|IsoID|EntryName Isoform IsoformName of ProteinName OS=OrganismName OX=OrganismIdentifier[ GN=GeneName]`
//!
//! ### Usage
//!
//! ```rust
//! let header = ">sp|Q4R572-2|1433B_MACFA Isoform Short of 14-3-3 protein beta/alpha OS=Macaca fascicularis OX=9541 GN=YWHAB".as_bytes();
//!
//! let parsed = uniprot_fasta_header::uniprotkb_iso(header).unwrap();
//!
//! let expected = uniprot_fasta_header::UniProtKBIsoform {
//!     database: uniprot_fasta_header::Database::SwissProt,
//!     identifier: "Q4R572".to_string(),
//!     isoform: "2".to_string(),
//!     entry_name: "1433B_MACFA".to_string(),
//!     protein_name: "Isoform Short of 14-3-3 protein beta/alpha".to_string(),
//!     organism_name: "Macaca fascicularis".to_string(),
//!     organism_identifier: "9541".to_string(),
//!     gene_name: Some("YWHAB".to_string()),
//! };
//!
//! assert_eq!(parsed, expected);
//! ```
//!
//! ## Reference
//!
//! - [UniProt Fasta header help page](https://www.uniprot.org/help/fasta-headers)
#![warn(missing_docs)]
#![allow(clippy::type_complexity)]

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;

mod error;
mod parser;
mod uniprotkb;
mod uniprotkb_isoform;

pub use uniprotkb::uniprotkb;
pub use uniprotkb::UniProtKB;
pub use uniprotkb_isoform::uniprotkb_iso;
pub use uniprotkb_isoform::UniProtKBIsoform;

/// UniProtKB database
#[derive(Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Database {
    /// UniProtKB/Swiss-Prot
    SwissProt,
    /// UniProtKB/TrEMBL
    TrEMBL,
}

/// Protein Existence types
///
/// See [Protein existence](https://www.uniprot.org/help/protein%5Fexistence).
#[derive(Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ProteinExistence {
    /// 1. Experimental evidence at protein level
    ExperimentalEvidenceProtein,
    /// 2. Experimental evidence at transcript level
    ExperimentalEvidenceTranscript,
    /// 3. Protein inferred from homology
    InferredHomology,
    /// 4. Protein predicted
    Predicted,
    /// 5. Protein uncertain
    Uncertain,
}

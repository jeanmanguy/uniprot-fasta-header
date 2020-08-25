# Uniprot-fasta-header

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

Parse [Uniprot fasta headers](https://www.uniprot.org/help/fasta-headers).

## Features

- UniProtKB header parsing
- UniProtKB isoform header parsing

## Usage

### UniProtKB

 ```rust
 use uniprot_fasta_header::{UniProtKB, UniprotFasta};

 let header = ">sp|P18355|YPFU_ECOLI Uncharacterized protein in traD-traI intergenic region OS=Escherichia coli (strain K12) OX=83333 PE=3 SV=1".as_bytes();

 let (_, parsed) = UniProtKB::parse(header).unwrap();

 let expected = UniProtKB {
     database: uniprot_fasta_header::Database::SwissProt,
     identifier: "P18355".to_string(),
     entry_name: "YPFU_ECOLI".to_string(),
     protein_name: "Uncharacterized protein in traD-traI intergenic region".to_string(),
     organism_name: "Escherichia coli (strain K12)".to_string(),
     organism_identifier: "83333".to_string(),
     gene_name: None,
     protein_existence: uniprot_fasta_header::ProteinExistence::InferredHomology,
     sequence_version: "1".to_string(),
 };

 assert_eq!(parsed, expected);
 ```

### UniProtKB isoform

 ```rust
 use uniprot_fasta_header::{UniProtKBIsoform, UniprotFasta};

 let header = ">sp|Q4R572-2|1433B_MACFA Isoform Short of 14-3-3 protein beta/alpha OS=Macaca fascicularis OX=9541 GN=YWHAB".as_bytes();

 let (_, parsed) = UniProtKBIsoform::parse(header).unwrap();

 let expected = UniProtKBIsoform {
     database: uniprot_fasta_header::Database::SwissProt,
     identifier: "Q4R572".to_string(),
     isoform: "2".to_string(),
     entry_name: "1433B_MACFA".to_string(),
     protein_name: "Isoform Short of 14-3-3 protein beta/alpha".to_string(),
     organism_name: "Macaca fascicularis".to_string(),
     organism_identifier: "9541".to_string(),
     gene_name: Some("YWHAB".to_string()),
 };
```

## Ideas & bugs

Please report new bugs and weird Fasta headers by creating a new issue on the [project repository](https://github.com/jeanmanguy/uniprot-fasta-header/issues). 

## License

Uniprot-fasta-header is distributed under the terms of the Apache License (Version 2.0). See [LICENSE](./LICENSE) for details.
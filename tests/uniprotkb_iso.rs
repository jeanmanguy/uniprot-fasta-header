use std::fs::File;
use std::io::{self, BufRead};

#[test]
fn agl27_arath() {
    let file = File::open("tests/assets/agl27_arath_iso.txt").unwrap();
    let reader = io::BufReader::new(file);
    for (i, line) in reader.lines().enumerate() {
        if let Ok(header) = line {
            if i == 0 {
                match uniprot_fasta_header::uniprotkb(header.as_bytes()) {
                    Ok(_) => {}
                    Err(e) => {
                        println!("{}", e);
                    }
                }
            } else {
                match uniprot_fasta_header::uniprotkb_iso(header.as_bytes()) {
                    Ok(_) => {}
                    Err(e) => {
                        println!("{}", e);
                    }
                }
            }
        }
    }
}

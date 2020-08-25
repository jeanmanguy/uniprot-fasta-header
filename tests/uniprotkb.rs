use elapsed::measure_time;
use std::fs::File;
use std::io::{self, BufRead};

#[test]
fn e_coli() {
    let file = File::open("tests/assets/E_coli_headers.txt").unwrap();
    let mut counter: usize = 0;

    let (elapsed, _) = measure_time(|| {
        for line in io::BufReader::new(file).lines() {
            if let Ok(header) = line {
                match uniprot_fasta_header::uniprotkb(header.as_bytes()) {
                    Ok(_) => {}
                    Err(e) => {
                        println!("{}", e);
                    }
                }
                counter += 1;
            }
        }
    });
    println!("parsed {} headers in = {}", counter, elapsed);
}

use crate::{Database, ProteinExistence};
use nom::{
    branch::alt,
    bytes::complete::{tag, take, take_until, take_while1, take_while_m_n},
    character::{is_alphanumeric, is_digit},
    combinator::opt,
    error::ErrorKind,
    sequence::{pair, preceded, separated_pair, tuple},
    IResult, Slice,
};
use once_cell::sync::OnceCell;
use regex::bytes::Regex;

// Tag pipe
pub fn pipe(input: &[u8]) -> IResult<&[u8], &[u8]> {
    tag("|")(input)
}

// Tag chevron
pub fn chevron(input: &[u8]) -> IResult<&[u8], &[u8]> {
    tag(">")(input)
}

// Tag space
// one or more space ot tab
pub fn space(input: &[u8]) -> IResult<&[u8], &[u8]> {
    take_while1(|c| c == b' ' || c == b'\t')(input)
}

// UniProtKB database
// sp or tr
pub fn db(input: &[u8]) -> IResult<&[u8], Database> {
    let out: IResult<&[u8], &[u8]> = alt((tag("sp"), tag("tr")))(input);

    match out {
        Ok((rest, database)) => match database {
            b"sp" => Ok((rest, Database::SwissProt)),
            b"tr" => Ok((rest, Database::TrEMBL)),
            _ => unreachable!(),
        },
        Err(e) => Err(e),
    }
}

// UniProt unique ID
// https://www.uniprot.org/help/accession%5Fnumbers
pub fn unique_id(input: &[u8]) -> IResult<&[u8], &[u8]> {
    static RE: OnceCell<Regex> = OnceCell::new();
    let re = RE.get_or_init(|| {
        Regex::new("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}").unwrap()
    });

    if let Some(c) = re.captures(&input) {
        let v: &[u8] = c
            .iter()
            .filter(|el| el.is_some())
            .map(|el| el.unwrap())
            .map(|m| input.slice(m.start()..m.end()))
            .collect::<Vec<&[u8]>>()
            .first()
            .unwrap();

        let offset = { v.as_ptr() as usize + v.len() - input.as_ptr() as usize };

        Ok((input.slice(offset..), v))
    } else {
        let res: IResult<_, _> = Err(nom::Err::Error(nom::error::make_error(
            input,
            ErrorKind::RegexpCapture,
        )));
        res
    }
}
// UniProt entry name
// https://www.uniprot.org/help/entry%5Fname
pub fn entry_name(input: &[u8]) -> IResult<&[u8], Vec<u8>> {
    let out: IResult<&[u8], (&[u8], &[u8], &[u8])> = tuple((
        take_while_m_n(2, 12, is_alphanumeric), // mnemonic protein
        tag("_"),
        take_while_m_n(2, 5, is_alphanumeric), // mnemonic species
    ))(input);

    match out {
        Ok((rest, entry)) => Ok((rest, [entry.0, entry.1, entry.2].concat())),
        Err(e) => Err(e),
    }
}

// Protein name : everything until we reach OS=
pub fn until_os(input: &[u8]) -> IResult<&[u8], &[u8]> {
    take_until(" OS=")(input)
}

// Organism name
// Not simply binomial species name (not valid for viruses and addional strain information)
pub fn os_until_ox(input: &[u8]) -> IResult<&[u8], &[u8]> {
    preceded(tag("OS="), take_until(" OX="))(input)
}

// The word Fragment between parenthesis
// Useless for now
// pub fn fragment(input: &[u8]) -> IResult<&[u8], &[u8]> {
//     tag("(Fragment)")(input)
// }

// NCBI taxonomy ID
// https://www.uniprot.org/help/taxonomic%5Fidentifier
pub fn organism_id(input: &[u8]) -> IResult<&[u8], &[u8]> {
    preceded(tag("OX="), take_while_m_n(1, 7, is_digit))(input)
}

// Anything after GN until the end or until PE=
// Really anything, there are some truly weird gene names on UniProtKB
#[inline]
fn gn(input: &[u8]) -> IResult<&[u8], &[u8]> {
    preceded(
        tag("GN="),
        alt((
            take_until(" PE="),
            take_while1(|c| is_alphanumeric(c) || c == b'-' || c == b'_' || c == b'.'),
        )),
    )(input)
}

// Parse Gene Name
// Anything after GN until the end or until PE=. Really anything, there are some truly weird gene names on UniProtKB. Also some headers don't even have gene names.
pub fn optional_gene_name(input: &[u8]) -> IResult<&[u8], (Option<&[u8]>, Option<&[u8]>)> {
    pair(opt(gn), opt(space))(input)
}

// Protein existence evidence
// 1 to 5
pub fn evidence(input: &[u8]) -> IResult<&[u8], ProteinExistence> {
    let out: IResult<&[u8], &[u8]> = preceded(tag("PE="), take_while1(is_digit))(input);

    match out {
        Ok((rest, database)) => match database {
            b"1" => Ok((rest, ProteinExistence::ExperimentalEvidenceProtein)),
            b"2" => Ok((rest, ProteinExistence::ExperimentalEvidenceTranscript)),
            b"3" => Ok((rest, ProteinExistence::InferredHomology)),
            b"4" => Ok((rest, ProteinExistence::Predicted)),
            b"5" => Ok((rest, ProteinExistence::Uncertain)),
            _ => unreachable!(),
        },
        Err(e) => Err(e),
    }
}

// Version
pub fn version(input: &[u8]) -> IResult<&[u8], &[u8]> {
    preceded(tag("SV="), take(1usize))(input)
}

// Isoform ID : acccession number + isoform number
pub fn iso_id(input: &[u8]) -> IResult<&[u8], (&[u8], &[u8])> {
    separated_pair(unique_id, tag("-"), take_while1(is_digit))(input)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use test_case::test_case;

    /* ---------------------------------- pipe ---------------------------------- */
    #[test]
    fn nom_pipe() {
        assert_eq!(pipe(b"|xxx").unwrap(), (&b"xxx"[..], &b"|"[..]))
    }

    /* --------------------------------- chevron -------------------------------- */
    #[test]
    fn nom_chevron() {
        assert_eq!(chevron(b">xxx").unwrap(), (&b"xxx"[..], &b">"[..]))
    }

    /* ---------------------------------- space --------------------------------- */

    #[test_case(b" ", b" " ; "1")]
    #[test_case(b"  ", b"  " ; "2")]
    #[test_case(b"\t", b"\t" ; "tab")]
    fn spaces(input: &[u8], expected: &[u8]) {
        let (_, parsed) = space(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    /* ----------------------------------- db ----------------------------------- */

    #[test_case(b"tr",  Database::TrEMBL ; "TrEMBL")]
    #[test_case(b"sp",  Database::SwissProt ; "Swissprot")]
    fn uniprot_db(input: &[u8], expected: Database) {
        let (_, parsed) = db(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    #[test]
    fn incorrect_uniprot_db() {
        let err = db(b"xx").unwrap_err();
        assert_eq!(
            nom::Err::Error((&[120u8, 120u8][..], nom::error::ErrorKind::Tag)),
            err
        );
    }

    /* -------------------------------- unique_id ------------------------------- */
    #[test_case(b"A0A023GPI8", b"A0A023GPI8"; "LECA_CANBL")]
    #[test_case(b"P12345", b"P12345"; "AATM_RABIT")]
    #[test_case(b"A2BC19", b"A2BC19"; "A2BC19_HELPX")]
    fn uniprot_id(input: &[u8], expected: &[u8]) {
        let (_, parsed) = unique_id(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    #[test]
    fn incorrect_uniprot_id() {
        assert_eq!(
            unique_id(b"xx").unwrap_err(),
            nom::Err::Error((&[120u8, 120u8][..], nom::error::ErrorKind::RegexpCapture))
        );

        assert_eq!(
            unique_id(b"!P0266").unwrap_err(),
            nom::Err::Error((
                &[33u8, 80u8, 48u8, 50u8, 54u8, 54u8][..],
                nom::error::ErrorKind::RegexpCapture
            ))
        );
    }

    /* ------------------------------- entry_name ------------------------------- */

    #[test_case(b"ACN2_ACAGO", b"ACN2_ACAGO"; "ACN2_ACAGO")]
    fn correct_entry_name(input: &[u8], expected: &[u8]) {
        let (_, parsed) = entry_name(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    #[test]
    fn incorrect_entry_name() {
        let empty: &[u8] = &[];

        // missing _
        assert_eq!(
            entry_name(b"xx").unwrap_err(),
            nom::Err::Error((empty, nom::error::ErrorKind::Tag))
        );

        // first part too long
        assert_eq!(
            entry_name(b"XXXXXXXXXXXXXXXXXXXXXXXXXXX").unwrap_err(),
            nom::Err::Error((
                &[
                    88u8, 88u8, 88u8, 88u8, 88u8, 88u8, 88u8, 88u8, 88u8, 88u8, 88u8, 88u8, 88u8,
                    88u8, 88u8
                ][..],
                nom::error::ErrorKind::Tag
            ))
        );

        // first part too short
        assert_eq!(
            entry_name(b"X_X").unwrap_err(),
            nom::Err::Error((&[88u8, 95u8, 88u8][..], nom::error::ErrorKind::TakeWhileMN))
        );

        // absent first part
        assert_eq!(
            entry_name(b"_X").unwrap_err(),
            nom::Err::Error((&[95u8, 88u8][..], nom::error::ErrorKind::TakeWhileMN))
        );

        // absent second part
        assert_eq!(
            entry_name(b"XXXXX_").unwrap_err(),
            nom::Err::Error((empty, nom::error::ErrorKind::TakeWhileMN))
        );

        // wrong separator
        assert_eq!(
            entry_name(b"XXXXX-XXXXX").unwrap_err(),
            nom::Err::Error((
                &[45u8, 88u8, 88u8, 88u8, 88u8, 88u8][..],
                nom::error::ErrorKind::Tag
            )),
        );
    }

    /* -------------------------------- until_os -------------------------------- */

    #[test_case(b"Acanthoscurrin-2 (Fragment) OS=Acanthoscurria", b"Acanthoscurrin-2 (Fragment)"; "ACN2_ACAGO")]
    fn nom_until_os(input: &[u8], expected: &[u8]) {
        let (_, parsed) = until_os(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    /* -------------------------------- fragment -------------------------------- */

    // #[test]
    // fn nom_fragment() {
    //     assert_eq!(
    //         fragment(b"(Fragment) ").unwrap(),
    //         (&b" "[..], &b"(Fragment)"[..])
    //     )
    // }

    /* ------------------------------- os_until_ox ------------------------------ */

    #[test_case(b"OS=Bacillus phage SPP1 OX=10724", b"Bacillus phage SPP1"; "TERS_BPSPP Isoform G1P*")]
    fn nom_os_until_ox(input: &[u8], expected: &[u8]) {
        let (_, parsed) = os_until_ox(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    /* ------------------------------- organism_id ------------------------------ */

    #[test_case(b"OX=10724", b"10724"; "Bacillus phage SPP1")]
    fn nom_organism_id(input: &[u8], expected: &[u8]) {
        let (_, parsed) = organism_id(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    /* ----------------------------------- gn ----------------------------------- */

    #[test_case(b"GN=acantho2 PE=1 SV=1", b"acantho2"; "acantho2")]
    #[test_case(b"GN=SA85-1.3 PE=2 SV=1", b"SA85-1.3"; "P18271")]
    #[test_case(b"GN=> PE=4 SV=1", b">"; "Q63813")]
    #[test_case(b"GN=CSN3 PE=1 SV=1", b"CSN3"; "CSN3")]
    #[test_case(b"GN=0 beta-2 globin PE=3 SV=1", b"0 beta-2 globin"; "Q62670")]
    #[test_case(b"GN=orf304 = ymf42 PE=4 SV=1", b"orf304 = ymf42"; "Q35688")]
    #[test_case(b"GN=YWHAB", b"YWHAB"; "Q4R572-2")]
    fn gene_name(input: &[u8], expected: &[u8]) {
        let (_, parsed) = gn(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    /* ------------------------------ opt_gene_name ----------------------------- */

    #[test_case(b"GN=acantho2 PE=1 SV=1", (Some(b"acantho2"), Some(b" ")); "acantho2")]
    #[test_case(b"PE=2 SV=1", (None, None); "no gene name")]
    fn opt_gene_name(input: &[u8], expected: (Option<&[u8]>, Option<&[u8]>)) {
        let (_, parsed) = optional_gene_name(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    /* -------------------------------- evidence -------------------------------- */

    #[test_case(b"PE=1",  ProteinExistence::ExperimentalEvidenceProtein ; "ExperimentalEvidenceProtein")]
    #[test_case(b"PE=2",  ProteinExistence::ExperimentalEvidenceTranscript ; "ExperimentalEvidenceTranscript")]
    #[test_case(b"PE=3",  ProteinExistence::InferredHomology ; "InferredHomology")]
    #[test_case(b"PE=4",  ProteinExistence::Predicted ; "Predicted")]
    #[test_case(b"PE=5",  ProteinExistence::Uncertain ; "Uncertain")]
    fn existence(input: &[u8], expected: ProteinExistence) {
        let (_, parsed) = evidence(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    #[test]
    fn incorrect_existence() {
        assert_eq!(
            evidence(b"xx").unwrap_err(),
            nom::Err::Error((&[120u8, 120u8][..], nom::error::ErrorKind::Tag))
        );
    }

    /* --------------------------------- version -------------------------------- */

    #[test_case(b"SV=1", b"1"; "acantho2")]
    fn versions(input: &[u8], expected: &[u8]) {
        let (_, parsed) = version(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }

    /* --------------------------------- iso_id --------------------------------- */

    #[test_case(b"P54307-2|TERS_BPSPP Isoform G1P* of Terminase small subunit OS=Bacillus phage SPP1 OX=10724 GN=1", (b"P54307", b"2"); "TERS_BPSPP Isoform G1P*")]
    fn isoform_id(input: &[u8], expected: (&[u8], &[u8])) {
        let (_, parsed) = iso_id(input).unwrap();
        pretty_assertions::assert_eq!(parsed, expected);
    }
}

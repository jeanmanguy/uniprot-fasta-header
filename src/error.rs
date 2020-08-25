use displaydoc::Display;
use nom::error::{ErrorKind, ParseError};
use thiserror::Error;

/// # Errors
#[derive(PartialEq, Display, Error, Debug)]
#[non_exhaustive]
pub enum UniProtHeaderError {
    /// Failed to parse `{0}` :
    ParsingError(String, String),
    /// Incomplete
    Incomplete,
}

impl ParseError<&[u8]> for UniProtHeaderError {
    fn from_error_kind(input: &[u8], kind: ErrorKind) -> Self {
        UniProtHeaderError::ParsingError(
            String::from_utf8_lossy(input).to_string(),
            kind.description().to_string(),
        )
    }
    fn append(_input: &[u8], _kind: ErrorKind, other: Self) -> Self {
        other
    }
}

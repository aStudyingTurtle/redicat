use anyhow::Result;
use csv;
use grep_cli::stdout;
use gzp::{deflate::Gzip, BgzfSyncReader, Compression, ZBuilder};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;
use termcolor::ColorChoice;

/// Build a CSV reader for optional file/stdin sources.
pub fn get_reader<P: AsRef<Path>>(
    path: &Option<P>,
    has_headers: bool,
    bgzipped: bool,
) -> Result<csv::Reader<Box<dyn Read>>> {
    let raw_reader: Box<dyn Read> = match path {
        Some(path) if path.as_ref().to_str().unwrap() != "-" => {
            let reader = BufReader::new(File::open(path)?);
            if bgzipped {
                Box::new(BgzfSyncReader::new(reader))
            } else {
                Box::new(reader)
            }
        }
        _ => {
            let reader = io::stdin();
            if bgzipped {
                Box::new(BgzfSyncReader::new(reader))
            } else {
                Box::new(reader)
            }
        }
    };

    Ok(csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(has_headers)
        .from_reader(raw_reader))
}

/// Build a CSV writer targeting a file or stdout with optional BGZF compression.
pub fn get_writer<P: AsRef<Path>>(
    path: &Option<P>,
    gzipped: bool,
    write_headers: bool,
    threads: usize,
    compression_level: u32,
) -> Result<csv::Writer<Box<dyn Write>>> {
    let raw_writer: Box<dyn Write> = match path {
        Some(path) if path.as_ref().to_str().unwrap() != "-" => {
            let writer = BufWriter::new(File::create(path)?);
            if gzipped {
                Box::new(
                    ZBuilder::<Gzip, _>::new()
                        .num_threads(threads)
                        .compression_level(Compression::new(compression_level))
                        .from_writer(writer),
                )
            } else {
                Box::new(writer)
            }
        }
        _ => {
            let writer = stdout(ColorChoice::Never);
            if gzipped {
                Box::new(
                    ZBuilder::<Gzip, _>::new()
                        .num_threads(threads)
                        .compression_level(Compression::new(compression_level))
                        .from_writer(writer),
                )
            } else {
                Box::new(writer)
            }
        }
    };

    Ok(csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(write_headers)
        .from_writer(raw_writer))
}

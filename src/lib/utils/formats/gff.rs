use anyhow::bail;
use flate2::read::MultiGzDecoder;
use noodles_gff as gff;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

pub fn open<P>(src: P) -> anyhow::Result<gff::Reader<Box<dyn BufRead>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::open(path);
    match path.extension().and_then(|x| x.to_str()) {
        Some("gz") => {
            let reader = file.map(MultiGzDecoder::new).map(BufReader::new)?;
            Ok(gff::Reader::new(Box::new(reader)))
        }
        Some("gff") | Some("gff3") => {
            let reader = file.map(BufReader::new)?;
            Ok(gff::Reader::new(Box::new(reader)))
        }
        _ => bail!("Unknown extension for GFF file."),
    }
}
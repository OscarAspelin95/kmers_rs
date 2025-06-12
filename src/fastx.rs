use bio::io::{fasta, fastq};
use flate2::read::MultiGzDecoder;
use std::path::PathBuf;
use std::{fs::File, io::BufReader};

pub enum FileType {
    Plain,
    Gzip,
}
pub enum FastxFile {
    Fastq(FileType),
    Fasta(FileType),
}

fn file_is_gzip(fastx: &PathBuf) -> bool {
    return fastx.extension().unwrap() == "gz";
}

fn file_is_fastq(fastx: &PathBuf) -> Option<FastxFile> {
    let fastq_endings: [&str; 2] = [".fastq.gz", ".fq.gz"];

    let is_gzip = file_is_gzip(&fastx);

    let valid_fastq_ending = fastq_endings
        .iter()
        .any(|ending| fastx.to_str().unwrap().ends_with(ending));

    match (valid_fastq_ending, is_gzip) {
        (true, true) => return Some(FastxFile::Fastq(FileType::Gzip)),
        (true, false) => return Some(FastxFile::Fasta(FileType::Plain)),
        _ => None,
    }
}

fn file_is_fasta(fastx: &PathBuf) -> Option<FastxFile> {
    let fasta_endings: [&str; 8] = [
        ".fasta",
        ".fna",
        ".fsa",
        ".fa",
        ".fasta.gz",
        ".fna.gz",
        ".fsa.gz",
        ".fa.gz",
    ];

    let is_gzip = file_is_gzip(&fastx);

    let valid_fasta_ending = fasta_endings
        .iter()
        .any(|ending| fastx.to_str().unwrap().ends_with(ending));

    match (valid_fasta_ending, is_gzip) {
        (true, true) => return Some(FastxFile::Fasta(FileType::Gzip)),
        (true, false) => return Some(FastxFile::Fasta(FileType::Plain)),
        _ => None,
    }
}

pub fn fastx_file_type(fastx: PathBuf) -> FastxFile {
    if let Some(fastxfile) = file_is_fastq(&fastx) {
        return fastxfile;
    }

    if let Some(fastxfile) = file_is_fasta(&fastx) {
        return fastxfile;
    }

    panic!("Invalid file format/extension {:?}", fastx);
}

pub fn fastq_plain_reader(fastq: PathBuf) -> fastq::Reader<BufReader<File>> {
    let reader = fastq::Reader::from_file(fastq).unwrap();
    return reader;
}

pub fn fasta_plain_reader(fasta: PathBuf) -> fasta::Reader<BufReader<File>> {
    let reader = fasta::Reader::from_file(fasta).unwrap();
    return reader;
}

pub fn fastq_gz_reader(fastq: PathBuf) -> fastq::Reader<BufReader<MultiGzDecoder<File>>> {
    let f = File::open(fastq).unwrap();
    let reader = fastq::Reader::from_bufread(BufReader::new(MultiGzDecoder::new(f)));
    return reader;
}

pub fn fasta_gz_reader(fasta: PathBuf) -> fasta::Reader<BufReader<MultiGzDecoder<File>>> {
    let f = File::open(fasta).unwrap();
    let reader = fasta::Reader::from_bufread(BufReader::new(MultiGzDecoder::new(f)));

    return reader;
}

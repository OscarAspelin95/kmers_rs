#![feature(portable_simd)]

mod fastx;
mod kmers;
mod kmers_simd;
mod kmers_utils;
mod kmers_vanilla;

use clap::Parser;

use dashmap::{DashMap, DashSet};
use fastx::{fasta_plain_reader, fastq_gz_reader};
use kmers::generate_kmers;
use log::info;
use rayon::prelude::*;
use simple_logger::SimpleLogger;
use std::{collections::HashSet, path::PathBuf};

#[derive(Parser, Debug)]
#[command(version = "0.0.1", about = "Placeholder.", long_about = "Placeholder.")]
struct Args {
    #[arg(long)]
    fastq: PathBuf,

    #[arg(long)]
    fasta: PathBuf,

    #[arg(short, long, default_value_t = 21, value_parser = clap::value_parser!(u8).range(7..=31))]
    kmer_size: u8,

    #[arg(short, long, default_value_t = 1)]
    ds_factor: u64,
}

fn main() {
    SimpleLogger::new().init().unwrap();

    let args = Args::parse();

    let kmer_size: usize = args.kmer_size as usize;
    let ds_factor: u64 = args.ds_factor;

    // BufReader for gzipped fastq.
    let fasta_reader = fasta_plain_reader(args.fasta);
    let fastq_reader = fastq_gz_reader(args.fastq);

    info!("Generating fasta kmers...");
    let fasta_dashmap: DashMap<String, HashSet<u64>> = DashMap::with_capacity(100_000);
    fasta_reader.records().par_bridge().for_each(|record| {
        let rec = record.unwrap();

        if rec.seq().len() < 200 {
            return;
        }
        let nt_string = rec.seq();
        let fasta_kmers = generate_kmers(nt_string, kmer_size, ds_factor);

        fasta_dashmap.insert(rec.id().to_owned(), fasta_kmers);
    });

    // Process reads in parallel.
    info!("Generating fastq kmers...");
    let fastq_dashmap: DashMap<String, HashSet<u64>> = DashMap::with_capacity(100_000);
    fastq_reader.records().par_bridge().for_each(|record| {
        let rec = record.unwrap();

        if rec.seq().len() < 1000 {
            return;
        }
        let nt_string = rec.seq();
        let fastq_kmers = generate_kmers(nt_string, kmer_size, ds_factor);

        fastq_dashmap.insert(rec.id().to_owned(), fastq_kmers);
    });
}


# Introduction
🚧 Work in progress for generating (MinFracHash) kmers from very long sequences in fastq and fasta files.

It uses a rather interesting approach, inspired by [sylph](https://github.com/bluenote-1577/sylph), for generating kmers.
  1. Uses rayon for parallel processing of reads/contigs.
  2. Each read/contig is chopped into 64 equal size chunks.
  3. These 64 chunks are processed in parallel using Rust [SIMD](https://doc.rust-lang.org/std/simd/type.u64x64.html).
  4. Since sequences are not always divisible by 64, the residual chunk is processed with a non-SIMD implementation.

# Note
Required Rust Nightly since it uses the simd module.

# TODO
* Currently generates kmers, but does not store them anywhere.
* Want some kind of containment ID between fasta contigs and fastq reads.
* ...

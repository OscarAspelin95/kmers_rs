use crate::kmers_simd::simd_u64_64_encoding;
use crate::kmers_utils::chunk_pos;
use crate::kmers_vanilla::kmerize;
use std::collections::HashSet;

#[inline]
fn generate_simd_kmers(
    nt_string: &[u8],
    simd_pos_chunks: &[(usize, usize)],
    kmer_size: usize,
    ds_factor: u64,
) -> HashSet<u64> {
    let simd_nt_chunks: Vec<&[u8]> = simd_pos_chunks
        .iter()
        .map(|range| {
            let (start, end) = *range;

            return &nt_string[start..end];
        })
        .collect();

    assert_eq!(simd_nt_chunks.len(), 64);

    let kmer_hash_set = simd_u64_64_encoding(kmer_size, &simd_nt_chunks, ds_factor);

    return kmer_hash_set;
}
#[inline]
fn generate_residual_kmers(
    nt_string: &[u8],
    last_pos_chunk: &(usize, usize),
    kmer_size: usize,
    ds_factor: u64,
) -> HashSet<u64> {
    // There is probably a better way of extracting two elements for a length two vector.
    let (last_start, last_end) = *last_pos_chunk;

    // This is the last part of nt_string that we cannot use SIMD for.
    let last_nt_string = &nt_string[last_start..last_end];

    let num_residual_kmers = kmerize(kmer_size, ds_factor, last_nt_string);

    return num_residual_kmers;
}
/// At the moment, we only allow splitting a nt string into 65 pieces.
/// This gives us 64 chunks of equal length to kmerize using SIMD and
/// a final chunk of arbitrary length, which we kmerize separately.
///
/// A neat thing would be to calculate for each nt_string what SIMD lane count we want to use.
/// Currently, there is a minimum read length requirement, dictated by the kmer size because
/// we slice the nt string into 64 chunks. If the nt string is too short, we cannot mathematically
/// create 64 chunks with lengths > kmer_size.
#[inline]
pub fn generate_kmers(nt_string: &[u8], kmer_size: usize, ds_factor: u64) -> HashSet<u64> {
    let chunks = chunk_pos(nt_string.len(), 64, kmer_size);
    assert_eq!(chunks.len(), 65);

    // ---- SIMD kmers.
    let simd_pos_chunks = &chunks[0..64];
    assert_eq!(simd_pos_chunks.len(), 64);

    // ---- SIMD kmer hashes.
    let mut simd_kmers = generate_simd_kmers(nt_string, simd_pos_chunks, kmer_size, ds_factor);

    // ---- Residual kmers.
    let last_pos_chunk: &(usize, usize) = chunks.last().unwrap();
    let residual_kmers = generate_residual_kmers(nt_string, last_pos_chunk, kmer_size, ds_factor);

    // ---- Extend our simd kmers with the residual ones.
    simd_kmers.extend(residual_kmers);
    return simd_kmers;
}

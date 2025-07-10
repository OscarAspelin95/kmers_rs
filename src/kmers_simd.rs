use crate::kmers_utils::BYTE_TO_SEQ;
use std::{
    collections::HashSet,
    ops::{BitAnd, BitOr, Sub},
    simd::{Simd, cmp::SimdOrd, u64x64},
    u64,
};

/// Should be benchmarked against the
/// original sylph implementation but seems to work.
#[inline]
fn mm_hash64_simd(kmer: Simd<u64, 64>) -> Simd<u64, 64> {
    let mut key = kmer;

    let all_ones = Simd::splat(u64::MAX); // equivalent to !0
    key = all_ones ^ (key + (key << Simd::splat(21))); // Bitwise NOT via XOR
    key ^= key >> Simd::splat(24);
    key = (key + (key << Simd::splat(3))) + (key << Simd::splat(8)); // key * 265
    key ^= key >> Simd::splat(14);
    key = (key + (key << Simd::splat(2))) + (key << Simd::splat(4)); // key * 21
    key ^= key >> Simd::splat(28);
    key += key << Simd::splat(31);

    return key;
}

/// SIMD implementation of the MinFracHash algorithm.
/// Currently only supports 64 lane SIMD
#[inline]
pub fn simd_u64_64_encoding(
    kmer_size: usize,
    simd_nt_chunks: Vec<&[u8]>,
    ds_factor: u64,
) -> HashSet<u64> {
    // Initialize empty forward and reverse kmer.
    let mut kmer_fwd: std::simd::Simd<u64, 64> = u64x64::splat(0);
    let mut kmer_reverse = u64x64::splat(0);

    // Used for reverse kmer.
    let nbits = kmer_size << 1;
    let mask: Simd<u64, 64> = u64x64::splat((1 << nbits) - 1);
    let shift = ((kmer_size - 1) * 2) as u64;

    // NOTE - this is a temp solution to get the chunk size.
    let chunk_len = simd_nt_chunks.first().unwrap().len();

    // Store the kmers in a hash set. For capacity, e.g., if the
    // chunk_len is 200bp, we have 64 nt sequences, each of length 200bp.
    // The number of kmers we'll generate is 64 * (chunk_len + kmer_size - 1).
    // We check both forward and reverse strand, but only keep one of them (smallest one).
    let hash_capacity = 64 * (chunk_len + kmer_size - 1);
    let mut hash_set = HashSet::with_capacity(hash_capacity);

    for i in 0..chunk_len {
        // Houston, we have a problem. The SIMD implementation is very fast. In fact, it is
        // so fast that this closure becomes the bottleneck.
        //
        // One better way could be to pre-compute all SIMDs beforehand.
        //
        // We have some type casting u8 to usize to u64 that we want to either
        // * Do in parallel
        // * Avoid all together (not sure this is possible).
        let nt_array: Vec<u64> = simd_nt_chunks
            .iter()
            .map(|chunk| {
                return BYTE_TO_SEQ[chunk[i] as usize] as u64;
            })
            .collect();

        // Forward.
        let nt_simd: Simd<u64, 64> = u64x64::from_slice(nt_array.as_slice());

        // Reverse.
        let nt_rev_simd = u64x64::splat(3).sub(nt_simd);

        // Forward kmer
        kmer_fwd = ((kmer_fwd << 2).bitor(nt_simd)).bitand(mask);

        // Reverse kmer
        kmer_reverse = kmer_reverse >> 2;
        kmer_reverse = kmer_reverse.bitor(nt_rev_simd << shift);

        // We have a valid kmer.
        // HOWEVER, we need to be able to handle N.
        if i >= kmer_size - 1 {
            // Get the lexicographically smallest kmer out of fws/rev.
            let smallest = kmer_fwd.simd_min(kmer_reverse);

            // Apply hashing to each SIMD lane.
            let simd_hash = mm_hash64_simd(smallest);

            simd_hash.to_array().iter().for_each(|hash| {
                if *hash <= u64::MAX / ds_factor {
                    hash_set.insert(*hash);
                }
            });
        }
    }
    return hash_set;
}

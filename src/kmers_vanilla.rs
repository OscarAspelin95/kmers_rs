use crate::kmers_utils::BYTE_TO_SEQ;

#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    let mut key = kmer;
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

pub fn kmerize(k: usize, ds_factor: u64, nt_string: &[u8]) -> usize {
    let shift = ((k - 1) * 2) as u64;

    let mut last_kmer_forward: u64 = 0;
    let mut last_kmer_reverse: u64 = 0;

    let mut num_residual_kmers: usize = 0;

    for i in 0..nt_string.len() {
        let nt = BYTE_TO_SEQ[nt_string[i] as usize] as u64;
        let nt_rev = 3 - nt;

        last_kmer_forward = (last_kmer_forward << 2) | nt;

        last_kmer_reverse = last_kmer_reverse >> 2;
        last_kmer_reverse = last_kmer_reverse | (nt_rev << shift);

        if i >= k - 1 {
            let canonical;

            if last_kmer_forward < last_kmer_reverse {
                canonical = last_kmer_forward;
            } else {
                canonical = last_kmer_reverse;
            }

            let canonical_hash = mm_hash64(canonical);

            if canonical_hash <= (u64::MAX / ds_factor) {
                num_residual_kmers += 1;
            }
        }
    }

    return num_residual_kmers;
}

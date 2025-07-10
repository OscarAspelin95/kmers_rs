/// Lookup table for ASCII nucleotides (a,A,c,C,g,G,t,T,u,U)
/// We map nucleotides accordingly:
///     * a/A => 0
///     * c/C => 1
///     * g/G => 2
///     * t/T/u/U => 3
///
/// This enables us to do easy reverse complement through 3-nt.
pub const BYTE_TO_SEQ: [u8; 256] = [
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

#[inline]
fn decode(byte: u64) -> u8 {
    match byte {
        0 => return b'A',
        1 => return b'C',
        2 => return b'G',
        3 => return b'T',
        _ => panic!("Invalid nucleotide."),
    };
}

/// Print a u64 encoded nucleotide with some bit manipulation.
pub fn print_nt_string(kmer: u64, k: usize) {
    let mut result = String::with_capacity(k);
    for i in 0..k {
        // Shift to extract the 2 bits corresponding to the current nucleotide
        let shift = 2 * (k - i - 1);
        let bits = (kmer >> shift) & 0b11;
        let base = match bits {
            0b00 => 'A',
            0b01 => 'C',
            0b10 => 'G',
            0b11 => 'T',
            _ => unreachable!(),
        };
        result.push(base);
    }
    println!("{}", result);
}

/// In order to fully use the SIMD implementation, we need to take a nt string and chop it up into
/// equal length substrings that we can kmerize in parallel. However, in order to make sure we
/// get all kmers, we need to:
/// * Make sure substrings overlap by kmer_size-1 nucleotides.
/// * If we cannot evently divide the sequence into equal length parts, we need to
///     generate num_chunks + 1 chunks, where the last chunk can have arbitrary length > 0
///     and all previous chunks have the same length.
///
/// Example, we have a nt string we want to chop up into 64 parts with kmer_size-1 overlap
/// to extract every single kmer. If the nt string length cannot be evenly chopped up into
/// 64 subsequences of equal length, we need 65 chunks. 64 with equal length and the remaining
/// chunk to complete the length of nt string. Then, we can run the 64 equal length chunks
/// with SIMD kmerization and the last chunk with a non-SIMD kmerization.
#[inline]
pub fn chunk_pos(seq_len: usize, num_chunks: usize, kmer_size: usize) -> Vec<(usize, usize)> {
    let chunk_len = (seq_len + ((num_chunks - 1) * (kmer_size - 1))) / num_chunks;

    if kmer_size > chunk_len {
        panic!(
            "Kmer size {} cannot be larger than chunk length {}. Sequence is too short to chunk given the kmer size.",
            num_chunks, kmer_size
        );
    }

    let mut start = 0;
    let mut end = chunk_len;

    let mut chunk_index: Vec<(usize, usize)> = Vec::with_capacity(num_chunks + 1);

    for _ in 0..num_chunks {
        chunk_index.push((start, end));

        start += chunk_len - kmer_size + 1;
        end = chunk_len + start;
    }

    chunk_index.push((start, seq_len));
    assert_eq!(chunk_index.len(), num_chunks + 1);

    return chunk_index;
}

use kmers_rs::kmers::generate_kmers;

use rstest::rstest;

#[rstest]
#[case(7, 1, &[b'G'; 1000])]
#[case(21, 1, &[b'G'; 867])]
#[case(31, 1, &[b'G'; 999])]
#[case(31, 1, &[b'T'; 786])]
#[case(31, 1, &[b'A'; 843])]
#[case(31, 1, &[b'C'; 1933])]

/// Test that we generate the correct number of kmers.
fn test_kmers(#[case] kmer_size: usize, #[case] ds_factor: u64, #[case] nt_string: &[u8]) {
    let num_kmers = generate_kmers(nt_string, kmer_size, ds_factor);

    let expected_num_hashes = nt_string.len() - kmer_size + 1;

    assert_eq!(num_kmers, expected_num_hashes);
}

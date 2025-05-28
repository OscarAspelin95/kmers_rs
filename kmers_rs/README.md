# README #

Repository for generating kmers efficiently with Rust. NOTE - requires Rust nightly due to the use of the std::simd module.



### How do I get set up? ###

* Requires Rust (check out rustup) Nightly:
    * Run "rustup toolchain install nightly"
    * Run "rustup override set nightly" in this directory.
* Run "cargo build" for generating a debug binary.
* Run "cargo build --release" for generating the release binary.

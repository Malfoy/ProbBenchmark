# ProbBenchmark

Rust command-line benchmarks for k-mer indexing and querying on FASTA files using multiple Bloom filter implementations.

The project currently provides eight binaries:

- `fastbloom`: benchmark based on [`fastbloom`](https://crates.io/crates/fastbloom)
- `classic_bloom`: benchmark based on [`bloom-filters`](https://crates.io/crates/bloom-filters)
- `roaring_bloom`: benchmark based on [`roaring-bloom-filter`](https://crates.io/crates/roaring-bloom-filter)
- `bloom_filter_rs`: benchmark based on [`bloom-filter-rs`](https://crates.io/crates/bloom-filter-rs)
- `bloomfx`: benchmark based on [`bloomfx`](https://crates.io/crates/bloomfx)
- `bloom_rs`: benchmark based on [`bloom_rs`](https://crates.io/crates/bloom_rs)
- `generic_bloom`: benchmark based on [`generic-bloom`](https://crates.io/crates/generic-bloom)
- `friendly-d-kmers`: generator for Friendly-D k-mers from an indexed FASTA file

Each benchmark binary:

- reads one FASTA file to index
- extracts all valid DNA k-mers (`A`, `C`, `G`, `T`)
- inserts indexed k-mers into a Bloom filter
- reads a second FASTA file to query
- checks query k-mers against the Bloom filter
- reports timing, CPU usage, memory usage, and k-mer counts

## Build

The repository is configured for aggressive release optimization:

- `opt-level = 3`
- `lto = "fat"`
- `codegen-units = 1`
- `panic = "abort"` for release
- `strip = "symbols"`
- default `-C target-cpu=native`

Because `target-cpu=native` is enabled in [.cargo/config.toml](/home/nadine/Code/ProbBenchmark/.cargo/config.toml), release binaries are optimized for the machine that builds them and may be less portable to older CPUs.

Build all release binaries:

```bash
cargo build -r
```

Build a specific binary:

```bash
cargo build -r --bin fastbloom
cargo build -r --bin classic_bloom
cargo build -r --bin roaring_bloom
cargo build -r --bin bloom_filter_rs
cargo build -r --bin bloomfx
cargo build -r --bin bloom_rs
cargo build -r --bin generic_bloom
cargo build -r --bin friendly-d-kmers
```

Run tests:

```bash
cargo test
```

## Output Binaries

After `cargo build -r`, the binaries are available in:

```bash
target/release/fastbloom
target/release/classic_bloom
target/release/roaring_bloom
target/release/bloom_filter_rs
target/release/bloomfx
target/release/bloom_rs
target/release/generic_bloom
target/release/friendly-d-kmers
```

## Benchmark Usage

The seven benchmark binaries use the same CLI.

Required arguments:

- `--index-fasta <FASTA>`
- `--query-fasta <FASTA>`
- `-k, --kmer-size <INT>`

Optional Bloom filter parameters:

- `--bloom-bits <INT>`
- `--bloom-bytes <INT>`
- `--hashes <INT>`
- `--false-positive-rate <FLOAT>` default: `0.01`

Optional runtime parameters:

- `--threads <INT>`
- `--batch-bases <INT>` default: `33554432`

Example with explicit Bloom size and hash count:

```bash
./target/release/fastbloom \
  --index-fasta ref.fa \
  --query-fasta query.fa \
  --kmer-size 31 \
  --bloom-bits 1000000000 \
  --hashes 7 \
  --threads 16
```

Example with automatic sizing:

```bash
./target/release/classic_bloom \
  --index-fasta ref.fa \
  --query-fasta query.fa \
  --kmer-size 31 \
  --false-positive-rate 0.001 \
  --threads 16
```

Example with the roaring implementation:

```bash
./target/release/roaring_bloom \
  --index-fasta ref.fa \
  --query-fasta query.fa \
  --kmer-size 31 \
  --bloom-bytes 134217728 \
  --hashes 6
```

Example with `bloom-filter-rs`:

```bash
./target/release/bloom_filter_rs \
  --index-fasta ref.fa \
  --query-fasta query.fa \
  --kmer-size 31 \
  --bloom-bits 1000000000 \
  --hashes 7
```

Example with `bloomfx`:

```bash
./target/release/bloomfx \
  --index-fasta ref.fa \
  --query-fasta query.fa \
  --kmer-size 31 \
  --bloom-bits 1000000000 \
  --hashes 7
```

Example with `bloom_rs`:

```bash
./target/release/bloom_rs \
  --index-fasta ref.fa \
  --query-fasta query.fa \
  --kmer-size 31 \
  --bloom-bits 1000000000 \
  --hashes 7
```

Example with `generic-bloom`:

```bash
./target/release/generic_bloom \
  --index-fasta ref.fa \
  --query-fasta query.fa \
  --kmer-size 31 \
  --bloom-bits 1000000000 \
  --hashes 7
```

## Reported Metrics

Each run prints tab-separated key/value pairs:

- `index_wall_time_s`
- `query_wall_time_s`
- `index_cpu_time_s`
- `query_cpu_time_s`
- `indexed_kmers`
- `queried_kmers`
- `query_positive_kmers`
- `bloom_bits`
- `bloom_bytes`
- `bloom_hashes`
- `max_ram_bytes`
- `max_ram_mib`
- `threads`
- `precount_pass`

Notes:

- `max_ram_bytes` and `max_ram_mib` are the process peak RSS for the whole run, not separate per-phase peaks.
- `precount_pass=true` means the program performed an extra pass over the indexed FASTA to estimate the number of k-mers because `--bloom-bits` or `--hashes` was not fully specified.
- only k-mers composed entirely of `A`, `C`, `G`, and `T` are indexed or queried

## Implementation Notes

### `fastbloom`

- Uses `fastbloom::AtomicBloomFilter`
- Supports concurrent insertions directly
- Best option in this repository for parallel indexing throughput

### `classic_bloom`

- Uses the classic Bloom filter from `bloom-filters`
- Indexing uses one shared Bloom filter protected by a mutex
- Avoids duplicate filter copies, but mutation is synchronized

### `roaring_bloom`

- Uses `roaring-bloom-filter::StableBloomFilter`
- Indexing also uses one shared Bloom filter protected by a mutex
- This crate requires sized values for insertion and lookup, so k-mers are hashed to `u64` keys before insertion and query

### `bloom_filter_rs`

- Uses `bloom-filter-rs::BloomFilter`
- Supports direct manual configuration of bit-array size and hash count
- Operates directly on k-mer byte slices

### `bloomfx`

- Uses `bloomfx::BloomFilter`
- Supports direct manual configuration of bit count and hash count
- This benchmark hashes k-mers to `u64` keys before insertion and query

### `bloom_rs`

- Uses `bloom_rs::BloomFilter`
- Supports direct manual configuration of bit count and hash count
- This benchmark hashes k-mers to `u64` keys before insertion and query

### `generic_bloom`

- Uses `generic_bloom::SimpleBloomFilter`

## Friendly-D Generator

The repository also includes the Friendly-D k-mer generator previously kept in `FDP/`.

Build only that binary:

```bash
cargo build -r --bin friendly-d-kmers
```

Run it:

```bash
./target/release/friendly-d-kmers <input.fasta> <k> <n> <d> [--output output.fa] [--seed value]
```

Arguments:

- `input.fasta`: input FASTA file to index
- `k`: k-mer size, from `1` to `31`
- `n`: number of Friendly-D k-mers to generate
- `d`: mutation distance from the start or end, from `1` to `k`
- `--output output.fa`: optional output FASTA path; if omitted, FASTA is written to stdout
- `--seed value`: optional random seed
- Supports direct manual configuration of counter count and hash count
- Uses a plain bitset-backed binary Bloom filter configuration

## Crate Origins

### Bloom filter crates

- `fastbloom`
  - Crate: <https://crates.io/crates/fastbloom>
  - Documentation: <https://docs.rs/fastbloom/0.17.0>
  - Repository: <https://github.com/tomtomwombat/fastbloom>

- `bloom-filters`
  - Crate: <https://crates.io/crates/bloom-filters>
  - Documentation: <https://docs.rs/bloom-filters/0.1.2>
  - Repository: <https://github.com/nervosnetwork/bloom-filters>
  - Notes: Rust port of BoomFilters, as stated by the crate metadata

- `roaring-bloom-filter`
  - Crate: <https://crates.io/crates/roaring-bloom-filter>
  - Documentation: <https://docs.rs/roaring-bloom-filter/0.2.0>
  - Repository: <https://github.com/oliverdding/roaring-bloom-filter-rs>
  - License note: this crate is published under `AGPL-3.0`

- `bloom-filter-rs`
  - Crate: <https://crates.io/crates/bloom-filter-rs>
  - Documentation: <https://docs.rs/bloom-filter-rs/0.1.0>
  - Repository: <https://github.com/sagalasan/bloom-filter>

- `bloomfx`
  - Crate: <https://crates.io/crates/bloomfx>
  - Documentation: <https://docs.rs/bloomfx/0.1.1>
  - Repository: <https://github.com/jdockerty/bloomfx>

- `bloom_rs`
  - Crate: <https://crates.io/crates/bloom_rs>
  - Documentation: <https://docs.rs/bloom_rs/0.1.0/bloom_rs>
  - Repository: <https://github.com/DavidCai1993/bloom.rs>

- `generic-bloom`
  - Crate: <https://crates.io/crates/generic-bloom>
  - Documentation: <https://docs.rs/generic-bloom/0.1.0>
  - Repository: <https://github.com/goose121/generic-bloom-rs>
  - License note: this crate is published under `AGPL-3.0-or-later`

### Supporting crates used by this project

- `clap`
  - Crate: <https://crates.io/crates/clap>
  - Repository: <https://github.com/clap-rs/clap-rs>
  - Purpose: CLI argument parsing

- `rayon`
  - Crate: <https://crates.io/crates/rayon>
  - Repository: <https://github.com/rayon-rs/rayon>
  - Purpose: parallel batch and sequence processing

- `anyhow`
  - Crate: <https://crates.io/crates/anyhow>
  - Repository: <https://github.com/dtolnay/anyhow>
  - Purpose: ergonomic error handling

- `libc`
  - Crate: <https://crates.io/crates/libc>
  - Repository: <https://github.com/rust-lang/libc>
  - Purpose: `getrusage` access for CPU and peak RAM reporting

## Project Layout

```text
src/lib.rs                  Shared FASTA parsing, k-mer traversal, metrics, reporting
src/bin/fastbloom.rs        fastbloom benchmark binary
src/bin/classic_bloom.rs    bloom-filters benchmark binary
src/bin/roaring_bloom.rs    roaring-bloom-filter benchmark binary
src/bin/bloom_filter_rs.rs  bloom-filter-rs benchmark binary
src/bin/bloomfx.rs          bloomfx benchmark binary
src/bin/bloom_rs.rs         bloom_rs benchmark binary
src/bin/generic_bloom.rs    generic-bloom benchmark binary
.cargo/config.toml          default target-cpu=native configuration
```

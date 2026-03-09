# Friendly-D k-mers

This folder is now integrated into the root `probbenchmark` Cargo package.
The generator source lives in `FDP/src/main.rs`, but it is built from the repository root.

Build:

```bash
cargo build --release --bin friendly-d-kmers
```

Run:

```bash
./target/release/friendly-d-kmers <input.fasta> <k> <n> <d> [--output output.fa] [--seed value]
```

The binary indexes k-mers from an input FASTA file and then randomly generates `N` unique Friendly-D k-mers.

## Concept

Given:

- a k-mer length `k`
- a distance `D`
- a set of indexed source k-mers extracted from a FASTA file

a **Friendly-D k-mer** is a k-mer that:

1. is obtained from an indexed source k-mer,
2. differs by exactly one nucleotide,
3. changes the nucleotide at distance `D` from the start or at distance `D` from the end,
4. is **not** present in the original indexed k-mer set.

Examples with `k = 5`:

- `D = 1` means mutate the first or last nucleotide.
- `D = 2` means mutate the second or the fourth nucleotide.

If `D` points to the same position from both sides, only that single position is used.



## Usage

```bash
cargo run --release --bin friendly-d-kmers -- <input.fasta> <k> <n> <d> [--output output.fa] [--seed value]
```

Arguments:

- `input.fasta`: input FASTA file to index
- `k`: k-mer size, from `1` to `31`
- `n`: number of Friendly-D k-mers to generate
- `d`: mutation distance from the start or end, from `1` to `k`
- `--output output.fa`: optional output FASTA path; if omitted, FASTA is written to stdout
- `--seed value`: optional random seed

Example:

```bash
cargo run --release --bin friendly-d-kmers -- reads.fa 31 1000 1 --output friendly.fa
```

## Output format

The output is FASTA, with one record per generated Friendly-D k-mer:

```fasta
>friendly_1 d=1 source=ACGT
GCGT
>friendly_2 d=1 source=GTAC
CTAC
```

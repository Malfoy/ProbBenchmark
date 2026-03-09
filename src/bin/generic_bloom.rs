use anyhow::{Context, Result};
use bitvec::prelude::{BitBox, Lsb0};
use clap::Parser;
use generic_bloom::{BloomFilter, SimpleBloomFilter};
use probbenchmark::{
    BenchmarkReport, CommonArgs, QueryStats, configure_threads, count_fasta_kmers,
    cpu_delta_seconds, optimal_num_hashes, print_report, resolve_bloom_bits, scan_fasta_batches,
    usage_snapshot, validate_common_args, visit_sequence_kmers,
};
use rayon::prelude::*;
use std::collections::hash_map::RandomState;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::Instant;

type GenericFilter = SimpleBloomFilter<BitBox<usize, Lsb0>, RandomState, Box<[RandomState]>>;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Indexes k-mers from one FASTA into a Bloom filter and queries k-mers from another FASTA with generic-bloom."
)]
struct Cli {
    #[command(flatten)]
    common: CommonArgs,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    validate_common_args(&cli.common)?;
    configure_threads(cli.common.threads)?;

    let report = run_benchmark(&cli.common)?;
    print_report(&report);
    Ok(())
}

fn run_benchmark(args: &CommonArgs) -> Result<BenchmarkReport> {
    let threads = rayon::current_num_threads();
    let needs_precount =
        args.hashes.is_none() || (args.bloom_bits.is_none() && args.bloom_bytes.is_none());

    let index_usage_start = usage_snapshot()?;
    let index_wall_start = Instant::now();

    let indexed_kmers = if needs_precount {
        count_fasta_kmers(&args.index_fasta, args.kmer_size, args.batch_bases)?
    } else {
        0
    };

    let bloom_bits = resolve_bloom_bits(args, indexed_kmers)?;
    let bloom_hashes = match args.hashes {
        Some(hashes) => hashes,
        None => {
            let expected_items = usize::try_from(indexed_kmers.max(1))
                .context("indexed k-mer count does not fit in usize on this platform")?;
            optimal_num_hashes(bloom_bits, expected_items)
        }
    };
    let mut bloom = GenericFilter::new(bloom_hashes as usize, bloom_bits);

    let indexed_kmers = if needs_precount {
        let inserted_kmers = index_fasta(
            &args.index_fasta,
            args.kmer_size,
            args.batch_bases,
            &mut bloom,
        )?;
        debug_assert_eq!(inserted_kmers, indexed_kmers);
        indexed_kmers
    } else {
        index_fasta(
            &args.index_fasta,
            args.kmer_size,
            args.batch_bases,
            &mut bloom,
        )?
    };

    let index_usage_end = usage_snapshot()?;
    let index_wall_time_s = index_wall_start.elapsed().as_secs_f64();
    let index_cpu_time_s = cpu_delta_seconds(index_usage_start, index_usage_end);

    let query_usage_start = usage_snapshot()?;
    let query_wall_start = Instant::now();
    let query_stats = query_fasta(&args.query_fasta, args.kmer_size, args.batch_bases, &bloom)?;
    let query_usage_end = usage_snapshot()?;

    Ok(BenchmarkReport {
        index_wall_time_s,
        query_wall_time_s: query_wall_start.elapsed().as_secs_f64(),
        index_cpu_time_s,
        query_cpu_time_s: cpu_delta_seconds(query_usage_start, query_usage_end),
        max_ram_bytes: query_usage_end.max_rss_bytes,
        indexed_kmers,
        queried_kmers: query_stats.total_kmers,
        query_positive_kmers: query_stats.positive_kmers,
        bloom_bits,
        bloom_hashes,
        threads,
        precount_pass: needs_precount,
    })
}

fn index_fasta(
    path: &Path,
    k: usize,
    batch_bases: usize,
    global_filter: &mut GenericFilter,
) -> Result<u64> {
    let mut total = 0_u64;
    let shared_filter = Arc::new(Mutex::new(std::mem::replace(
        global_filter,
        GenericFilter::new(1, 1),
    )));

    scan_fasta_batches(path, batch_bases, |batch| {
        let batch_total: u64 = batch
            .par_iter()
            .map(|seq| {
                let mut local_kmers = Vec::new();
                let indexed_kmers = visit_sequence_kmers(seq, k, |kmer| {
                    local_kmers.push(kmer.to_vec());
                });

                let mut filter = shared_filter.lock().expect("generic-bloom mutex poisoned");
                for kmer in &local_kmers {
                    filter.insert(kmer);
                }
                indexed_kmers
            })
            .sum();

        total += batch_total;
        Ok(())
    })?;

    *global_filter = Arc::into_inner(shared_filter)
        .expect("generic-bloom filter still has multiple owners")
        .into_inner()
        .expect("generic-bloom mutex poisoned");

    Ok(total)
}

fn query_fasta(
    path: &Path,
    k: usize,
    batch_bases: usize,
    filter: &GenericFilter,
) -> Result<QueryStats> {
    let mut totals = QueryStats::default();
    scan_fasta_batches(path, batch_bases, |batch| {
        let batch_stats = batch
            .par_iter()
            .map(|seq| {
                let mut stats = QueryStats::default();
                visit_sequence_kmers(seq, k, |kmer| {
                    stats.total_kmers += 1;
                    if filter.contains(&kmer) {
                        stats.positive_kmers += 1;
                    }
                });
                stats
            })
            .reduce(QueryStats::default, QueryStats::combine);
        totals = totals.combine(batch_stats);
        Ok(())
    })?;
    Ok(totals)
}

#[cfg(test)]
mod tests {
    use super::GenericFilter;
    use generic_bloom::BloomFilter;
    use probbenchmark::visit_sequence_kmers;

    #[test]
    fn generic_bloom_reports_inserted_hits() {
        let mut filter = GenericFilter::new(3, 1024);
        visit_sequence_kmers(b"ACGT", 3, |kmer| {
            filter.insert(&kmer.to_vec());
        });

        let mut positives = 0_u64;
        let mut total = 0_u64;
        visit_sequence_kmers(b"ACGT", 3, |kmer| {
            total += 1;
            if filter.contains(&kmer.to_vec()) {
                positives += 1;
            }
        });

        assert_eq!(total, 2);
        assert_eq!(positives, 2);
    }
}

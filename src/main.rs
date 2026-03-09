use anyhow::{Context, Result};
use clap::Parser;
use fastbloom::AtomicBloomFilter;
use probbenchmark::{
    BenchmarkReport, CommonArgs, QueryStats, configure_threads, count_fasta_kmers,
    cpu_delta_seconds, optimal_num_hashes, print_report, resolve_bloom_bits, scan_fasta_batches,
    usage_snapshot, validate_common_args, visit_sequence_kmers,
};
use rayon::prelude::*;
use std::time::Instant;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Indexes k-mers from one FASTA into a Bloom filter and queries k-mers from another FASTA with fastbloom."
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
    let needs_precount = args.hashes.is_none() || (args.bloom_bits.is_none() && args.bloom_bytes.is_none());

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
    let bloom = AtomicBloomFilter::with_num_bits(bloom_bits).hashes(bloom_hashes);

    let indexed_kmers = if needs_precount {
        let inserted_kmers = index_fasta(&args.index_fasta, args.kmer_size, args.batch_bases, &bloom)?;
        debug_assert_eq!(inserted_kmers, indexed_kmers);
        indexed_kmers
    } else {
        index_fasta(&args.index_fasta, args.kmer_size, args.batch_bases, &bloom)?
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

fn index_fasta(path: &std::path::Path, k: usize, batch_bases: usize, filter: &AtomicBloomFilter) -> Result<u64> {
    let mut total = 0_u64;
    scan_fasta_batches(path, batch_bases, |batch| {
        let batch_total: u64 = batch
            .par_iter()
            .map(|seq| visit_sequence_kmers(seq, k, |kmer| {
                filter.insert(kmer);
            }))
            .sum();
        total += batch_total;
        Ok(())
    })?;
    Ok(total)
}

fn query_fasta(
    path: &std::path::Path,
    k: usize,
    batch_bases: usize,
    filter: &AtomicBloomFilter,
) -> Result<QueryStats> {
    let mut totals = QueryStats::default();
    scan_fasta_batches(path, batch_bases, |batch| {
        let batch_stats = batch
            .par_iter()
            .map(|seq| {
                let mut stats = QueryStats::default();
                visit_sequence_kmers(seq, k, |kmer| {
                    stats.total_kmers += 1;
                    if filter.contains(kmer) {
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
    use super::query_fasta;
    use fastbloom::AtomicBloomFilter;
    use probbenchmark::{QueryStats, visit_sequence_kmers};
    use std::path::Path;

    #[test]
    fn visit_and_query_counts_hits() {
        let filter = AtomicBloomFilter::with_num_bits(1024).hashes(3);
        visit_sequence_kmers(b"ACGT", 3, |kmer| {
            filter.insert(kmer);
        });

        let mut stats = QueryStats::default();
        visit_sequence_kmers(b"ACGT", 3, |kmer| {
            stats.total_kmers += 1;
            if filter.contains(kmer) {
                stats.positive_kmers += 1;
            }
        });

        assert_eq!(stats.total_kmers, 2);
        assert_eq!(stats.positive_kmers, 2);
    }

    #[test]
    fn query_on_empty_fasta_is_zero() {
        let filter = AtomicBloomFilter::with_num_bits(128).hashes(2);
        let result = query_fasta(Path::new("/dev/null"), 3, 1024, &filter);
        assert!(result.is_ok());
        let stats = result.unwrap_or_default();
        assert_eq!(stats.total_kmers, 0);
        assert_eq!(stats.positive_kmers, 0);
    }
}

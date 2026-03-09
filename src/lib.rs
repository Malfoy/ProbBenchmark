use anyhow::{Context, Result, bail};
use clap::Args;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::collections::hash_map::DefaultHasher;
use std::f64::consts::LN_2;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

pub const DEFAULT_BATCH_BASES: usize = 32 * 1024 * 1024;

#[derive(Debug, Clone, Args)]
pub struct CommonArgs {
    #[arg(long, value_name = "FASTA")]
    pub index_fasta: PathBuf,

    #[arg(long, value_name = "FASTA")]
    pub query_fasta: PathBuf,

    #[arg(short = 'k', long)]
    pub kmer_size: usize,

    #[arg(long, conflicts_with = "bloom_bytes")]
    pub bloom_bits: Option<usize>,

    #[arg(long, conflicts_with = "bloom_bits")]
    pub bloom_bytes: Option<usize>,

    #[arg(long)]
    pub hashes: Option<u32>,

    #[arg(long, default_value_t = 0.01)]
    pub false_positive_rate: f64,

    #[arg(long)]
    pub threads: Option<usize>,

    #[arg(long, default_value_t = DEFAULT_BATCH_BASES)]
    pub batch_bases: usize,
}

#[derive(Debug)]
pub struct BenchmarkReport {
    pub index_wall_time_s: f64,
    pub query_wall_time_s: f64,
    pub index_cpu_time_s: f64,
    pub query_cpu_time_s: f64,
    pub max_ram_bytes: u64,
    pub indexed_kmers: u64,
    pub queried_kmers: u64,
    pub query_positive_kmers: u64,
    pub bloom_bits: usize,
    pub bloom_hashes: u32,
    pub threads: usize,
    pub precount_pass: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct UsageSnapshot {
    pub user_seconds: f64,
    pub system_seconds: f64,
    pub max_rss_bytes: u64,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct QueryStats {
    pub total_kmers: u64,
    pub positive_kmers: u64,
}

impl QueryStats {
    pub fn combine(self, other: Self) -> Self {
        Self {
            total_kmers: self.total_kmers + other.total_kmers,
            positive_kmers: self.positive_kmers + other.positive_kmers,
        }
    }
}

pub fn validate_common_args(args: &CommonArgs) -> Result<()> {
    if args.kmer_size == 0 {
        bail!("k-mer size must be greater than 0");
    }
    if args.batch_bases == 0 {
        bail!("batch size must be greater than 0");
    }
    if let Some(bits) = args.bloom_bits
        && bits == 0
    {
        bail!("bloom_bits must be greater than 0");
    }
    if let Some(bytes) = args.bloom_bytes
        && bytes == 0
    {
        bail!("bloom_bytes must be greater than 0");
    }
    if let Some(hashes) = args.hashes
        && hashes == 0
    {
        bail!("hashes must be greater than 0");
    }
    if args.false_positive_rate <= 0.0 || args.false_positive_rate >= 1.0 {
        bail!("false_positive_rate must be between 0 and 1");
    }
    if let Some(threads) = args.threads
        && threads == 0
    {
        bail!("threads must be greater than 0");
    }
    Ok(())
}

pub fn configure_threads(threads: Option<usize>) -> Result<()> {
    if let Some(threads) = threads {
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .map_err(|err| anyhow::anyhow!("failed to configure Rayon thread pool: {err}"))?;
    }
    Ok(())
}

pub fn count_fasta_kmers(path: &Path, k: usize, batch_bases: usize) -> Result<u64> {
    let mut total = 0_u64;
    scan_fasta_batches(path, batch_bases, |batch| {
        let batch_total: u64 = batch
            .par_iter()
            .map(|seq| count_sequence_kmers(seq, k))
            .sum();
        total += batch_total;
        Ok(())
    })?;
    Ok(total)
}

pub fn count_sequence_kmers(seq: &[u8], k: usize) -> u64 {
    let mut valid_bases = 0_usize;
    let mut kmers = 0_u64;

    for &base in seq {
        if is_valid_base(base) {
            valid_bases += 1;
            if valid_bases >= k {
                kmers += 1;
            }
        } else {
            valid_bases = 0;
        }
    }

    kmers
}

pub fn visit_sequence_kmers<F>(seq: &[u8], k: usize, mut visitor: F) -> u64
where
    F: FnMut(&[u8]),
{
    let mut valid_bases = 0_usize;
    let mut kmers = 0_u64;

    for (end, &base) in seq.iter().enumerate() {
        if is_valid_base(base) {
            valid_bases += 1;
            if valid_bases >= k {
                let start = end + 1 - k;
                visitor(&seq[start..=end]);
                kmers += 1;
            }
        } else {
            valid_bases = 0;
        }
    }

    kmers
}

pub fn scan_fasta_batches<F>(path: &Path, batch_bases: usize, mut on_batch: F) -> Result<()>
where
    F: FnMut(Vec<Vec<u8>>) -> Result<()>,
{
    let file = File::open(path)
        .with_context(|| format!("failed to open FASTA file {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut line = Vec::with_capacity(16 * 1024);
    let mut batch = Vec::new();
    let mut batch_base_count = 0_usize;
    let mut current_seq = Vec::new();
    let mut saw_header = false;

    loop {
        line.clear();
        let bytes_read = reader
            .read_until(b'\n', &mut line)
            .with_context(|| format!("failed while reading FASTA file {}", path.display()))?;
        if bytes_read == 0 {
            break;
        }

        trim_line_end(&mut line);
        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if saw_header {
                push_record(
                    &mut current_seq,
                    &mut batch,
                    &mut batch_base_count,
                    batch_bases,
                    &mut on_batch,
                )?;
            } else {
                saw_header = true;
            }
            continue;
        }

        if !saw_header {
            bail!(
                "FASTA file {} is malformed: sequence data found before the first header",
                path.display()
            );
        }

        append_sequence_line(&mut current_seq, &line);
    }

    if saw_header {
        push_record(
            &mut current_seq,
            &mut batch,
            &mut batch_base_count,
            batch_bases,
            &mut on_batch,
        )?;
    }

    if !batch.is_empty() {
        on_batch(batch)?;
    }

    Ok(())
}

pub fn resolve_bloom_bits(args: &CommonArgs, indexed_kmers: u64) -> Result<usize> {
    if let Some(bits) = args.bloom_bits {
        return Ok(bits);
    }

    if let Some(bytes) = args.bloom_bytes {
        return bytes
            .checked_mul(8)
            .context("bloom_bytes overflowed when converted to bits");
    }

    let expected_items = usize::try_from(indexed_kmers.max(1))
        .context("indexed k-mer count does not fit in usize on this platform")?;
    Ok(optimal_num_bits(expected_items, args.false_positive_rate))
}

pub fn optimal_num_bits(expected_items: usize, false_positive_rate: f64) -> usize {
    let expected_items = expected_items.max(1) as f64;
    let ln2_sq = LN_2 * LN_2;
    ((expected_items * false_positive_rate.ln().abs()) / ln2_sq).ceil() as usize
}

pub fn optimal_num_hashes(num_bits: usize, expected_items: usize) -> u32 {
    let expected_items = expected_items.max(1) as f64;
    let hashes = ((num_bits as f64 / expected_items) * LN_2).ceil() as u32;
    hashes.max(1)
}

pub fn hash_bytes_default(bytes: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    bytes.hash(&mut hasher);
    hasher.finish()
}

pub fn usage_snapshot() -> Result<UsageSnapshot> {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    // SAFETY: getrusage initializes the provided rusage struct on success.
    let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
    if rc != 0 {
        return Err(std::io::Error::last_os_error()).context("getrusage failed");
    }

    // SAFETY: the previous getrusage call succeeded.
    let usage = unsafe { usage.assume_init() };
    Ok(UsageSnapshot {
        user_seconds: timeval_to_seconds(usage.ru_utime),
        system_seconds: timeval_to_seconds(usage.ru_stime),
        max_rss_bytes: max_rss_to_bytes(usage.ru_maxrss),
    })
}

pub fn cpu_delta_seconds(start: UsageSnapshot, end: UsageSnapshot) -> f64 {
    (end.user_seconds + end.system_seconds) - (start.user_seconds + start.system_seconds)
}

pub fn print_report(report: &BenchmarkReport) {
    println!("index_wall_time_s\t{:.6}", report.index_wall_time_s);
    println!("query_wall_time_s\t{:.6}", report.query_wall_time_s);
    println!("index_cpu_time_s\t{:.6}", report.index_cpu_time_s);
    println!("query_cpu_time_s\t{:.6}", report.query_cpu_time_s);
    println!("indexed_kmers\t{}", report.indexed_kmers);
    println!("queried_kmers\t{}", report.queried_kmers);
    println!("query_positive_kmers\t{}", report.query_positive_kmers);
    println!("bloom_bits\t{}", report.bloom_bits);
    println!("bloom_bytes\t{}", report.bloom_bits.div_ceil(8));
    println!("bloom_hashes\t{}", report.bloom_hashes);
    println!("max_ram_bytes\t{}", report.max_ram_bytes);
    println!(
        "max_ram_mib\t{:.3}",
        report.max_ram_bytes as f64 / (1024.0 * 1024.0)
    );
    println!("threads\t{}", report.threads);
    println!("precount_pass\t{}", report.precount_pass);
}

fn push_record<F>(
    current_seq: &mut Vec<u8>,
    batch: &mut Vec<Vec<u8>>,
    batch_base_count: &mut usize,
    batch_bases: usize,
    on_batch: &mut F,
) -> Result<()>
where
    F: FnMut(Vec<Vec<u8>>) -> Result<()>,
{
    *batch_base_count += current_seq.len();
    batch.push(std::mem::take(current_seq));

    if *batch_base_count >= batch_bases {
        let ready = std::mem::take(batch);
        *batch_base_count = 0;
        on_batch(ready)?;
    }

    Ok(())
}

fn append_sequence_line(sequence: &mut Vec<u8>, line: &[u8]) {
    sequence.reserve(line.len());
    for &byte in line {
        if !byte.is_ascii_whitespace() {
            sequence.push(byte.to_ascii_uppercase());
        }
    }
}

fn trim_line_end(line: &mut Vec<u8>) {
    while let Some(&last) = line.last() {
        if last == b'\n' || last == b'\r' {
            line.pop();
        } else {
            break;
        }
    }
}

#[inline]
fn is_valid_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

fn timeval_to_seconds(timeval: libc::timeval) -> f64 {
    timeval.tv_sec as f64 + (timeval.tv_usec as f64 / 1_000_000.0)
}

#[cfg(target_os = "macos")]
fn max_rss_to_bytes(max_rss: libc::c_long) -> u64 {
    max_rss as u64
}

#[cfg(not(target_os = "macos"))]
fn max_rss_to_bytes(max_rss: libc::c_long) -> u64 {
    (max_rss as u64) * 1024
}

#[cfg(test)]
mod tests {
    use super::{count_sequence_kmers, visit_sequence_kmers};

    #[test]
    fn counts_only_valid_windows() {
        assert_eq!(count_sequence_kmers(b"ACGTAC", 3), 4);
        assert_eq!(count_sequence_kmers(b"ACNTAC", 3), 1);
        assert_eq!(count_sequence_kmers(b"AAAANAAA", 4), 1);
    }

    #[test]
    fn visits_expected_kmers() {
        let mut kmers = Vec::new();
        let total = visit_sequence_kmers(b"ACGTNACGT", 4, |kmer| kmers.push(kmer.to_vec()));

        assert_eq!(total, 2);
        assert_eq!(kmers, vec![b"ACGT".to_vec(), b"ACGT".to_vec()]);
    }
}

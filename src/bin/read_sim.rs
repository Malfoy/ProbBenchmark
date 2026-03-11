use anyhow::{Context, Result, bail};
use clap::Parser;
use probbenchmark::configure_threads;
use rand::rngs::SmallRng;
use rand::{Rng, RngCore, SeedableRng};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

const INPUT_BUFFER_BYTES: usize = 8 * 1024 * 1024;
const OUTPUT_BUFFER_BYTES: usize = 8 * 1024 * 1024;
const THREAD_FLUSH_BYTES: usize = 4 * 1024 * 1024;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Simulate FASTA reads from a reference with substitutions at a fixed error rate."
)]
struct Cli {
    #[arg(long, value_name = "FASTA")]
    reference_fasta: PathBuf,

    #[arg(long, value_name = "FLOAT")]
    error_rate: f64,

    #[arg(long, value_name = "INT")]
    length: usize,

    #[arg(long, value_name = "FLOAT")]
    depth: f64,

    #[arg(long)]
    threads: Option<usize>,

    #[arg(long)]
    seed: Option<u64>,

    #[arg(long, value_name = "FASTA")]
    output: Option<PathBuf>,
}

struct ReferenceIndex {
    fragments: Vec<Vec<u8>>,
    cumulative_starts: Vec<u64>,
    total_bases: u64,
    total_starts: u64,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    validate_args(&cli)?;
    configure_threads(cli.threads)?;

    let index = load_reference(&cli.reference_fasta, cli.length)?;
    if index.total_starts == 0 {
        bail!(
            "reference {} has no valid ACGT window of length {}",
            cli.reference_fasta.display(),
            cli.length
        );
    }

    let read_count = compute_read_count(index.total_bases, cli.length, cli.depth)?;
    simulate_reads(
        &index,
        cli.length,
        read_count,
        cli.error_rate,
        cli.seed,
        cli.output,
    )
}

fn validate_args(args: &Cli) -> Result<()> {
    if args.length == 0 {
        bail!("length must be greater than 0");
    }
    if args.depth <= 0.0 {
        bail!("depth must be greater than 0");
    }
    if !args.depth.is_finite() {
        bail!("depth must be finite");
    }
    if !(0.0..=1.0).contains(&args.error_rate) || !args.error_rate.is_finite() {
        bail!("error_rate must be a finite value between 0 and 1");
    }
    if let Some(threads) = args.threads
        && threads == 0
    {
        bail!("threads must be greater than 0");
    }
    Ok(())
}

fn load_reference(path: &PathBuf, read_len: usize) -> Result<ReferenceIndex> {
    let file =
        File::open(path).with_context(|| format!("failed to open reference {}", path.display()))?;
    let mut reader = BufReader::with_capacity(INPUT_BUFFER_BYTES, file);
    let mut line = Vec::with_capacity(16 * 1024);
    let mut current = Vec::new();
    let mut fragments = Vec::new();
    let mut saw_header = false;

    loop {
        line.clear();
        let bytes_read = reader
            .read_until(b'\n', &mut line)
            .with_context(|| format!("failed while reading reference {}", path.display()))?;
        if bytes_read == 0 {
            break;
        }

        trim_line_end(&mut line);
        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if !saw_header {
                saw_header = true;
            }
            if !current.is_empty() {
                fragments.push(std::mem::take(&mut current));
            }
            continue;
        }

        if !saw_header {
            bail!(
                "reference {} is malformed: sequence data found before first header",
                path.display()
            );
        }

        for &byte in &line {
            match normalize_base(byte) {
                Some(base) => current.push(base),
                None => {
                    if byte.is_ascii_whitespace() {
                        continue;
                    }
                    if !current.is_empty() {
                        fragments.push(std::mem::take(&mut current));
                    }
                }
            }
        }
    }

    if !current.is_empty() {
        fragments.push(current);
    }

    if fragments.is_empty() {
        bail!(
            "reference {} does not contain any valid A/C/G/T bases",
            path.display()
        );
    }

    let mut filtered = Vec::with_capacity(fragments.len());
    let mut cumulative_starts = Vec::with_capacity(fragments.len());
    let mut total_bases = 0_u64;
    let mut total_starts = 0_u64;

    for fragment in fragments {
        if fragment.len() < read_len {
            continue;
        }

        total_bases = total_bases
            .checked_add(fragment.len() as u64)
            .context("reference base count overflowed")?;
        let starts = (fragment.len() - read_len + 1) as u64;
        total_starts = total_starts
            .checked_add(starts)
            .context("reference start count overflowed")?;

        filtered.push(fragment);
        cumulative_starts.push(total_starts);
    }

    Ok(ReferenceIndex {
        fragments: filtered,
        cumulative_starts,
        total_bases,
        total_starts,
    })
}

fn compute_read_count(total_bases: u64, read_len: usize, depth: f64) -> Result<u64> {
    let reads = (total_bases as f64 * depth) / read_len as f64;
    if !reads.is_finite() || reads > u64::MAX as f64 {
        bail!("requested depth produces an invalid read count");
    }
    let read_count = reads.ceil() as u64;
    if read_count == 0 {
        bail!("requested depth produces zero reads");
    }
    Ok(read_count)
}

fn simulate_reads(
    index: &ReferenceIndex,
    read_len: usize,
    read_count: u64,
    error_rate: f64,
    seed: Option<u64>,
    output: Option<PathBuf>,
) -> Result<()> {
    let writer_target: Box<dyn Write + Send> = match output {
        Some(path) => Box::new(
            File::create(&path)
                .with_context(|| format!("failed to create output {}", path.display()))?,
        ),
        None => Box::new(io::stdout()),
    };
    let writer = Arc::new(Mutex::new(BufWriter::with_capacity(
        OUTPUT_BUFFER_BYTES,
        writer_target,
    )));

    let base_seed = seed.unwrap_or_else(rand::random::<u64>);
    let threshold = error_threshold(error_rate);
    let thread_count = rayon::current_num_threads().max(1) as u64;

    (0..thread_count)
        .into_par_iter()
        .try_for_each(|worker_id| -> Result<()> {
            let start = (read_count * worker_id) / thread_count;
            let end = (read_count * (worker_id + 1)) / thread_count;
            if start == end {
                return Ok(());
            }

            let mut rng = SmallRng::seed_from_u64(splitmix64(base_seed ^ worker_id));
            let mut read_buf = vec![0_u8; read_len];
            let mut out_buf = Vec::with_capacity(THREAD_FLUSH_BYTES + read_len + 64);

            for read_id in start..end {
                let (fragment, offset) = sample_start(index, &mut rng);
                read_buf.copy_from_slice(&fragment[offset..offset + read_len]);
                apply_substitution_errors(&mut read_buf, threshold, &mut rng);

                out_buf.extend_from_slice(b">read_");
                append_u64(&mut out_buf, read_id + 1);
                out_buf.push(b'\n');
                out_buf.extend_from_slice(&read_buf);
                out_buf.push(b'\n');

                if out_buf.len() >= THREAD_FLUSH_BYTES {
                    flush_thread_buffer(&writer, &mut out_buf)?;
                }
            }

            flush_thread_buffer(&writer, &mut out_buf)?;
            Ok(())
        })?;

    let mut guard = writer
        .lock()
        .map_err(|_| anyhow::anyhow!("writer mutex poisoned"))?;
    guard.flush().context("failed to flush output")?;
    Ok(())
}

fn sample_start<'a>(index: &'a ReferenceIndex, rng: &mut SmallRng) -> (&'a [u8], usize) {
    let ticket = rng.gen_range(0..index.total_starts);
    let fragment_idx = index
        .cumulative_starts
        .partition_point(|&value| value <= ticket);
    let previous = if fragment_idx == 0 {
        0
    } else {
        index.cumulative_starts[fragment_idx - 1]
    };
    let offset = (ticket - previous) as usize;
    (&index.fragments[fragment_idx], offset)
}

fn apply_substitution_errors(read: &mut [u8], threshold: Option<u64>, rng: &mut SmallRng) {
    let Some(threshold) = threshold else {
        return;
    };

    for base in read {
        if rng.next_u64() <= threshold {
            let choice = (rng.next_u32() % 3) as u8;
            *base = substitute_base(*base, choice);
        }
    }
}

fn substitute_base(base: u8, choice: u8) -> u8 {
    match base {
        b'A' => [b'C', b'G', b'T'][choice as usize],
        b'C' => [b'A', b'G', b'T'][choice as usize],
        b'G' => [b'A', b'C', b'T'][choice as usize],
        b'T' => [b'A', b'C', b'G'][choice as usize],
        _ => base,
    }
}

fn error_threshold(error_rate: f64) -> Option<u64> {
    if error_rate <= 0.0 {
        return None;
    }
    if error_rate >= 1.0 {
        return Some(u64::MAX);
    }
    Some((error_rate * u64::MAX as f64) as u64)
}

fn flush_thread_buffer(
    writer: &Arc<Mutex<BufWriter<Box<dyn Write + Send>>>>,
    buffer: &mut Vec<u8>,
) -> Result<()> {
    if buffer.is_empty() {
        return Ok(());
    }
    let mut guard = writer
        .lock()
        .map_err(|_| anyhow::anyhow!("writer mutex poisoned"))?;
    guard
        .write_all(buffer)
        .context("failed while writing simulated reads")?;
    buffer.clear();
    Ok(())
}

fn append_u64(buf: &mut Vec<u8>, value: u64) {
    let mut tmp = [0_u8; 20];
    let mut value = value;
    let mut idx = tmp.len();

    loop {
        idx -= 1;
        tmp[idx] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }

    buf.extend_from_slice(&tmp[idx..]);
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

fn normalize_base(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(b'A'),
        b'C' => Some(b'C'),
        b'G' => Some(b'G'),
        b'T' => Some(b'T'),
        _ => None,
    }
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

#[cfg(test)]
mod tests {
    use super::{append_u64, compute_read_count, substitute_base};

    #[test]
    fn read_count_uses_coverage_formula() {
        let reads = compute_read_count(10_000, 100, 20.0).expect("coverage should be valid");
        assert_eq!(reads, 2_000);
    }

    #[test]
    fn append_u64_writes_ascii_digits() {
        let mut buf = Vec::new();
        append_u64(&mut buf, 12_340_001);
        assert_eq!(buf, b"12340001");
    }

    #[test]
    fn substitutions_never_return_original_base() {
        for &base in b"ACGT" {
            for choice in 0..3 {
                assert_ne!(substitute_base(base, choice), base);
            }
        }
    }
}

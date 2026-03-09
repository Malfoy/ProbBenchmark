use ahash::AHashSet;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

type PackedKmer = u64;

#[derive(Clone, Debug, Eq, PartialEq)]
struct Config {
    input: PathBuf,
    k: usize,
    n: usize,
    d: usize,
    output: Option<PathBuf>,
    seed: Option<u64>,
}

#[derive(Clone, Copy, Debug)]
struct FriendlyPair {
    friendly: PackedKmer,
    source: PackedKmer,
}

struct SearchState {
    seen: AHashSet<PackedKmer>,
    pairs: Vec<FriendlyPair>,
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().collect();
    if let Some(program) = args.first() {
        if args.iter().skip(1).any(|arg| arg == "-h" || arg == "--help") {
            println!("{}", usage(program));
            return ExitCode::SUCCESS;
        }
    }

    match run(args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            let _ = writeln!(io::stderr().lock(), "error: {message}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: Vec<String>) -> Result<(), String> {
    let config = parse_args(args)?;
    validate_config(&config)?;

    let (source_set, source_kmers) =
        index_fasta(&config.input, config.k).map_err(|err| format!("failed to index FASTA: {err}"))?;
    if source_kmers.is_empty() {
        return Err("no valid kmers were found in the FASTA file".to_string());
    }

    let selected =
        sample_friendly_kmers_parallel(source_set, source_kmers, config.k, config.d, config.n, config.seed)?;
    write_output(config.output.as_deref(), &selected, config.k, config.d)
        .map_err(|err| format!("failed to write output: {err}"))?;

    Ok(())
}

fn parse_args<I>(args: I) -> Result<Config, String>
where
    I: IntoIterator<Item = String>,
{
    let mut args = args.into_iter();
    let program = args
        .next()
        .unwrap_or_else(|| "friendly-d-kmers".to_string());
    let values: Vec<String> = args.collect();

    if values.is_empty() || values.iter().any(|arg| arg == "-h" || arg == "--help") {
        return Err(usage(&program));
    }

    if values.len() < 4 {
        return Err(format!("missing required arguments\n\n{}", usage(&program)));
    }

    let input = PathBuf::from(&values[0]);
    let k = parse_usize(&values[1], "k")?;
    let n = parse_usize(&values[2], "n")?;
    let d = parse_usize(&values[3], "d")?;

    let mut output = None;
    let mut seed = None;
    let mut index = 4;

    while index < values.len() {
        match values[index].as_str() {
            "--output" => {
                index += 1;
                let value = values
                    .get(index)
                    .ok_or_else(|| "--output requires a file path".to_string())?;
                output = Some(PathBuf::from(value));
            }
            "--seed" => {
                index += 1;
                let value = values
                    .get(index)
                    .ok_or_else(|| "--seed requires an integer value".to_string())?;
                seed = Some(
                    value
                        .parse::<u64>()
                        .map_err(|_| format!("invalid seed value: {value}"))?,
                );
            }
            flag => {
                return Err(format!("unknown argument: {flag}\n\n{}", usage(&program)));
            }
        }
        index += 1;
    }

    Ok(Config {
        input,
        k,
        n,
        d,
        output,
        seed,
    })
}

fn parse_usize(value: &str, name: &str) -> Result<usize, String> {
    value
        .parse::<usize>()
        .map_err(|_| format!("invalid {name} value: {value}"))
}

fn usage(program: &str) -> String {
    format!(
        "Usage: {program} <input.fasta> <k> <n> <d> [--output output.fa] [--seed value]\n\
         \n\
         Generates N unique Friendly-D k-mers and writes them in FASTA format.\n\
         Each output sequence is absent from the original index and differs from\n\
         an indexed source k-mer at distance D from the start or the end.\n\
         \n\
         Constraints:\n\
         - 1 <= k <= 31\n\
         - 1 <= d <= k"
    )
}

fn validate_config(config: &Config) -> Result<(), String> {
    if config.k == 0 || config.k > 31 {
        return Err("k must be between 1 and 31 to fit the packed u64 2-bit encoding".to_string());
    }
    if config.n == 0 {
        return Err("n must be greater than zero".to_string());
    }
    if config.d == 0 || config.d > config.k {
        return Err("d must be between 1 and k".to_string());
    }
    Ok(())
}

fn index_fasta(path: &Path, k: usize) -> io::Result<(AHashSet<PackedKmer>, Vec<PackedKmer>)> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let mut line = String::new();
    let mut seen = AHashSet::new();
    let mut ordered = Vec::new();

    let mask = mask_for_k(k);
    let mut rolling = 0u64;
    let mut valid_len = 0usize;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        if line.starts_with('>') {
            rolling = 0;
            valid_len = 0;
            continue;
        }

        for &byte in line.as_bytes() {
            match encode_base(byte) {
                Some(bits) => {
                    rolling = ((rolling << 2) | u64::from(bits)) & mask;
                    valid_len = valid_len.saturating_add(1).min(k);
                    if valid_len == k && seen.insert(rolling) {
                        ordered.push(rolling);
                    }
                }
                None => {
                    if byte != b'\n' && byte != b'\r' {
                        rolling = 0;
                        valid_len = 0;
                    }
                }
            }
        }
    }

    Ok((seen, ordered))
}

fn sample_friendly_kmers_parallel(
    source_set: AHashSet<PackedKmer>,
    source_kmers: Vec<PackedKmer>,
    k: usize,
    d: usize,
    n: usize,
    seed: Option<u64>,
) -> Result<Vec<FriendlyPair>, String> {
    let left = d - 1;
    let right = k - d;
    let source_set = Arc::new(source_set);
    let source_kmers = Arc::new(source_kmers);
    let attempts = Arc::new(AtomicUsize::new(0));
    let max_attempts = max_attempts(source_kmers.len(), n);
    let state = Arc::new(Mutex::new(SearchState {
        seen: AHashSet::with_capacity(n.saturating_mul(2)),
        pairs: Vec::with_capacity(n),
    }));

    (0..rayon::current_num_threads())
        .into_par_iter()
        .for_each(|worker_id| {
            let mut rng = make_worker_rng(seed, worker_id);

            loop {
                {
                    let guard = state.lock().expect("search state mutex poisoned");
                    if guard.pairs.len() >= n {
                        break;
                    }
                }

                let attempt = attempts.fetch_add(1, Ordering::Relaxed);
                if attempt >= max_attempts {
                    break;
                }

                let source = source_kmers[rng.gen_range(0..source_kmers.len())];
                let position = pick_position(left, right, &mut rng);
                let friendly = mutate_at_position(source, k, position, &mut rng);

                if source_set.contains(&friendly) {
                    continue;
                }

                let mut guard = state.lock().expect("search state mutex poisoned");
                if guard.pairs.len() >= n {
                    break;
                }
                if guard.seen.insert(friendly) {
                    guard.pairs.push(FriendlyPair { friendly, source });
                }
            }
        });

    let mut guard = state.lock().expect("search state mutex poisoned");
    if guard.pairs.len() < n {
        return Err(format!(
            "unable to find {} unique friendly k-mers after {} random attempts; found {}",
            n,
            max_attempts,
            guard.pairs.len()
        ));
    }

    guard.pairs.sort_unstable_by_key(|pair| pair.friendly);
    Ok(guard.pairs.drain(..n).collect())
}

fn write_output(path: Option<&Path>, selected: &[FriendlyPair], k: usize, d: usize) -> io::Result<()> {
    match path {
        Some(path) => {
            let file = File::create(path)?;
            let writer = BufWriter::with_capacity(1024 * 1024, file);
            write_fasta(writer, selected, k, d)
        }
        None => {
            let stdout = io::stdout();
            let writer = BufWriter::with_capacity(1024 * 1024, stdout.lock());
            write_fasta(writer, selected, k, d)
        }
    }
}

fn write_fasta<W: Write>(mut writer: W, selected: &[FriendlyPair], k: usize, d: usize) -> io::Result<()> {
    for (index, pair) in selected.iter().enumerate() {
        let friendly = decode_kmer(pair.friendly, k);
        let source = decode_kmer(pair.source, k);
        writeln!(
            writer,
            ">friendly_{} d={} source={}",
            index + 1,
            d,
            source
        )?;
        writeln!(writer, "{friendly}")?;
    }
    writer.flush()
}

fn max_attempts(source_count: usize, n: usize) -> usize {
    source_count
        .saturating_mul(128)
        .max(n.saturating_mul(200_000))
        .max(1_000_000)
}

fn make_worker_rng(seed: Option<u64>, worker_id: usize) -> SmallRng {
    match seed {
        Some(seed) => SmallRng::seed_from_u64(
            seed ^ (0x9E37_79B9_7F4A_7C15_u64.wrapping_mul(worker_id as u64 + 1)),
        ),
        None => SmallRng::from_entropy(),
    }
}

fn pick_position<R: Rng>(left: usize, right: usize, rng: &mut R) -> usize {
    if left == right || rng.gen_bool(0.5) {
        left
    } else {
        right
    }
}

fn mutate_at_position<R: Rng>(source: PackedKmer, k: usize, position: usize, rng: &mut R) -> PackedKmer {
    let shift = 2 * (k - 1 - position);
    let original_base = ((source >> shift) & 0b11) as u8;
    let mut alt_base = rng.gen_range(0..3) as u8;
    if alt_base >= original_base {
        alt_base += 1;
    }
    (source & !(0b11u64 << shift)) | (u64::from(alt_base) << shift)
}

fn mask_for_k(k: usize) -> u64 {
    (1u64 << (2 * k)) - 1
}

fn encode_base(byte: u8) -> Option<u8> {
    match byte.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn decode_kmer(mut packed: PackedKmer, k: usize) -> String {
    let mut bytes = vec![b'A'; k];
    for index in (0..k).rev() {
        bytes[index] = decode_base((packed & 0b11) as u8);
        packed >>= 2;
    }
    String::from_utf8(bytes).expect("kmer contains only ASCII nucleotides")
}

fn decode_base(bits: u8) -> u8 {
    match bits {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => unreachable!("packed kmer uses only 2-bit nucleotides"),
    }
}

#[cfg(test)]
mod tests {
    use super::{decode_kmer, mask_for_k, mutate_at_position, sample_friendly_kmers_parallel, PackedKmer};
    use ahash::AHashSet;
    use rand::rngs::SmallRng;
    use rand::SeedableRng;

    fn pack(sequence: &str) -> PackedKmer {
        sequence.bytes().fold(0u64, |acc, byte| {
            let bits = match byte {
                b'A' => 0u64,
                b'C' => 1u64,
                b'G' => 2u64,
                b'T' => 3u64,
                _ => panic!("invalid base"),
            };
            (acc << 2) | bits
        })
    }

    #[test]
    fn mask_supports_k_30() {
        assert_eq!(mask_for_k(30), (1u64 << 60) - 1);
    }

    #[test]
    fn random_search_finds_valid_friendly_kmers() {
        let source = pack("ACGT");
        let source_set = AHashSet::from_iter([source]);
        let hits = sample_friendly_kmers_parallel(source_set, vec![source], 4, 1, 4, Some(7))
            .expect("friendly search should succeed");

        assert_eq!(hits.len(), 4);
        let rendered = AHashSet::from_iter(hits.iter().map(|pair| decode_kmer(pair.friendly, 4)));
        assert_eq!(rendered.len(), 4);
        assert!(!rendered.contains("ACGT"));
    }

    #[test]
    fn center_position_mutation_changes_only_middle_base() {
        let source = pack("AAA");
        let mut seen = AHashSet::new();
        for _ in 0..64 {
            let mut rng = SmallRng::seed_from_u64(9);
            let friendly = mutate_at_position(source, 3, 1, &mut rng);
            seen.insert(decode_kmer(friendly, 3));
        }

        for kmer in seen {
            assert_eq!(&kmer[0..1], "A");
            assert_eq!(&kmer[2..3], "A");
            assert_ne!(&kmer[1..2], "A");
        }
    }
}

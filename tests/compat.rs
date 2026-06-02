use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

fn ours_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-multi-replace"))
}

fn golden(name: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), name)
}

fn oracle_script() -> String {
    format!("{}/tests/oracle_skbio.py", env!("CARGO_MANIFEST_DIR"))
}

fn parse_matrix(text: &str) -> HashMap<String, Vec<f64>> {
    text.lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let mut f = l.split('\t');
            let sample = f.next().unwrap().to_string();
            let row: Vec<f64> = f.map(|c| c.trim().parse().unwrap()).collect();
            (sample, row)
        })
        .collect()
}

fn ours(table: &str, delta: Option<&str>) -> HashMap<String, Vec<f64>> {
    let mut cmd = Command::new(ours_bin());
    cmd.args(["--input", &golden(table)]);
    if let Some(d) = delta {
        cmd.args(["--delta", d]);
    }
    let out = cmd.output().expect("run rsomics-multi-replace");
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    parse_matrix(&String::from_utf8(out.stdout).unwrap())
}

const TOL: f64 = 1e-9;

fn assert_close(a: &HashMap<String, Vec<f64>>, b: &HashMap<String, Vec<f64>>, ctx: &str) {
    assert_eq!(a.len(), b.len(), "{ctx} sample count");
    for (sample, arow) in a {
        let brow = b
            .get(sample)
            .unwrap_or_else(|| panic!("{ctx}: missing sample {sample}"));
        assert_eq!(arow.len(), brow.len(), "{ctx} {sample} width");
        for (j, (&x, &y)) in arow.iter().zip(brow).enumerate() {
            assert!(
                (x - y).abs() <= TOL,
                "{ctx} {sample}[{j}]: {x} vs {y} (|Δ|={})",
                (x - y).abs()
            );
        }
    }
}

/// Committed skbio-captured matrix — always-on regression gate, runs with no
/// scikit-bio present.
fn check_committed(table: &str, delta: Option<&str>, expected_file: &str) {
    let expected = parse_matrix(&std::fs::read_to_string(golden(expected_file)).unwrap());
    let got = ours(table, delta);
    assert_close(&got, &expected, expected_file);
}

#[test]
fn committed_small() {
    check_committed("small_table.tsv", None, "small_expected.tsv");
}

#[test]
fn committed_med() {
    check_committed("med_table.tsv", None, "med_expected.tsv");
}

#[test]
fn committed_delta() {
    check_committed("small_table.tsv", Some("0.01"), "small_expected_d01.tsv");
}

/// scikit-bio is the named oracle; loud-skip if it (or python) is unavailable.
/// `RSOMICS_SKBIO_PYTHON` overrides the interpreter (e.g. an isolated venv).
fn skbio_python() -> Option<String> {
    let mut candidates = Vec::new();
    if let Ok(p) = std::env::var("RSOMICS_SKBIO_PYTHON") {
        candidates.push(p);
    }
    candidates.push("python3".into());
    candidates.push("python".into());
    for py in candidates {
        let probe = Command::new(&py)
            .args(["-c", "import skbio.stats.composition"])
            .output();
        if let Ok(out) = probe
            && out.status.success()
        {
            return Some(py);
        }
    }
    eprintln!("SKIP: scikit-bio not importable — install `scikit-bio` to run the differential");
    None
}

fn oracle(py: &str, table: &str, delta: Option<&str>) -> HashMap<String, Vec<f64>> {
    let mut cmd = Command::new(py);
    cmd.arg(oracle_script()).arg(golden(table));
    if let Some(d) = delta {
        cmd.arg(d);
    }
    let out = cmd.output().expect("run scikit-bio oracle");
    assert!(
        out.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    parse_matrix(&String::from_utf8(out.stdout).unwrap())
}

/// Live differential: value-exact (~1e-9) vs skbio multi_replace().
fn differential(table: &str, delta: Option<&str>) {
    let Some(py) = skbio_python() else { return };
    let o = ours(table, delta);
    let t = oracle(&py, table, delta);
    assert_close(&o, &t, table);
}

#[test]
fn matches_skbio_small() {
    differential("small_table.tsv", None);
}

#[test]
fn matches_skbio_med() {
    differential("med_table.tsv", None);
}

#[test]
fn matches_skbio_delta() {
    differential("small_table.tsv", Some("0.01"));
}

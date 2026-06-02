use std::io::{BufRead, Write};

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};

mod table;

pub use table::Table;

pub struct ReplacedMatrix {
    pub samples: Vec<String>,
    pub components: Vec<String>,
    pub data: Vec<f64>,
}

/// Multiplicative replacement of zeros (Martin-Fernandez 2003), matching
/// `skbio.stats.composition.multi_replace`: close each row to sum 1, then set
/// every zero to `delta` and scale the non-zeros by `1 - n_zeros*delta` so the
/// row still sums to 1. `delta` defaults to `(1/D)^2` where `D` is the number
/// of components.
///
/// # Errors
/// A non-finite or negative input, an all-zero row (skbio's `closure` rejects
/// these), or a `delta` so large it drives a non-zero proportion below 0.
pub fn multi_replace(table: &Table, delta: Option<f64>) -> Result<ReplacedMatrix> {
    let m = table.n_components();
    let delta = delta.unwrap_or_else(|| (1.0 / m as f64).powi(2));

    let mut data = vec![0.0_f64; table.data.len()];
    let bad: Option<RowError> = data
        .par_chunks_mut(m)
        .zip(table.data.par_chunks(m))
        .map(|(out, row)| {
            let mut sum = 0.0;
            let mut zeros = 0usize;
            for &x in row {
                if !x.is_finite() {
                    return Some(RowError::NonFinite(x));
                }
                if x < 0.0 {
                    return Some(RowError::Negative(x));
                }
                sum += x;
            }
            if sum == 0.0 {
                return Some(RowError::AllZero);
            }
            for &x in row {
                if x == 0.0 {
                    zeros += 1;
                }
            }
            let zcnts = 1.0 - zeros as f64 * delta;
            if zcnts < 0.0 {
                return Some(RowError::NegativeProportion);
            }
            for (o, &x) in out.iter_mut().zip(row) {
                let closed = x / sum;
                *o = if closed == 0.0 { delta } else { zcnts * closed };
            }
            None
        })
        .reduce(|| None, |a, b| a.or(b));

    if let Some(e) = bad {
        return Err(e.into_err(delta));
    }

    Ok(ReplacedMatrix {
        samples: table.samples.clone(),
        components: table.components.clone(),
        data,
    })
}

enum RowError {
    NonFinite(f64),
    Negative(f64),
    AllZero,
    NegativeProportion,
}

impl RowError {
    fn into_err(self, delta: f64) -> RsomicsError {
        let msg = match self {
            RowError::NonFinite(v) => format!("input matrix has a non-finite value: {v}"),
            RowError::Negative(v) => format!("input matrix has a negative component: {v}"),
            RowError::AllZero => "input matrix has a composition with all zeros".into(),
            RowError::NegativeProportion => {
                format!("delta {delta} created a negative proportion — use a smaller delta")
            }
        };
        RsomicsError::InvalidInput(msg)
    }
}

impl ReplacedMatrix {
    /// Write as a TSV: a `sample` header of component IDs, then one
    /// `sample_id<delim>value...` line per sample, floats shortest-round-trip.
    ///
    /// # Errors
    /// Propagates write errors.
    pub fn write_tsv<W: Write>(&self, mut out: W, delim: char) -> Result<()> {
        let m = self.components.len();
        let mut line = String::with_capacity(16 * (m + 1));
        line.push_str("sample");
        for c in &self.components {
            line.push(delim);
            line.push_str(c);
        }
        line.push('\n');
        out.write_all(line.as_bytes()).map_err(RsomicsError::Io)?;

        let mut buf = FloatBuf::new();
        for (i, sample) in self.samples.iter().enumerate() {
            line.clear();
            line.push_str(sample);
            for &v in &self.data[i * m..(i + 1) * m] {
                line.push(delim);
                line.push_str(buf.format(v));
            }
            line.push('\n');
            out.write_all(line.as_bytes()).map_err(RsomicsError::Io)?;
        }
        Ok(())
    }
}

/// std's `{}` already emits the shortest decimal that round-trips; reuse one
/// growable buffer to avoid a per-value allocation.
struct FloatBuf {
    s: String,
}

impl FloatBuf {
    fn new() -> Self {
        Self {
            s: String::with_capacity(24),
        }
    }

    fn format(&mut self, v: f64) -> &str {
        use std::fmt::Write;
        self.s.clear();
        let _ = write!(self.s, "{v}");
        &self.s
    }
}

/// # Errors
/// Propagates parse, compute, and write errors.
pub fn run<W: Write>(
    table_reader: impl BufRead,
    out: W,
    delim: char,
    delta: Option<f64>,
) -> Result<()> {
    let table = Table::parse(table_reader, delim)?;
    let res = multi_replace(&table, delta)?;
    res.write_tsv(out, delim)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(txt: &str) -> Table {
        Table::parse(txt.as_bytes(), '\t').unwrap()
    }

    fn rows(r: &ReplacedMatrix) -> Vec<Vec<f64>> {
        let m = r.components.len();
        r.data.chunks(m).map(|c| c.to_vec()).collect()
    }

    #[test]
    fn docstring_example() {
        let t = parse("\tc1\tc2\tc3\tc4\ns1\t0.2\t0.4\t0.4\t0\ns2\t0\t0.5\t0.5\t0\n");
        let r = multi_replace(&t, None).unwrap();
        let want = [
            [0.1875, 0.375, 0.375, 0.0625],
            [0.0625, 0.4375, 0.4375, 0.0625],
        ];
        for (g, w) in rows(&r).iter().zip(&want) {
            for (a, b) in g.iter().zip(w) {
                assert!((a - b).abs() < 1e-12, "{a} vs {b}");
            }
        }
    }

    #[test]
    fn counts_get_closed_first() {
        let t = parse("\tc1\tc2\tc3\tc4\ns1\t2\t2\t6\t0\n");
        let r = multi_replace(&t, None).unwrap();
        let want = [0.1875, 0.1875, 0.5625, 0.0625];
        for (a, b) in r.data.iter().zip(&want) {
            assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        }
    }

    #[test]
    fn rows_sum_to_one() {
        let r = multi_replace(&parse("\tc1\tc2\tc3\ns1\t5\t0\t9\ns2\t0\t8\t0\n"), None).unwrap();
        for chunk in r.data.chunks(3) {
            assert!((chunk.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn no_zeros_remain() {
        let r = multi_replace(&parse("\tc1\tc2\tc3\ns1\t0\t0\t9\n"), None).unwrap();
        assert!(r.data.iter().all(|&v| v > 0.0));
    }

    #[test]
    fn custom_delta() {
        let r = multi_replace(
            &parse("\tc1\tc2\tc3\tc4\ns1\t0.2\t0.4\t0.4\t0\n"),
            Some(0.01),
        )
        .unwrap();
        let want = [0.198, 0.396, 0.396, 0.01];
        for (a, b) in r.data.iter().zip(&want) {
            assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        }
    }

    #[test]
    fn large_delta_bails() {
        // Two zeros, delta 0.6 → zcnts = 1 - 1.2 < 0.
        assert!(multi_replace(&parse("\tc1\tc2\tc3\tc4\ns1\t2\t4\t0\t0\n"), Some(0.6)).is_err());
    }

    #[test]
    fn all_zero_row_bails() {
        assert!(multi_replace(&parse("\tc1\tc2\ns1\t0\t0\n"), None).is_err());
    }

    #[test]
    fn negative_bails() {
        assert!(multi_replace(&parse("\tc1\tc2\ns1\t-1\t5\n"), None).is_err());
    }
}

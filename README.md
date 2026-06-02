# rsomics-multi-replace

Multiplicative zero-replacement of a compositional matrix, as a single fast CLI.
Equivalent to `skbio.stats.composition.multi_replace` (formerly
`multiplicative_replacement`).

Compositional log-ratio transforms (CLR, ILR) are undefined on zeros. This is
the standard preprocessing step that removes them: each row is closed to sum 1,
every zero is set to a small `delta`, and the surviving non-zeros are scaled by
`1 - n_zeros·delta` so the row still sums to 1.

```
rsomics-multi-replace --input table.tsv [--delta F] [-o replaced.tsv]
```

- `table.tsv` — composition table: header row of component IDs (corner cell
  ignored), then one `sample_id  value...` line per sample. Counts or
  proportions; each row is closed internally.
- `--delta` — value substituted for each zero. Default `(1/D)^2` where `D` is
  the number of components, matching skbio. A `delta` large enough to push a
  non-zero proportion below 0 is rejected.
- `--csv` — comma-separated I/O instead of tab.

Output is the replaced matrix: a `sample` header of component IDs, then one
`sample_id<TAB>value...` line per sample. Pipe it straight into `rsomics-clr` or
`rsomics-ilr`.

Rows are independent, so `-t` parallelises across samples.

## Origin

This crate is an independent Rust reimplementation of
`skbio.stats.composition.multi_replace` based on:

- Martín-Fernández, J. A., Barceló-Vidal, C., Pawlowsky-Glahn, V., "Dealing With
  Zeros and Missing Values in Compositional Data Sets Using Nonparametric
  Imputation", *Mathematical Geology* 35(3), 2003. DOI: 10.1023/A:1023866030544
- The scikit-bio implementation (Modified BSD License), read and cited:
  `multi_replace` calls `closure(mat)`, marks zeros, then returns
  `where(zero, delta, (1 - n_zeros·delta)·closed)` with `delta = (1/D)²` by
  default, raising if any resulting proportion is negative.

Values are bit-reproducible vs scikit-bio (pure closed-form arithmetic, no RNG,
no iteration); `tests/compat.rs` diffs a committed skbio-captured golden and a
live `skbio.stats.composition.multi_replace` run to ~1e-9.

License: MIT OR Apache-2.0.
Upstream credit: scikit-bio https://scikit-bio.org/ (Modified BSD License).

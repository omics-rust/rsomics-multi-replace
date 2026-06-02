use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_multi_replace::{Table, multi_replace};

fn synth(n_samples: usize, n_components: usize) -> Table {
    let mut txt = String::new();
    for j in 0..n_components {
        txt.push('\t');
        txt.push_str(&format!("c{j}"));
    }
    txt.push('\n');
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % 100
    };
    for s in 0..n_samples {
        txt.push_str(&format!("s{s}"));
        for _ in 0..n_components {
            txt.push('\t');
            // ~20% zeros, the rest small counts.
            txt.push_str(&format!("{}", next().saturating_sub(20)));
        }
        txt.push('\n');
    }
    Table::parse(txt.as_bytes(), '\t').unwrap()
}

fn bench(c: &mut Criterion) {
    let table = synth(500, 2000);
    c.bench_function("multi_replace_500x2000", |b| {
        b.iter(|| multi_replace(black_box(&table), None).unwrap())
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);

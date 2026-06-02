#!/usr/bin/env python3
"""skbio.stats.composition.multi_replace oracle.

Usage: oracle_skbio.py table.tsv [delta] > replaced.tsv
Emits the zero-replaced matrix as `sample<TAB>value...` matching rsomics-multi-replace.
"""
import sys

import numpy as np
import pandas as pd
from skbio.stats.composition import multi_replace

table_path = sys.argv[1]
delta = float(sys.argv[2]) if len(sys.argv) > 2 else None

table = pd.read_csv(table_path, sep="\t", index_col=0)
mat = table.to_numpy(dtype=float)
out = np.atleast_2d(multi_replace(mat, delta=delta))

cols = list(table.columns)
sys.stdout.write("sample\t" + "\t".join(cols) + "\n")
for sample, row in zip(table.index, out):
    sys.stdout.write(str(sample) + "\t" + "\t".join(repr(float(v)) for v in row) + "\n")

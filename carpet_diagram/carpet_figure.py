"""
Carpet diagram (portrait plot) of model-vs-reconstruction temperature RMSE.

An update of IPCC AR5 WG1 Fig. 9.12 for a single variable (annual-mean
temperature). Columns = models, rows = (period, reconstruction compilation),
cells coloured by RMSE in degrees C.

Reads every `output/rmse_long_<period>.csv` produced by the benchmarking
notebook, so adding further periods (midHolocene, lgm, midPliocene-eoi400)
later makes new rows appear automatically — no edits needed here.

Run with the `my-cli-py` conda env, from the carpet_diagram/ directory.
"""

import os
import glob
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import Normalize

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, 'output')

# Preferred period ordering (oldest→youngest is not meaningful here; this is
# just the row order). Any period found but not listed is appended after these.
PERIOD_ORDER = ['midPliocene-eoi400', 'lgm', 'lig127k', 'midHolocene']
COMPILATION_ORDER = ['Hoffman', 'Capron']

# Cells backed by fewer than this many proxy points are flagged (small samples
# give noisy RMSE — e.g. the Capron lig127k subset has only ~7 points).
MIN_POINTS = 10


def load_rmse_tables():
    """Concatenate every rmse_long_<period>.csv into one tidy frame."""
    rows = []
    for f in sorted(glob.glob(os.path.join(OUT_DIR, 'rmse_long_*.csv'))):
        period = re.match(r'rmse_long_(.+)\.csv', os.path.basename(f)).group(1)
        df = pd.read_csv(f)
        df['period'] = period
        rows.append(df)
    if not rows:
        raise FileNotFoundError(f'no rmse_long_*.csv in {OUT_DIR}')
    return pd.concat(rows, ignore_index=True)


def order_rows(df):
    """Return ordered list of (period, compilation) row keys present in df."""
    periods = ([p for p in PERIOD_ORDER if p in set(df['period'])]
               + sorted(set(df['period']) - set(PERIOD_ORDER)))
    comps = ([c for c in COMPILATION_ORDER if c in set(df['compilation'])]
             + sorted(set(df['compilation']) - set(COMPILATION_ORDER)))
    return [(p, c) for p in periods for c in comps
            if not df[(df.period == p) & (df.compilation == c)].empty]


def make_figure():
    df = load_rmse_tables()

    row_keys = order_rows(df)
    # Columns: models sorted by mean RMSE across all rows (best on the left).
    model_order = (df.groupby('model')['rmse'].mean().sort_values().index.tolist())

    rmse = df.pivot_table(index=['period', 'compilation'], columns='model', values='rmse')
    npts = df.pivot_table(index=['period', 'compilation'], columns='model', values='n_points')
    rmse = rmse.reindex(index=row_keys, columns=model_order)
    npts = npts.reindex(index=row_keys, columns=model_order)

    M = rmse.values
    nrows, ncols = M.shape

    fig_w = max(8.0, 0.55 * ncols + 3.0)
    fig_h = max(2.6, 0.6 * nrows + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    norm = Normalize(vmin=np.nanmin(M), vmax=np.nanmax(M))
    cmap = plt.get_cmap('YlOrRd').copy()
    cmap.set_bad('#e8e8e8')

    im = ax.imshow(M, aspect='auto', cmap=cmap, norm=norm)

    # Cell annotations: RMSE value, with a '*' where the sample is small.
    for i in range(nrows):
        for j in range(ncols):
            v = M[i, j]
            if np.isnan(v):
                continue
            n = npts.values[i, j]
            txt = f'{v:.1f}' + ('*' if np.isfinite(n) and n < MIN_POINTS else '')
            # White text on dark cells, black on light.
            lum = 0.0 if np.isnan(v) else norm(v)
            color = 'white' if lum > 0.6 else '#222222'
            ax.text(j, i, txt, ha='center', va='center', fontsize=7.5, color=color)

    # Axes / ticks
    ax.set_xticks(range(ncols))
    ax.set_xticklabels(model_order, rotation=55, ha='right', fontsize=8)
    ax.set_yticks(range(nrows))
    ax.set_yticklabels([f'{p}\n{c}' for p, c in row_keys], fontsize=8.5)
    ax.set_xticks(np.arange(-.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-.5, nrows, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1.5)
    ax.tick_params(which='minor', length=0)
    ax.tick_params(which='major', length=0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label('RMSE of annual-mean temperature anomaly (°C)', fontsize=8.5)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title('Model–reconstruction temperature mismatch (carpet diagram)',
                 fontsize=12, fontweight='bold', pad=10)
    if (npts.values[np.isfinite(npts.values)] < MIN_POINTS).any():
        fig.text(0.01, 0.01, f'* fewer than {MIN_POINTS} proxy points — RMSE is noisy',
                 fontsize=7.5, color='#555555', ha='left')

    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, 'carpet_diagram.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    print('Saved →', out_path)


if __name__ == '__main__':
    make_figure()

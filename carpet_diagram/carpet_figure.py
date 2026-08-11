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
from matplotlib.colors import Normalize, LinearSegmentedColormap

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, 'output')


def better_worse_cmap():
    """Diverging colormap for RMSE-relative-to-row-mean.

    Below-mean (better than average) → warm YlOrRd (pale yellow → red); above-mean
    (worse) → cool BuPu (pale → blue → purple); the mean sits at the pale centre.
    """
    n = 128
    better = plt.get_cmap('YlOrRd')(np.linspace(1.0, 0.0, n))  # red → pale (best → mean)
    worse = plt.get_cmap('BuPu')(np.linspace(0.0, 1.0, n))     # pale → purple (mean → worst)
    return LinearSegmentedColormap.from_list('better_worse', np.vstack([better, worse]))

# Explicit row order, in two groups separated by a double rule in the figure:
# the proxy compilations first, then the gridded data assimilation products.
# Any (period, compilation) present in the data but not listed here is appended
# as a further group rather than silently dropped.
ROW_GROUPS = [
    [('midHolocene', 'Bartlein'),
     ('midHolocene', 'Temp12k'),
     ('lgm', 'Bartlein'),
     ('lig127k', 'Hoffman'),
     ('lig127k', 'Capron'),
     ('midPliocene-eoi400', 'Foley-Dowsett')],
    # Cleator is a data assimilation product too (it assimilates the Bartlein
    # 21 ka pollen synthesis above into PMIP3 output), so it belongs here.
    [('midHolocene', 'Erb'),
     ('lgm', 'Cleator'),
     ('lgm', 'Osman'),
     ('midPliocene-eoi400', 'Tierney')],
]

# Cells backed by fewer than this many proxy points are flagged (small samples
# give noisy RMSE — e.g. the Capron lig127k subset has only ~7 points).
MIN_POINTS = 10

# The CMIP6-era PMIP models, from the PMIP data holdings table at
# https://pcmdi.llnl.gov/CMIP6/ArchiveStatistics/esgf_data_holdings/PMIP/index.html
# (EC-Earth3-Veg and EC-Earth3-Veg-LR are listed there but have no CVDP output).
CMIP6_PMIP_MODELS = [
    'ACCESS-ESM1-5', 'AWI-ESM-1-1-LR', 'CESM2', 'CESM2-FV2', 'CESM2-WACCM-FV2',
    'CNRM-CM6-1', 'EC-Earth3', 'EC-Earth3-LR', 'EC-Earth3-Veg', 'EC-Earth3-Veg-LR',
    'FGOALS-f3-L', 'FGOALS-g3', 'GISS-E2-1-G', 'HadGEM3-GC31-LL', 'INM-CM4-8',
    'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NESM3',
    'NorESM1-F', 'NorESM2-LM',
]
# Kept as their own columns in the collapsed figure; everything else is averaged.
KEEP_MODELS = frozenset(CMIP6_PMIP_MODELS) | {'UofT-CCSM-4'}
COLLAPSED_LABEL = 'PMIP3'


def collapse_old_models(df):
    """Average every model outside KEEP_MODELS into a single COLLAPSED_LABEL column.

    The average is the plain mean of those models' RMSE within each row — not the
    RMSE of their ensemble-mean anomaly field, which would be lower because model
    errors partly cancel. `n_points` is carried through as the median so the
    small-sample '*' flag still behaves.
    """
    keep = df[df.model.isin(KEEP_MODELS)].copy()
    old = df[~df.model.isin(KEEP_MODELS)]
    if old.empty:
        return keep, []

    agg = (old.groupby(['period', 'compilation'])
              .agg(rmse=('rmse', 'mean'), n_points=('n_points', 'median'))
              .reset_index())
    agg['model'] = COLLAPSED_LABEL
    return pd.concat([keep, agg], ignore_index=True), sorted(old.model.unique())


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
    """Return (row_keys, boundaries) for the (period, compilation) rows in df.

    `row_keys` is the flat ordered row list; `boundaries` holds the row indices
    after which a group ends (used to draw the double rule between groups).
    """
    present = set(map(tuple, df[['period', 'compilation']].drop_duplicates().values))
    groups = [[k for k in g if k in present] for g in ROW_GROUPS]

    listed = {k for g in ROW_GROUPS for k in g}
    extra = sorted(present - listed)
    if extra:
        print('note: rows not in ROW_GROUPS, appended at the end:', extra)
        groups.append(extra)

    groups = [g for g in groups if g]
    row_keys = [k for g in groups for k in g]
    boundaries = np.cumsum([len(g) for g in groups])[:-1].tolist()
    return row_keys, boundaries


def make_figure(out_name='carpet_diagram.png', collapse=False, title=None):
    df = load_rmse_tables()

    collapsed_models = []
    if collapse:
        df, collapsed_models = collapse_old_models(df)

    row_keys, boundaries = order_rows(df)

    rmse = df.pivot_table(index=['period', 'compilation'], columns='model', values='rmse')
    npts = df.pivot_table(index=['period', 'compilation'], columns='model', values='n_points')
    rmse = rmse.reindex(index=row_keys)

    # Per-row normalisation relative to the row mean: each cell is the model's signed
    # deviation from that (period, compilation) row's average RMSE, scaled by the row's
    # largest absolute deviation so the most extreme model reaches ±1. Better-than-
    # average (lower RMSE) is negative → warm (yellow/red); worse-than-average is
    # positive → cool (blue/purple); the mean sits at the pale centre. This shows model
    # ranking regardless of each row's absolute RMSE magnitude. Cells keep absolute °C.
    def relative_to_mean(r):
        dev = r - r.mean()
        scale = np.nanmax(np.abs(dev))
        return dev / scale if np.isfinite(scale) and scale > 0 else dev * 0.0

    rel = rmse.apply(relative_to_mean, axis=1)
    # Columns: models sorted by mean relative RMSE (best-on-average on the left).
    # The collapsed multi-model column is not a model, so it is pinned to the far
    # right instead of competing in the ranking.
    model_order = rel.mean(axis=0).sort_values().index.tolist()
    if collapse and COLLAPSED_LABEL in model_order:
        model_order = [m for m in model_order if m != COLLAPSED_LABEL] + [COLLAPSED_LABEL]
    rmse = rmse.reindex(columns=model_order)
    rel = rel.reindex(columns=model_order)
    npts = npts.reindex(index=row_keys, columns=model_order)

    M = rmse.values
    Mn = rel.values
    nrows, ncols = M.shape

    fig_w = max(8.0, 0.55 * ncols + 3.0)
    fig_h = max(2.6, 0.6 * nrows + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    norm = Normalize(vmin=-1.0, vmax=1.0)
    cmap = better_worse_cmap()
    cmap.set_bad('#e8e8e8')

    im = ax.imshow(Mn, aspect='auto', cmap=cmap, norm=norm)

    # Cell annotations: absolute RMSE value, with a '*' where the sample is small.
    for i in range(nrows):
        for j in range(ncols):
            v = M[i, j]
            if np.isnan(v):
                continue
            n = npts.values[i, j]
            txt = f'{v:.1f}' + ('*' if np.isfinite(n) and n < MIN_POINTS else '')
            # White text on saturated (far-from-mean) cells, black near the pale centre.
            color = 'white' if abs(Mn[i, j]) > 0.6 else '#222222'
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

    # Double rule between row groups (proxy compilations above, data
    # assimilation products below).
    for b in boundaries:
        for dy in (-0.085, 0.085):
            ax.axhline(b - 0.5 + dy, color='#222222', linewidth=1.2, zorder=5)

    # Matching double rule before the collapsed multi-model column.
    if collapse and COLLAPSED_LABEL in model_order:
        xb = model_order.index(COLLAPSED_LABEL) - 0.5
        for dx in (-0.085, 0.085):
            ax.axvline(xb + dx, color='#222222', linewidth=1.2, zorder=5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, ticks=[-1, 0, 1])
    cbar.set_label('RMSE vs row mean', fontsize=8.5)
    cbar.ax.set_yticklabels(['better', 'mean', 'worse'])
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(title or 'Model–reconstruction temperature mismatch (carpet diagram)',
                 fontsize=12, fontweight='bold', pad=10)
    notes = ['cell numbers are absolute RMSE (°C); colour is each model vs its row-mean RMSE']
    if (npts.values[np.isfinite(npts.values)] < MIN_POINTS).any():
        notes.append(f'* fewer than {MIN_POINTS} proxy points — RMSE is noisy')
    fig.text(0.01, 0.01, '   '.join(notes), fontsize=7.5, color='#555555', ha='left')
    if collapsed_models:
        fig.text(0.01, -0.02,
                 f'{COLLAPSED_LABEL} column = mean of the RMSEs of whichever of these '
                 f'{len(collapsed_models)} pre-CMIP6 models ran that period (mean of '
                 'RMSEs, not RMSE of the ensemble mean): ' + ', '.join(collapsed_models),
                 fontsize=6.5, color='#555555', ha='left', wrap=True)

    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, out_name)
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    print('Saved →', out_path)


if __name__ == '__main__':
    make_figure()
    make_figure(out_name='carpet_diagram_cmip6.png', collapse=True,
                title='Model–reconstruction temperature mismatch — CMIP6-era PMIP models')

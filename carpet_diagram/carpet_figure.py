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

# Ensemble-mean summary columns, in the order they are pinned to the right-hand
# edge. PMIP4 averages the CMIP6 archive models only — UofT-CCSM-4 keeps its own
# column and is deliberately *not* folded in, since it is not in the CMIP6
# holdings. Move it into the PMIP4 membership test below if you want it counted.
COMPOSITE_LABELS = ['PMIP4', 'PMIP3']


def composite_membership(df):
    """Map each summary column label -> the models it averages, for this data."""
    have = set(df.model.unique())
    return {
        'PMIP4': sorted(have & set(CMIP6_PMIP_MODELS)),
        'PMIP3': sorted(have - KEEP_MODELS),
    }


def add_group_means(df):
    """Keep KEEP_MODELS as columns and append the PMIP4 / PMIP3 mean columns.

    Each summary is the plain mean of its members' RMSE within a row — not the
    RMSE of their ensemble-mean anomaly field, which would be lower because model
    errors partly cancel. `n_points` is carried through as the median so the
    small-sample '*' flag still behaves.
    """
    members = composite_membership(df)
    frames = [df[df.model.isin(KEEP_MODELS)].copy()]
    for label in COMPOSITE_LABELS:
        sub = df[df.model.isin(members[label])]
        if sub.empty:
            continue
        agg = (sub.groupby(['period', 'compilation'])
                  .agg(rmse=('rmse', 'mean'), n_points=('n_points', 'median'))
                  .reset_index())
        agg['model'] = label
        frames.append(agg)
    return pd.concat(frames, ignore_index=True), members


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

    members = {}
    if collapse:
        df, members = add_group_means(df)
    composites = [c for c in COMPOSITE_LABELS if c in set(df.model)]

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
    #
    # The baseline is the *individual models only*: the PMIP4/PMIP3 summary columns are
    # placed on that scale but excluded from defining it, so the colours mean "vs the
    # typical model" rather than "vs a mix of models and ensemble means". A consequence
    # worth knowing when reading the figure: PMIP4 is by construction the mean of most
    # of those same models, so it sits near the pale centre in every row — that is
    # arithmetic, not a result. PMIP3's colour is the informative one.
    individual = [m for m in rmse.columns if m not in composites]

    def relative_to_mean(r):
        dev = r - r[individual].mean()
        scale = np.nanmax(np.abs(dev[individual]))
        return dev / scale if np.isfinite(scale) and scale > 0 else dev * 0.0

    rel = rmse.apply(relative_to_mean, axis=1)
    # Columns: models sorted by mean relative RMSE (best-on-average on the left).
    # The summary columns are not models, so they are pinned to the far right in
    # COMPOSITE_LABELS order instead of competing in the ranking.
    model_order = rel[individual].mean(axis=0).sort_values().index.tolist() + composites
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

    # Matching double rule before the block of summary columns.
    if composites:
        xb = model_order.index(composites[0]) - 0.5
        for dx in (-0.085, 0.085):
            ax.axvline(xb + dx, color='#222222', linewidth=1.2, zorder=5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, ticks=[-1, 0, 1])
    cbar.set_label('RMSE vs row mean', fontsize=8.5)
    cbar.ax.set_yticklabels(['better', 'mean', 'worse'])
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(title or 'Model–reconstruction temperature mismatch (carpet diagram)',
                 fontsize=12, fontweight='bold', pad=10)
    notes = ['cell numbers are absolute RMSE (°C); colour is each model vs the row-mean RMSE '
             'of the individual models']
    if (npts.values[np.isfinite(npts.values)] < MIN_POINTS).any():
        notes.append(f'* fewer than {MIN_POINTS} proxy points — RMSE is noisy')
    fig.text(0.01, 0.01, '   '.join(notes), fontsize=7.5, color='#555555', ha='left')

    if composites:
        lines = ['summary columns are the mean of their members\' RMSEs in each row (mean of '
                 'RMSEs, not RMSE of the ensemble mean), over whichever members ran that period;'
                 ' they do not contribute to the colour scale']
        for label in composites:
            lines.append(f'{label} ({len(members[label])}): ' + ', '.join(members[label]))
        for k, line in enumerate(lines):
            fig.text(0.01, -0.018 - 0.019 * k, line, fontsize=6.5, color='#555555',
                     ha='left', wrap=True)

    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, out_name)
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    print('Saved →', out_path)


if __name__ == '__main__':
    # The headline figure is the reduced one: a column per CMIP6-era PMIP model
    # (plus UofT-CCSM-4), with PMIP4 and PMIP3 ensemble-mean columns on the right.
    make_figure(out_name='carpet_diagram.png', collapse=True)
    # The every-model version is kept alongside it.
    make_figure(out_name='carpet_diagram_all_models.png',
                title='Model–reconstruction temperature mismatch — all models')

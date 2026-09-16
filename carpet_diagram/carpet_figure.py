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
    """Diverging colormap for RMSE relative to the reference column.

    Better than the reference → warm YlOrRd (white → yellow → red); worse → cool
    BuPu (white → blue → purple). The inner end of each ramp is faded to pure
    white so the reference value sits on an exactly white centre, rather than on
    the pale-yellow/pale-blue seam the two colormaps would otherwise meet at.
    """
    n = 128
    better = plt.get_cmap('YlOrRd')(np.linspace(1.0, 0.0, n))  # red → pale (best → centre)
    worse = plt.get_cmap('BuPu')(np.linspace(0.0, 1.0, n))     # pale → purple (centre → worst)

    white = np.array([1.0, 1.0, 1.0, 1.0])
    k = n // 5                                   # fade the innermost fifth to white
    t = np.linspace(0.0, 1.0, k)[:, None]
    better[n - k:] = better[n - k:] * (1 - t) + white * t
    worse[:k] = white * (1 - t) + worse[:k] * t

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
     ('lgm', 'Annan'),
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


def add_group_means(df, labels=COMPOSITE_LABELS, keep_only=True):
    """Append the requested ensemble-mean columns, optionally dropping old models.

    `labels` selects which of COMPOSITE_LABELS to build; `keep_only` restricts the
    individual columns to KEEP_MODELS (the reduced figure) rather than showing
    every model (the all-models figure). Membership is always worked out from the
    full table, so PMIP4 averages the same models either way.

    Each summary is the plain mean of its members' RMSE within a row — not the
    RMSE of their ensemble-mean anomaly field, which would be lower because model
    errors partly cancel. `n_points` is carried through as the median so the
    small-sample '*' flag still behaves.
    """
    members = composite_membership(df)
    frames = [df[df.model.isin(KEEP_MODELS)].copy() if keep_only else df.copy()]
    for label in labels:
        sub = df[df.model.isin(members[label])]
        if sub.empty:
            continue
        agg = (sub.groupby(['period', 'compilation'])
                  .agg(rmse=('rmse', 'mean'), n_points=('n_points', 'median'))
                  .reset_index())
        agg['model'] = label
        frames.append(agg)
    return pd.concat(frames, ignore_index=True), {k: members[k] for k in labels}


def load_rmse_tables(rmse_dir=None):
    """Concatenate every rmse_long_<period>.csv into one tidy frame.

    `rmse_dir` defaults to output/; pass a subdirectory to draw the same figure
    from a subsetted table (see carpet_regional.py), which keeps those files out
    of the glob that builds the headline figure.
    """
    rmse_dir = rmse_dir or OUT_DIR
    rows = []
    for f in sorted(glob.glob(os.path.join(rmse_dir, 'rmse_long_*.csv'))):
        period = re.match(r'rmse_long_(.+)\.csv', os.path.basename(f)).group(1)
        df = pd.read_csv(f)
        df['period'] = period
        rows.append(df)
    if not rows:
        raise FileNotFoundError(f'no rmse_long_*.csv in {rmse_dir}')
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


def make_figure(out_name='carpet_diagram.png', labels=COMPOSITE_LABELS,
                keep_only=True, title=None, rmse_dir=None, note=None):
    df = load_rmse_tables(rmse_dir)

    df, members = add_group_means(df, labels=labels, keep_only=keep_only)
    composites = [c for c in COMPOSITE_LABELS if c in set(df.model)]

    row_keys, boundaries = order_rows(df)

    rmse = df.pivot_table(index=['period', 'compilation'], columns='model', values='rmse')
    npts = df.pivot_table(index=['period', 'compilation'], columns='model', values='n_points')
    rmse = rmse.reindex(index=row_keys)

    # Per-row normalisation: each cell is its signed deviation from a reference RMSE,
    # scaled by the row's largest absolute deviation so the most extreme column reaches
    # ±1. Below the reference (better) is negative → warm (yellow/red); above (worse) is
    # positive → cool (blue/purple); the reference itself is the pale centre. This shows
    # model ranking regardless of each row's absolute RMSE. Cells keep absolute °C.
    #
    # The reference is the PMIP4 ensemble mean where that column exists, so every cell
    # reads directly as "better or worse than the CMIP6 multi-model mean" and PMIP4 is
    # white by definition. The all-models figure has no PMIP4 column and falls back to
    # the row mean. The scale spans every column including the summaries, so nothing
    # clips even when PMIP3 is the most extreme value in a row.
    individual = [m for m in rmse.columns if m not in composites]
    reference = 'PMIP4' if 'PMIP4' in rmse.columns else None

    def relative_to_reference(r):
        centre = r['PMIP4'] if reference else r[individual].mean()
        if not np.isfinite(centre):
            return r * np.nan
        dev = r - centre
        scale = np.nanmax(np.abs(dev))
        return dev / scale if np.isfinite(scale) and scale > 0 else dev * 0.0

    rel = rmse.apply(relative_to_reference, axis=1)
    # Columns: the summary columns first, at the left-hand edge in COMPOSITE_LABELS
    # order, then the individual models alphabetically.
    model_order = composites + sorted(individual)
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

    # Matching double rule after the block of summary columns on the left.
    if composites:
        xb = len(composites) - 0.5
        for dx in (-0.085, 0.085):
            ax.axvline(xb + dx, color='#222222', linewidth=1.2, zorder=5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, ticks=[-1, 0, 1])
    cbar.set_label(f'RMSE vs {reference} mean' if reference else 'RMSE vs row mean',
                   fontsize=8.5)
    cbar.ax.set_yticklabels(['better', reference or 'mean', 'worse'])
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(title or 'Model–reconstruction temperature mismatch (carpet diagram)',
                 fontsize=12, fontweight='bold', pad=10)
    centre_txt = (f'the {reference} ensemble-mean RMSE (so {reference} is white by definition)'
                  if reference else 'the row-mean RMSE of the individual models')
    notes = [f'cell numbers are absolute RMSE (°C); colour is each column vs {centre_txt}']
    if note:
        notes.insert(0, note)
    if (npts.values[np.isfinite(npts.values)] < MIN_POINTS).any():
        notes.append(f'* fewer than {MIN_POINTS} proxy points — RMSE is noisy')
    fig.text(0.01, 0.01, '   '.join(notes), fontsize=7.5, color='#555555', ha='left')

    if composites:
        lines = ['summary columns are the mean of their members\' RMSEs in each row (mean of '
                 'RMSEs, not RMSE of the ensemble mean), over whichever members ran that period']
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
    # (plus UofT-CCSM-4), with PMIP4 and PMIP3 ensemble-mean columns on the left.
    make_figure(out_name='carpet_diagram.png')
    # The every-model version keeps the PMIP4 reference column so its colours are
    # centred the same way; PMIP3 is left off because all of its members are
    # already shown as their own columns. Add 'PMIP3' to labels if you want it.
    make_figure(out_name='carpet_diagram_all_models.png',
                labels=['PMIP4'], keep_only=False,
                title='Model–reconstruction temperature mismatch — all models')

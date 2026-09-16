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

    Better than the reference (lower RMSE) → green; worse → red; the reference
    itself sits on an exactly white centre. Each half is one sequential ramp
    (Greens, Reds) with its inner end faded to pure white, so the two meet at
    white rather than at the pale seam the raw colormaps would give.
    """
    n = 128
    better = plt.get_cmap('Greens')(np.linspace(1.0, 0.0, n))  # dark green → pale (best → centre)
    worse = plt.get_cmap('Reds')(np.linspace(0.0, 1.0, n))     # pale → dark red (centre → worst)

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
    # The P2F recommended-SST records, kept in their own group below the other
    # proxy compilations and above the data assimilation products.
    [('midHolocene', 'P2F'),
     ('lgm', 'P2F')],
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
    # Every per-model value column is averaged the same way; `model_mean` (the
    # regional figure's annotation) rides along with `rmse` when it is present.
    value_cols = [c for c in ('rmse', 'model_mean', 'bias') if c in df.columns]
    for label in labels:
        sub = df[df.model.isin(members[label])]
        if sub.empty:
            continue
        agg = {c: (c, 'mean') for c in value_cols}
        agg['n_points'] = ('n_points', 'median')
        agg = sub.groupby(['period', 'compilation']).agg(**agg).reset_index()
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
                keep_only=True, title=None, rmse_dir=None, note=None,
                value_col='rmse', value_label=None, recon_col=False):
    """Draw the carpet diagram.

    `value_col` chooses what the printed cell numbers are; the colour is always
    the RMSE relative to the reference column, whatever is printed. The regional
    figure prints `model_mean` (the weighted mean model anomaly over that row's
    reconstruction points) instead of the RMSE. `recon_col` adds a column at the
    right-hand edge holding the matching weighted mean of the reconstructions,
    which is a property of the row rather than of any model.
    """
    df = load_rmse_tables(rmse_dir)
    if value_col not in df.columns:
        raise KeyError(f'{value_col!r} is not a column of the rmse_long tables '
                       f'(have: {", ".join(df.columns)})')
    recon_by_row = (df.groupby(['period', 'compilation'])['recon_mean'].first()
                    if recon_col else None)
    if recon_col and recon_by_row is None:
        raise KeyError('recon_col=True needs a recon_mean column')

    df, members = add_group_means(df, labels=labels, keep_only=keep_only)
    composites = [c for c in COMPOSITE_LABELS if c in set(df.model)]

    row_keys, boundaries = order_rows(df)

    rmse = df.pivot_table(index=['period', 'compilation'], columns='model', values='rmse')
    npts = df.pivot_table(index=['period', 'compilation'], columns='model', values='n_points')
    vals = (rmse if value_col == 'rmse' else
            df.pivot_table(index=['period', 'compilation'], columns='model', values=value_col))
    rmse = rmse.reindex(index=row_keys)
    vals = vals.reindex(index=row_keys)

    # Per-row normalisation: each cell is its signed deviation from a reference RMSE,
    # scaled by the row's largest absolute deviation so the most extreme column reaches
    # ±1. Below the reference (better) is negative → green; above (worse) is
    # positive → red; the reference itself is the white centre. This shows
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
    vals = vals.reindex(columns=model_order)
    npts = npts.reindex(index=row_keys, columns=model_order)

    M = vals.values
    Mn = rel.values
    nrows, ncols = M.shape

    fig_w = max(8.0, 0.55 * (ncols + bool(recon_col)) + 3.0)
    fig_h = max(2.6, 0.6 * nrows + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    norm = Normalize(vmin=-1.0, vmax=1.0)
    cmap = better_worse_cmap()
    cmap.set_bad('#e8e8e8')

    im = ax.imshow(Mn, aspect='auto', cmap=cmap, norm=norm)

    # Cell annotations: the printed value (RMSE by default, the mean model anomaly
    # in the regional figure), with a '*' where the sample is small.
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

    # The reconstruction's own weighted mean, in a gutter column past the right
    # edge of the heatmap. It is a property of the row, not of any model, so it
    # gets no cell colour — it is the target the model numbers are aiming at.
    if recon_col:
        for i, key in enumerate(row_keys):
            v = recon_by_row.get(key, np.nan)
            if np.isfinite(v):
                ax.text(ncols, i, f'{v:.1f}', ha='center', va='center',
                        fontsize=7.5, color='#222222', fontweight='bold')
        ax.set_xlim(-0.5, ncols + 0.5)
        for dx in (-0.085, 0.085):
            ax.axvline(ncols - 0.5 + dx, color='#222222', linewidth=1.2, zorder=5)

    # Axes / ticks
    ax.set_xticks(range(ncols + bool(recon_col)))
    ax.set_xticklabels(model_order + (['reconstruction'] if recon_col else []),
                       rotation=55, ha='right', fontsize=8)
    ax.set_yticks(range(nrows))
    ax.set_yticklabels([f'{p}\n{c}' for p, c in row_keys], fontsize=8.5)
    ax.set_xticks(np.arange(-.5, ncols, 1), minor=True)  # gridlines on the heatmap only
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
    value_txt = value_label or 'absolute RMSE (°C)'
    notes = [f'cell numbers are {value_txt}; colour is each column\'s RMSE vs {centre_txt}']
    if recon_col:
        notes.append('the right-hand column is the same weighted mean of the reconstructions '
                     '— the value the models are aiming at')
    if note:
        notes.insert(0, note)
    if (npts.values[np.isfinite(npts.values)] < MIN_POINTS).any():
        notes.append(f'* fewer than {MIN_POINTS} proxy points — RMSE is noisy')
    # One note per line: joining them end to end makes the saved bbox as wide as
    # the text, which stretches the whole figure.
    for k, line in enumerate(notes):
        fig.text(0.01, -0.006 - 0.021 * k, line, fontsize=7.5, color='#555555', ha='left')

    if composites:
        lines = ['summary columns are the mean of their members\' RMSEs in each row (mean of '
                 'RMSEs, not RMSE of the ensemble mean), over whichever members ran that period']
        for label in composites:
            lines.append(f'{label} ({len(members[label])}): ' + ', '.join(members[label]))
        base = -0.006 - 0.021 * len(notes) - 0.008
        for k, line in enumerate(lines):
            fig.text(0.01, base - 0.019 * k, line, fontsize=6.5, color='#555555',
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

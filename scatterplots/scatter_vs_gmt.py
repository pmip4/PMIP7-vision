"""
Scatter plots of climate-mode / precipitation responses against global-mean
temperature change, one point per simulation.

Reads `tidy_numbers_per_simulations.csv` (produced by compile_tidy_numbers.ncl,
one row per model_experiment CVDP file). Every quantity is converted to an
anomaly from that same model's piControl:

  - global-mean precipitation  -> percentage change, (exp-pi)/pi*100
  - all variability modes       -> plain difference,  exp-pi
  - x-axis (global-mean temp)   -> plain difference,  exp-pi   (degrees C)

Six panels: precipitation (with a 2 %/degC reference line), ENSO (Nino3.4
s.d.), Atlantic Nino (ATL3 s.d.), NAO, PDO and AMM. Points are coloured by
experiment.

Run with the `my-cli-py` conda env, from the scatterplots/ directory.
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(SCRIPT_DIR, 'tidy_numbers_per_simulations.csv')
OUT_DIR = os.path.join(SCRIPT_DIR, 'output')

# NCL writes missing data as one of these fill sentinels; anything this large
# is not a physical value.
FILL_THRESH = 1e15

# Experiments to plot (paleo periods + the 1pctCO2 future case), each paired
# with its model's piControl. The `-cal-adj` variants are deliberately excluded
# to avoid double-counting; abrupt4xCO2 rows are all fill values so they drop
# out on their own.
EXP_STYLE = {
    'midPliocene-eoi400': ('#b2182b', 'mid-Pliocene'),
    'lgm':                ('#2166ac', 'LGM'),
    'lig127k':            ('#ef8a62', 'LIG (127k)'),
    'midHolocene':        ('#1a9850', 'mid-Holocene'),
    '1pctCO2':            ('#762a83', '1% CO2'),
    'abrupt4xCO2':        ('#a40e4C', 'Abrupt 4x CO2')
}

# Panels: (column, title, y-axis label, mode). mode 'pct' -> percentage change,
# 'diff' -> exp minus piControl.
PANELS = [
    ('global_mean_precipitation', 'Global-mean precipitation', 'change (%)', 'pct'),
    ('Nino34 stddev',             'ENSO (Nino3.4 s.d.)',        r'$\Delta$ s.d. ($\degree$C)', 'diff'),
    ('AtlNino stddev',               'Atlantic Nino (ATL3 s.d.)',  r'$\Delta$ s.d. ($\degree$C)', 'diff'),
    ('NAO DJF',                   'NAO (DJF)',                  r'$\Delta$ s.d. (hPa)',        'diff'),
    ('IOD stddev',                'Indian Ocean Dipole (DMI)',  r'$\Delta$ s.d. ($\degree$C)', 'diff'),
    ('AMM AnnCycAmp',             'Atlantic Meridional Mode ',  r'$\Delta$ ($\degree$C)',      'diff'),
]

PRECIP_COL = 'global_mean_precipitation'
GMT_COL = 'GMTemp'


def load_clean():
    """Load the CSV, masking fill sentinels and physically absurd precip."""
    df = pd.read_csv(CSV)
    numeric = df.columns.difference(['model', 'experiment'])
    df[numeric] = df[numeric].apply(pd.to_numeric, errors='coerce')
    df[numeric] = df[numeric].mask(df[numeric].abs() >= FILL_THRESH)
    # A couple of rows carry a garbage precip value (e.g. 2.6e5) that is below
    # the fill threshold; global-mean precip is ~2-4 mm/day, so clip anything
    # outside a generous physical range.
    df[PRECIP_COL] = df[PRECIP_COL].where(df[PRECIP_COL].between(0, 100))
    return df


def anomalies(df):
    """Return a long frame of (model, experiment, dGMT, panel anomalies)."""
    pi = df[df.experiment == 'piControl'].drop_duplicates('model').set_index('model')
    records = []
    for _, row in df[df.experiment.isin(EXP_STYLE)].iterrows():
        model = row['model']
        if model not in pi.index:
            continue
        ctrl = pi.loc[model]
        dgmt = row[GMT_COL] - ctrl[GMT_COL]
        if not np.isfinite(dgmt):
            continue
        rec = {'model': model, 'experiment': row['experiment'], 'dGMT': dgmt}
        for col, _title, _ylab, mode in PANELS:
            e, c = row[col], ctrl[col]
            if not (np.isfinite(e) and np.isfinite(c)):
                rec[col] = np.nan
            elif mode == 'pct':
                rec[col] = (e - c) / c * 100.0 if c != 0 else np.nan
            else:
                rec[col] = e - c
        records.append(rec)
    return pd.DataFrame.from_records(records)


def make_figure():
    df = load_clean()
    an = anomalies(df)

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    axes = axes.ravel()

    for ax, (col, title, ylab, mode) in zip(axes, PANELS):
        ax.axhline(0, color='#bbbbbb', lw=0.8, zorder=0)
        ax.axvline(0, color='#bbbbbb', lw=0.8, zorder=0)

        for exp, (color, _label) in EXP_STYLE.items():
            sub = an[an.experiment == exp]
            ax.scatter(sub['dGMT'], sub[col], s=42, c=color,
                       edgecolor='white', linewidth=0.5, alpha=0.9, zorder=3)

        # 2 %/degC reference line on the precipitation panel.
        if mode == 'pct':
            xlim = ax.get_xlim()
            xs = np.array(xlim)
            ax.plot(xs, 2.0 * xs, color='#333333', ls='--', lw=1.2, zorder=2,
                    label='2 %/$\\degree$C')
            ax.set_xlim(xlim)
            ax.legend(loc='upper left', fontsize=8, frameon=False)

        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r'Global-mean $\Delta$T ($\degree$C)', fontsize=9)
        ax.set_ylabel(ylab, fontsize=9)
        ax.tick_params(labelsize=8)

    # Shared experiment legend along the bottom.
    handles = [plt.Line2D([0], [0], marker='o', ls='', mec='white', mew=0.5,
                          mfc=color, ms=8, label=label)
               for color, label in EXP_STYLE.values()]
    fig.legend(handles=handles, loc='lower center', ncol=len(EXP_STYLE),
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.01))

    fig.suptitle('Simulated responses vs global-mean temperature change '
                 '(anomalies from each model\'s piControl)',
                 fontsize=13, fontweight='bold')
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, 'scatter_vs_gmt.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    print('Saved ->', out_path)


if __name__ == '__main__':
    make_figure()

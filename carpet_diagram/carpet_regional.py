"""
Regional variant of the carpet diagram.

Same benchmarking as `carpet_figure.py`, but restricted to the reconstruction
points inside a latitude/longitude box — by default 45–66.33 °N, 30 °W–15 °E
(NW Europe and the NE Atlantic, the northern edge being the Arctic Circle).

Nothing needs re-running upstream: the notebooks already write
`output/model_anom_at_recon_<period>.csv`, which holds every proxy point with
its coordinates, its reconstructed anomaly, its weight and each model's anomaly
sampled at that point. This script re-does Step 4 (the weighted RMSE) over the
subset of rows inside the box, writes the resulting tables to
`output/region_<name>/`, and hands them to `carpet_figure.make_figure`.

The regional tables live in their own directory so that the `rmse_long_*.csv`
glob behind the headline figure never picks them up.

Run with the `my-cli-py` conda env, from the carpet_diagram/ directory:

    python carpet_regional.py                         # the default box
    python carpet_regional.py --lat 45 66.33 --lon -30 15 --name euroatlantic
"""

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

import carpet_figure as cf

OUT_DIR = cf.OUT_DIR

# Columns of model_anom_at_recon_<period>.csv that describe the point rather
# than a model; everything else in the header is a model column.
META_COLS = ['compilation', 'reference', 'site', 'Proxy', 'Latitude', 'Longitude',
             'recon_anom', 'weight']

DEFAULT_REGION = dict(name='euroatlantic', lat=(45.0, 66.33), lon=(-30.0, 15.0))


def fmt_bounds(lat, lon):
    """Human-readable box label, e.g. '45–66.33 °N, 30 °W–15 °E'."""
    def ns(v):
        return f'{abs(v):g} °{"N" if v >= 0 else "S"}'

    def ew(v):
        return f'{abs(v):g} °{"E" if v >= 0 else "W"}'

    return f'{lat[0]:g}–{ns(lat[1])}, {ew(lon[0])}–{ew(lon[1])}'


def weighted_mean(values, weights):
    """Weighted mean over the finite entries; NaN if nothing is finite."""
    ok = np.isfinite(values) & np.isfinite(weights)
    if not ok.any():
        return np.nan
    return float(np.sum(weights[ok] * values[ok]) / np.sum(weights[ok]))


def regional_rmse(path, lat, lon):
    """Weighted RMSE, bias and mean anomaly per (model, compilation) in the box.

    `rmse` and `bias` mirror Step 4 of the notebooks exactly — the same
    weight-weighted statistics over the points where both the sampled model
    value and the weight are finite — but over the boxed subset of rows.

    Two extra columns drive the regional figure, which prints anomalies rather
    than errors: `model_mean` is the weighted mean of the model's own anomaly
    over those same points, and `recon_mean` is the weighted mean of the
    reconstruction over every boxed point of that compilation. `recon_mean` is a
    property of the row, so it repeats down the model column.
    """
    df = pd.read_csv(path)
    models = [c for c in df.columns if c not in META_COLS]

    inside = (df.Latitude.between(*lat)) & (df.Longitude.between(*lon))
    sub = df[inside]

    recon_mean = {comp: weighted_mean(g['recon_anom'].values, g['weight'].values)
                  for comp, g in sub.groupby('compilation')}

    records = []
    for m in models:
        for comp, grp in sub.groupby('compilation'):
            diff = grp[m].values - grp['recon_anom'].values
            w = grp['weight'].values
            valid = np.isfinite(diff) & np.isfinite(w)
            n = int(valid.sum())
            if n:
                dv, wv = diff[valid], w[valid]
                rmse = float(np.sqrt(np.sum(wv * dv ** 2) / np.sum(wv)))
                bias = float(np.sum(wv * dv) / np.sum(wv))
                # The model mean uses the same points as the RMSE, so the two
                # always describe the same sample.
                model_mean = float(np.sum(wv * grp[m].values[valid]) / np.sum(wv))
            else:
                rmse = bias = model_mean = np.nan
            records.append({'model': m, 'compilation': comp, 'n_points': n,
                            'rmse': rmse, 'bias': bias, 'model_mean': model_mean,
                            'recon_mean': recon_mean[comp]})
    return pd.DataFrame(records), sorted(df.compilation.unique()), len(sub), len(df)


def build_tables(region_dir, lat, lon):
    """Write one rmse_long_<period>.csv per period into region_dir."""
    os.makedirs(region_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(OUT_DIR, 'model_anom_at_recon_*.csv')))
    if not files:
        raise FileNotFoundError(f'no model_anom_at_recon_*.csv in {OUT_DIR} — '
                                'run the notebooks first')
    for f in files:
        period = re.match(r'model_anom_at_recon_(.+)\.csv', os.path.basename(f)).group(1)
        rmse_long, all_comps, n_in, n_all = regional_rmse(f, lat, lon)
        counts = (rmse_long.groupby('compilation').n_points.max()
                  .sort_index().to_dict()) if len(rmse_long) else {}
        # A compilation with no points in the box has nothing to benchmark
        # against, so its row is dropped rather than drawn empty. Say so out
        # loud — a silently missing row is easy to misread as a missing period.
        empty = [c for c in all_comps if counts.get(c, 0) == 0]
        if empty:
            print(f'  {period}: dropped (no points in box): ' + ', '.join(empty))
            rmse_long = rmse_long[~rmse_long.compilation.isin(empty)]

        rmse_long.to_csv(os.path.join(region_dir, f'rmse_long_{period}.csv'), index=False)
        kept = {c: n for c, n in counts.items() if n}
        print(f'  {period}: {n_in}/{n_all} points in box — ' +
              ', '.join(f'{c} {n}' for c, n in kept.items()))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lat', nargs=2, type=float, metavar=('SOUTH', 'NORTH'),
                   default=DEFAULT_REGION['lat'], help='latitude bounds, °N')
    p.add_argument('--lon', nargs=2, type=float, metavar=('WEST', 'EAST'),
                   default=DEFAULT_REGION['lon'],
                   help='longitude bounds, °E (proxy convention, −180–180)')
    p.add_argument('--name', default=DEFAULT_REGION['name'],
                   help='slug used in the output directory and figure names')
    args = p.parse_args()

    lat, lon = tuple(args.lat), tuple(args.lon)
    region_dir = os.path.join(OUT_DIR, f'region_{args.name}')
    label = fmt_bounds(lat, lon)

    print(f'Region {args.name}: {label}')
    build_tables(region_dir, lat, lon)

    note = (f'restricted to the {label} reconstruction points '
            '(rows are the same compilations, subsetted — not re-fitted)')
    # The regional figure prints the mean anomaly rather than the RMSE: over a
    # single region the anomaly itself is the readable quantity, and the
    # reconstruction column gives it something to be read against. Colour still
    # comes from the RMSE, so the ranking is unchanged.
    shared = dict(rmse_dir=region_dir, note=note, value_col='model_mean',
                  value_label='the weighted mean model anomaly (°C) over that row\'s '
                              'reconstruction points', recon_col=True)
    cf.make_figure(out_name=f'carpet_diagram_{args.name}.png', **shared,
                   title=f'Model–reconstruction temperature anomalies — {label}')
    cf.make_figure(out_name=f'carpet_diagram_{args.name}_all_models.png',
                   labels=['PMIP4'], keep_only=False, **shared,
                   title=f'Model–reconstruction temperature anomalies, all models — {label}')


if __name__ == '__main__':
    main()

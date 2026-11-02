"""
Regional variant of the carpet diagram.

Same benchmarking as `carpet_figure.py`, but restricted to the reconstruction
points inside a region — either a latitude/longitude box (the default: 45–66.33
°N, 30 °W–15 °E, NW Europe and the NE Atlantic, the northern edge being the
Arctic Circle) or a named IPCC AR6 reference region, e.g. `--ar6 MED`.

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
    python carpet_regional.py --ar6 MED               # an IPCC AR6 region
"""

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
import regionmask
import shapely

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


class Region:
    """A named selection of points: `inside(lat, lon)` returns a boolean array."""

    def __init__(self, name, label, inside):
        self.name = name
        self.label = label
        self.inside = inside


def box_region(lat, lon, name):
    """Points inside a plain latitude/longitude box, bounds inclusive."""
    def inside(lat_v, lon_v):
        return ((lat_v >= lat[0]) & (lat_v <= lat[1]) &
                (lon_v >= lon[0]) & (lon_v <= lon[1]))

    return Region(name, fmt_bounds(lat, lon), inside)


def ar6_region(abbrev, name=None):
    """Points inside an IPCC AR6 reference region, by its abbreviation.

    The polygons come from `regionmask.defined_regions.ar6.all`, the 58 regions
    of Iturbide et al. (2020) used by the AR6 WGI Atlas — the same source the
    synthesis figure uses. `.all` rather than `.land`, so combined land-and-ocean
    regions such as MED keep their marine part and the SST compilations are not
    silently dropped. Longitudes there run −180–180, matching the proxy
    convention of the point CSVs.
    """
    regions = regionmask.defined_regions.ar6.all
    match = [r for r in regions if r.abbrev.upper() == abbrev.upper()]
    if not match:
        raise SystemExit(f'no AR6 region {abbrev!r}. Available: '
                         + ', '.join(sorted(r.abbrev for r in regions)))
    region = match[0]
    poly = region.polygon

    def inside(lat_v, lon_v):
        return shapely.contains_xy(poly, np.asarray(lon_v), np.asarray(lat_v))

    return Region(name or region.abbrev.lower(),
                  f'IPCC AR6 {region.abbrev} ({region.name})', inside)


def weighted_mean(values, weights):
    """Weighted mean over the finite entries; NaN if nothing is finite."""
    ok = np.isfinite(values) & np.isfinite(weights)
    if not ok.any():
        return np.nan
    return float(np.sum(weights[ok] * values[ok]) / np.sum(weights[ok]))


def regional_rmse(path, region):
    """Weighted RMSE, bias and mean anomaly per (model, compilation) in region.

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

    sub = df[region.inside(df.Latitude.values, df.Longitude.values)]

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


def build_tables(region_dir, region):
    """Write one rmse_long_<period>.csv per period into region_dir."""
    os.makedirs(region_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(OUT_DIR, 'model_anom_at_recon_*.csv')))
    if not files:
        raise FileNotFoundError(f'no model_anom_at_recon_*.csv in {OUT_DIR} — '
                                'run the notebooks first')
    for f in files:
        period = re.match(r'model_anom_at_recon_(.+)\.csv', os.path.basename(f)).group(1)
        rmse_long, all_comps, n_in, n_all = regional_rmse(f, region)
        counts = (rmse_long.groupby('compilation').n_points.max()
                  .sort_index().to_dict()) if len(rmse_long) else {}
        # A compilation with no points in the region has nothing to benchmark
        # against, so its row is dropped rather than drawn empty. Say so out
        # loud — a silently missing row is easy to misread as a missing period.
        empty = [c for c in all_comps if counts.get(c, 0) == 0]
        if empty:
            print(f'  {period}: dropped (no points in region): ' + ', '.join(empty))
            rmse_long = rmse_long[~rmse_long.compilation.isin(empty)]

        rmse_long.to_csv(os.path.join(region_dir, f'rmse_long_{period}.csv'), index=False)
        kept = {c: n for c, n in counts.items() if n}
        print(f'  {period}: {n_in}/{n_all} points in region — ' +
              ', '.join(f'{c} {n}' for c, n in kept.items()))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lat', nargs=2, type=float, metavar=('SOUTH', 'NORTH'),
                   default=DEFAULT_REGION['lat'], help='latitude bounds, °N')
    p.add_argument('--lon', nargs=2, type=float, metavar=('WEST', 'EAST'),
                   default=DEFAULT_REGION['lon'],
                   help='longitude bounds, °E (proxy convention, −180–180)')
    p.add_argument('--ar6', metavar='ABBREV',
                   help='select an IPCC AR6 reference region by abbreviation '
                        '(e.g. MED) instead of a box; --lat/--lon are ignored')
    p.add_argument('--name', default=None,
                   help='slug used in the output directory and figure names '
                        '(defaults to the AR6 abbreviation, or euroatlantic)')
    args = p.parse_args()

    if args.ar6:
        region = ar6_region(args.ar6, args.name)
    else:
        region = box_region(tuple(args.lat), tuple(args.lon),
                            args.name or DEFAULT_REGION['name'])
    region_dir = os.path.join(OUT_DIR, f'region_{region.name}')

    print(f'Region {region.name}: {region.label}')
    build_tables(region_dir, region)

    label = region.label
    note = (f'restricted to the reconstruction points inside {label} '
            '(rows are the same compilations, subsetted — not re-fitted)')
    # The regional figure prints the mean anomaly rather than the RMSE: over a
    # single region the anomaly itself is the readable quantity, and the
    # reconstruction column gives it something to be read against. Colour still
    # comes from the RMSE, so the ranking is unchanged.
    shared = dict(rmse_dir=region_dir, note=note, value_col='model_mean',
                  value_label='the weighted mean model anomaly (°C) over that row\'s '
                              'reconstruction points', recon_col=True)
    cf.make_figure(out_name=f'carpet_diagram_{region.name}.png', **shared,
                   title=f'Model–reconstruction temperature anomalies — {label}')
    cf.make_figure(out_name=f'carpet_diagram_{region.name}_all_models.png',
                   labels=['PMIP4'], keep_only=False, **shared,
                   title=f'Model–reconstruction temperature anomalies, all models — {label}')


if __name__ == '__main__':
    main()

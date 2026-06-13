"""Generate carpet_diagram_lig127k.ipynb via nbformat. Run once, then nbconvert --execute."""
import os
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))
co = lambda s: cells.append(nbf.v4.new_code_cell(s))

md("""# Carpet diagram — lig127k temperature benchmarking

First step toward an update of IPCC AR5 WG1 Fig. 9.12 (a *portrait plot* / *carpet
diagram*). For each PMIP lig127k simulation we:

1. Load the annual-mean near-surface temperature field `tas_spatialmean_ann` from the
   CVDP output, and the matching `piControl` run.
2. Compute the **lig127k − piControl anomaly** on each model's native grid.
3. **Sample** that anomaly field at the location of every proxy reconstruction.
4. Summarise the model–data mismatch as the **root-mean-squared error (RMSE)** against
   each of the two reconstruction compilations (Hoffman and Capron, from
   Otto-Bliesner et al. 2021).

Intermediate tables are written to `output/`. Run with the `my-cli-py` conda env.""")

co("""import os, glob, re
import numpy as np
import pandas as pd
import xarray as xr

ROOT      = os.getcwd()  # run the notebook from the carpet_diagram/ directory
CVDP_DIR  = os.path.join(ROOT, 'cvdp_output_by_experiment')
RECON_DIR = os.path.join(ROOT, 'recons', 'lig127k')
OUT_DIR   = os.path.join(ROOT, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

EXPERIMENT = 'lig127k'
VAR        = 'tas_spatialmean_ann'
print('CVDP   :', CVDP_DIR)
print('recons :', RECON_DIR)
print('output :', OUT_DIR)""")

md("""## Step 1 — Reconstruction compilations

The three CSVs from the Otto-Bliesner et al. (2021) supplement (Tables S2–S4) split the
proxies by region and have slightly different column layouts, but all share
`Compilation reference`, `Latitude`, `Longitude` and `Anom` (the central annual-mean
temperature anomaly estimate). We normalise them and split into the two compilations.""")

co("""recon_files = sorted(glob.glob(os.path.join(RECON_DIR, 'mat_*.csv')))
assert recon_files, 'no reconstruction CSVs found'

frames = []
for f in recon_files:
    df = pd.read_csv(f, header=2)
    site_col = 'Site Name' if 'Site Name' in df.columns else 'Core'
    sub = df[['Compilation reference', site_col, 'Latitude', 'Longitude', 'Proxy', 'Anom']].copy()
    sub = sub.rename(columns={'Compilation reference': 'reference', site_col: 'site'})
    sub['source_table'] = os.path.basename(f)
    frames.append(sub)

recon = pd.concat(frames, ignore_index=True)
# Keep only rows with a usable location and anomaly
recon = recon.dropna(subset=['Latitude', 'Longitude', 'Anom'])
recon[['Latitude', 'Longitude', 'Anom']] = recon[['Latitude', 'Longitude', 'Anom']].astype(float)

# Two compilations: Hoffman and Capron
recon['compilation'] = np.where(recon['reference'].str.contains('Hoffman', case=False),
                                'Hoffman',
                                np.where(recon['reference'].str.contains('Capron', case=False),
                                         'Capron', 'other'))
assert (recon['compilation'] != 'other').all(), recon.loc[recon.compilation=='other','reference'].unique()

print(recon.groupby('compilation').size())
recon_out = os.path.join(OUT_DIR, 'recon_points_lig127k.csv')
recon.to_csv(recon_out, index=False)
print('wrote', recon_out)
recon.head()""")

md("""## Step 2 — Model lig127k − piControl anomalies

Each model writes its CVDP field on its own native grid, so anomalies are computed
per model. We pair every `lig127k` file with the same model's `piControl` file
(ignoring the auxiliary `.monsoon.` / `.tas.indices.` files). If the two grids ever
differ, the control is bilinearly regridded onto the lig127k grid before differencing.""")

co("""def model_files(experiment):
    \"\"\"Map model name -> path for the main CVDP file of an experiment.\"\"\"
    out = {}
    pat = os.path.join(CVDP_DIR, experiment, f'*_{experiment}.cvdp_data.*.nc')
    for f in sorted(glob.glob(pat)):
        b = os.path.basename(f)
        if '.monsoon.' in b or '.tas.indices.' in b or b.endswith('.tmp'):
            continue
        model = b.split(f'_{experiment}.cvdp_data')[0]
        out[model] = f  # last (sorted) wins; there is one main file per model
    return out

lig_files = model_files('lig127k')
pi_files  = model_files('piControl')
models = sorted(set(lig_files) & set(pi_files))
print(f'{len(models)} models with both lig127k and piControl:')
print(models)
missing_pi = sorted(set(lig_files) - set(pi_files))
if missing_pi:
    print('lig127k models with no piControl (skipped):', missing_pi)""")

co("""def load_field(path):
    da = xr.open_dataset(path, decode_times=False)[VAR].sortby('lat').sortby('lon')
    # A few CVDP grids (e.g. LOVECLIM piControl) carry duplicate lon values,
    # which break interpolation; keep the first occurrence of each coordinate.
    for dim in ('lat', 'lon'):
        _, idx = np.unique(da[dim].values, return_index=True)
        if len(idx) != da.sizes[dim]:
            da = da.isel({dim: np.sort(idx)})
    return da

anomalies = {}
for m in models:
    lig = load_field(lig_files[m])
    pi  = load_field(pi_files[m])
    if lig.shape != pi.shape or not (np.allclose(lig.lat, pi.lat) and np.allclose(lig.lon, pi.lon)):
        pi = pi.interp(lat=lig.lat, lon=lig.lon)
    anomalies[m] = (lig - pi).rename('tas_anom')
    print(f'{m:18s} grid {lig.shape}  mean anom {float(anomalies[m].mean()):+.2f} C')""")

md("""## Step 3 — Sample model anomalies at reconstruction locations

Model longitudes run 0–360°, the proxy longitudes −180–180°, so targets are wrapped to
0–360 and the field is made cyclic in longitude before bilinear interpolation.""")

co("""def sample_points(field, lats, lons):
    \"\"\"Bilinearly sample a (lat, lon) field at scattered points; lon made cyclic.\"\"\"
    lon_cyc = np.append(field.lon.values, field.lon.values[0] + 360.0)
    fcyc = xr.concat([field, field.isel(lon=0)], dim='lon').assign_coords(lon=lon_cyc)
    ta = xr.DataArray(np.asarray(lats), dims='point')
    to = xr.DataArray(np.asarray(lons) % 360.0, dims='point')
    return fcyc.interp(lat=ta, lon=to).values

sampled = recon[['compilation', 'reference', 'site', 'Proxy',
                 'Latitude', 'Longitude', 'Anom']].copy()
sampled = sampled.rename(columns={'Anom': 'recon_anom'})
for m in models:
    sampled[m] = sample_points(anomalies[m], sampled['Latitude'].values, sampled['Longitude'].values)

sampled_out = os.path.join(OUT_DIR, 'model_anom_at_recon_lig127k.csv')
sampled.to_csv(sampled_out, index=False)
print('wrote', sampled_out, '  shape', sampled.shape)
sampled.head()""")

md("""## Step 4 — RMSE of each model against each compilation

For every model × compilation we take the model-minus-proxy difference across all proxy
points in that compilation and report the RMSE (and the number of points contributing,
since a point off the model grid edge can return NaN).""")

co("""records = []
for m in models:
    for comp, grp in sampled.groupby('compilation'):
        diff = grp[m].values - grp['recon_anom'].values
        valid = np.isfinite(diff)
        n = int(valid.sum())
        rmse = float(np.sqrt(np.mean(diff[valid] ** 2))) if n else np.nan
        bias = float(np.mean(diff[valid])) if n else np.nan
        records.append({'model': m, 'compilation': comp, 'n_points': n,
                        'rmse': rmse, 'bias': bias})

rmse_long = pd.DataFrame(records)
rmse_wide = rmse_long.pivot(index='model', columns='compilation', values='rmse')
rmse_wide.columns = [f'{c}_RMSE' for c in rmse_wide.columns]

rmse_long.to_csv(os.path.join(OUT_DIR, 'rmse_long_lig127k.csv'), index=False)
rmse_wide.to_csv(os.path.join(OUT_DIR, 'rmse_summary_lig127k.csv'))
print('wrote rmse_long_lig127k.csv and rmse_summary_lig127k.csv')
rmse_wide.sort_values(rmse_wide.columns[0])""")

md("""These RMSE values are the building blocks of the carpet diagram: one column per model,
one row per reconstruction compilation, coloured by RMSE. The next step (a later figure
script) will extend this to the other periods (midHolocene, lgm, midPliocene-eoi400) and
render the portrait plot.""")

nb['cells'] = cells
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'carpet_diagram_lig127k.ipynb')
nbf.write(nb, out)
print('wrote', out)

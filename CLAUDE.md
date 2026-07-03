# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PMIP7-vision produces synthesis figures for the PMIP7 (Paleoclimate Modelling Intercomparison Project phase 7). It contains two independent figure efforts:

- `synthesis_figure/` — an IPCC-style hexagonal region plot of PMIP4 climate anomalies across AR6 land regions.
- `carpet_diagram/` — a model-vs-reconstruction benchmarking "portrait plot" / carpet diagram (an update of IPCC AR5 WG1 Fig. 9.12), quantifying each simulation's RMSE from proxy reconstruction compilations.

## Running the Figure Script

Use the `my-cli-py` conda environment (the system python lacks the scientific stack):

```bash
source activate my-cli-py
cd synthesis_figure
python hex_figure.py
```

Output is written to `synthesis_figure/output/synthesis_hexfig.png`. Paths in the script are derived from the script's own location, so it can be run from anywhere.

## Key Dependencies

`numpy`, `matplotlib`, `xarray`, `shapely`, `regionmask`

AR6 region geometries come from `regionmask.defined_regions.ar6.land` (46 regions); regionmask downloads its shapefile to `~/.cache/regionmask/` automatically on first use.

## Architecture of `hex_figure.py`

Data flows through three stages:

1. **Region geometry** — `regionmask.defined_regions.ar6.land` supplies the 46 AR6 land/land-ocean regions and their abbreviations.

2. **Regional means** (`compute_regional_means`) — opens the NetCDF ensemble files in `synthesis_data_files/`, averages over the `model` dimension, applies the AR6 mask, and computes cosine-latitude-weighted means for each region × scenario × variable combination. Returns a nested dict `results[scenario][abbrev][var]`.

3. **Figure rendering** — places each region as a hexagon on a schematic geographic grid (`HEX_LAYOUT`), styled after the IPCC Interactive Atlas: hexes pack tightly into continent-shaped clusters with continent labels. Each hex has 8 wedges (left half = temperature, right half = precipitation); scenario pairs mirror about the vertical divider, ordered top→bottom mid-Pliocene / LIG / mid-Holocene / LGM. Wedges are angular sectors clipped to the hexagon with shapely (`wedge_polygon`), so rims follow the hex edges. Colour encodes anomaly direction and magnitude via `TwoSlopeNorm` diverging colormaps (`RdBu_r` for temperature, `BrBG` for precipitation).

## Data Files

Eight NetCDF files in `synthesis_figure/synthesis_data_files/`, one per variable (tas, pr) × scenario (lgm, midHolocene, lig127k, midPliocene-eoi400). Precipitation files are pre-computed percentage changes (`_percentage.nc`); temperature files are absolute anomalies in °C.

## Hex Layout

`HEX_LAYOUT` in `hex_figure.py` maps each AR6 abbreviation to `(col, row)` offset-grid coordinates. Odd rows are shifted right by half a hex width (pointy-top hexagons). Hexes are placed so each continent forms a contiguous touching cluster; `CONTINENT_LABELS` holds the label positions in plot coordinates.

## Open Design Questions

See `synthesis_figure/README.md` for unresolved decisions on: signal/confidence encoding (arrows vs. fill opacity), region subset, and where AR6 regional mean values will ultimately come from for PMIP7 (currently using PMIP4 data as a placeholder).

# carpet_diagram/

Benchmarks each model's annual-mean temperature anomaly against proxy reconstruction compilations and reports RMSE — the building block of the carpet/portrait plot. Implemented for all four periods: `midPliocene-eoi400`, `lgm`, `lig127k`, `midHolocene`.

## Building/Running the Notebooks

`build_notebook.py` is the source of truth — edit it, not the `.ipynb` files. It generates one notebook per period (`carpet_diagram_<period>.ipynb`) from a shared template: Steps 2–4 (model anomalies, sampling, RMSE) are identical across periods, and only Step 1 (reconstruction parsing) is period-specific, defined in the `PERIODS` list. Regenerate and execute with:

```bash
source activate my-cli-py
cd carpet_diagram
python build_notebook.py        # writes carpet_diagram_<period>.ipynb for all four periods
jupyter nbconvert --to notebook --execute --inplace carpet_diagram_midPliocene-eoi400.ipynb
jupyter nbconvert --to notebook --execute --inplace carpet_diagram_lgm.ipynb
jupyter nbconvert --to notebook --execute --inplace carpet_diagram_lig127k.ipynb
jupyter nbconvert --to notebook --execute --inplace carpet_diagram_midHolocene.ipynb
```

The notebooks read their working directory via `os.getcwd()`, so run them from `carpet_diagram/`. Add a new period by appending a config (with its Step 1 recon-parsing cell) to `PERIODS`.

## Data Sources

- **CVDP model output**: `cvdp_output_by_experiment/` is a **symlink** to `/data/oacd/CVDP/by_experiment/`. One NetCDF per simulation named `model_experiment.cvdp_data.years.nc`; the only field used is `tas_spatialmean_ann` (annual-mean near-surface temperature, °C). Each model is on its **own native grid**; the experiment and piControl share a model's grid. Auxiliary `.monsoon.` / `.tas.indices.` / `.tmp` files must be filtered out. Some grids (e.g. LOVECLIM piControl) carry duplicate longitude values that break interpolation — `load_field` de-duplicates coordinates. A few CVDP files lack `tas_spatialmean_ann` entirely (e.g. EC-EARTH-2-2, KCM1-2-2, NorESM2-LM midHolocene) — Step 2 checks for the variable and drops those models, logging which. Use the plain `midHolocene` CVDP subdirectory, **not** `midHolocene-cal-adj`.
- **Reconstructions**: each period has its own `recons/<period>/` directory; values are already anomalies, and each period's Step 1 cell normalises its files to the shared `compilation, reference, site, Proxy, Latitude, Longitude, Anom` schema.
  - `recons/lig127k/` holds three region-split CSVs from the Otto-Bliesner et al. (2021) supplement (Tables S2–S4), header row 2. Column layouts differ but all share `Compilation reference`, `Latitude`, `Longitude`, `Anom`. Rows split into two compilations by reference string: **Hoffman** (78 pts, global SST) and **Capron** (7 pts, high-latitude).
  - `recons/midHolocene/` holds two compilations, one CSV each, plain header: **Bartlein** (`Bartlein_mat.csv`, 635 pts, pollen MAT, anomaly column `mat_anm_mean`) and **Temp12k** (`temp12k_anom6k_mat.csv`, 335 pts, multiproxy, anomaly column `anom`). Both use plain `lat`/`lon` columns; some Bartlein headers carry leading spaces (stripped on load).
  - `recons/lgm/` holds two near-global **gridded** anomaly products (so Step 1 sets `weight = cos(latitude)` for an area-fair RMSE): **Cleator** (`cleator2020_recon.csv`, ~2214 land grid cells; a long metadata preamble precedes a `#`-prefixed header at line 77 read with `skiprows=77, header=None`; the MAT-anomaly column is `MAT`) and **Osman** (`Osman_LGMR_21ka_SAT_anom_climo.nc`, the LGMR ensemble-mean `sat` field, already a 21 ka − (0–1 ka) anomaly; the 96×144 grid is flattened to ~13824 points, dropping fill values). Note the Osman file was regenerated as an anomaly (`ncdiff` of the 21 ka and 0–1 ka climatologies) — the original `..._SAT_climo.nc` was absolute temperature.
  - `recons/midPliocene-eoi400/` holds one compilation: **Foley-Dowsett** (`FoleyDowsett2019_cs_mp_sst_data_30k_plus_NOAA.csv`, 39 Uk'37 SST sites, 37 usable). Following Haywood et al. (2020) the anomaly is the **NOAA ERSST5** column — the last, unnamed column (`Unnamed: 14`), i.e. SST minus the NOAA ERSST5 modern climatology, *not* the HadISST-based `Anom` column. The `Supplementary/`, `.zip`, and `.pdf` alongside it are the untouched supplement bundle the CSV was extracted from.

## Pipeline

1. Parse + normalise the recon files (Step 1, period-specific) into the shared `compilation, reference, site, Proxy, Latitude, Longitude, Anom` schema (plus an optional `weight` column). Sources vary by period: CSVs with region-split/preamble/plain headers, and one gridded NetCDF (Osman) flattened to points.
2. Pair each experiment model with its piControl, compute `<experiment> − piControl` anomaly on the native grid (regrid control if grids ever differ); drop models whose CVDP file lacks the tas field.
3. Bilinearly sample each model's anomaly field at every proxy location — model lon is 0–360°, proxy lon −180–180°, so targets are wrapped and the field is made cyclic in longitude first. A `weight` column is carried through (defaults to 1.0 when Step 1 omits it).
4. Compute per-model weight-weighted RMSE (and bias) against each compilation. With unit weights this is the ordinary RMSE (all periods except lgm); the gridded lgm compilations use `cos(latitude)` weights so the near-global RMSE is area-fair rather than pole-heavy.

Outputs land in `carpet_diagram/output/`, one set per period: `recon_points_<period>.csv`, `model_anom_at_recon_<period>.csv`, `rmse_long_<period>.csv`, `rmse_summary_<period>.csv`. `carpet_figure.py` auto-discovers every `rmse_long_*.csv` and renders all periods/compilations as rows — no edit needed when a period is added. Colour is **normalised per row relative to the row mean**: each cell is the model's signed deviation from that row's average RMSE (scaled by the row's largest absolute deviation), shown on a diverging map — better-than-average (lower RMSE) is warm yellow/red, worse-than-average is cool blue/purple, and the row mean is the pale centre (`better_worse_cmap`). This keeps model ranking legible despite very different absolute RMSE between periods; the printed cell numbers remain absolute RMSE in °C, and columns are ordered by each model's mean relative RMSE (best-on-average on the left).

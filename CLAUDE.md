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

Benchmarks each model's annual-mean temperature anomaly against proxy reconstruction compilations and reports RMSE — the building block of the carpet/portrait plot. Currently implemented for `lig127k` only; the long-term target also covers midHolocene, lgm, midPliocene-eoi400.

## Building/Running the Notebook

`carpet_diagram_lig127k.ipynb` is generated from `build_notebook.py` (the script is the source of truth — edit it, not the `.ipynb` directly). Regenerate and execute with:

```bash
source activate my-cli-py
cd carpet_diagram
python build_notebook.py        # writes carpet_diagram_lig127k.ipynb
jupyter nbconvert --to notebook --execute --inplace carpet_diagram_lig127k.ipynb
```

The notebook reads its working directory via `os.getcwd()`, so run it from `carpet_diagram/`.

## Data Sources

- **CVDP model output**: `cvdp_output_by_experiment/` is a **symlink** to `/data/oacd/CVDP/by_experiment/`. One NetCDF per simulation named `model_experiment.cvdp_data.years.nc`; the only field used is `tas_spatialmean_ann` (annual-mean near-surface temperature, °C). Each model is on its **own native grid**; lig127k and piControl share a model's grid. Auxiliary `.monsoon.` / `.tas.indices.` / `.tmp` files must be filtered out. Some grids (e.g. LOVECLIM piControl) carry duplicate longitude values that break interpolation — `load_field` de-duplicates coordinates.
- **Reconstructions**: `recons/lig127k/` holds three region-split CSVs from the Otto-Bliesner et al. (2021) supplement (Tables S2–S4). Their column layouts differ but all share `Compilation reference`, `Latitude`, `Longitude`, `Anom`. Rows split into two compilations by reference string: **Hoffman** (78 pts, global SST) and **Capron** (7 pts, high-latitude). Values are already anomalies.

## Pipeline

1. Parse + normalise the recon CSVs (header row 2), drop blank trailing rows, split into Hoffman/Capron.
2. Pair each lig127k model with its piControl, compute `lig127k − piControl` anomaly on the native grid (regrid control if grids ever differ).
3. Bilinearly sample each model's anomaly field at every proxy location — model lon is 0–360°, proxy lon −180–180°, so targets are wrapped and the field is made cyclic in longitude first.
4. Compute per-model RMSE (and bias) against each compilation.

Outputs land in `carpet_diagram/output/`: `recon_points_lig127k.csv`, `model_anom_at_recon_lig127k.csv`, `rmse_long_lig127k.csv`, `rmse_summary_lig127k.csv`.

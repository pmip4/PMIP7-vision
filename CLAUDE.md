# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PMIP7-vision produces synthesis figures for the PMIP7 (Paleoclimate Modelling Intercomparison Project phase 7). The current deliverable is an IPCC-style hexagonal region plot showing PMIP4 climate anomalies across AR6 land regions.

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

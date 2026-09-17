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

- **CVDP model output**: `cvdp_output_by_experiment/` is a **symlink** to the current CVDP release (as of Aug 2026, `/data/oacd/CVDP/v6.1/by_expt_rs1/`; previously `/data/oacd/CVDP/by_experiment/`). The symlink target changes when the CVDP output is refreshed — the code only ever reads through the symlink. One main NetCDF per simulation named `model_experiment.cvdp_data.startyear-endyear.nc`; the only field used is `tas_spatialmean_ann` (annual-mean near-surface temperature, °C). Each model is on its **own native grid**; the experiment and piControl share a model's grid. Auxiliary per-variable outputs (`.siconc.`, `.zos.`, `.monsoon.`, `.tas.indices.`, `.tmp`) sit in the same directory and carry **no tas field** — `model_files` selects the main file by regex (an extra token before the year range means auxiliary), rather than by exclusion list, because each CVDP release adds new auxiliary kinds. Observational/reanalysis files (BEST, ERA5, HadISST…) lack the `_experiment` suffix and so never match. Some grids (e.g. LOVECLIM piControl) carry duplicate longitude values that break interpolation — `load_field` de-duplicates coordinates. Step 2 still checks the file actually carries `tas_spatialmean_ann` and drops+logs any model that doesn't (nothing is dropped under v6.1). Use the plain `midHolocene` CVDP subdirectory, **not** `midHolocene-cal-adj`.
- **Reconstructions**: each period has its own `recons/<period>/` directory; values are already anomalies, and each period's Step 1 cell normalises its files to the shared `compilation, reference, site, Proxy, Latitude, Longitude, Anom` schema.
  - `recons/lig127k/` holds three region-split CSVs from the Otto-Bliesner et al. (2021) supplement (Tables S2–S4), header row 2. Column layouts differ but all share `Compilation reference`, `Latitude`, `Longitude`, `Anom`. Rows split into two compilations by reference string: **Hoffman** (78 pts, global SST) and **Capron** (7 pts, high-latitude).
  - `recons/midHolocene/` holds four compilations. Two scattered-site CSVs with plain headers and unit weights: **Bartlein** (`Bartlein_mat.csv`, 635 pts, pollen MAT, anomaly column `mat_anm_mean`) and **Temp12k** (`temp12k_anom6k_mat.csv`, 335 pts, multiproxy, anomaly column `anom`); both use plain `lat`/`lon` columns and some Bartlein headers carry leading spaces (stripped on load). Plus one **gridded** data assimilation product, **Erb** (Holocene Reconstruction; Erb et al. 2022), flattened to points with `cos(latitude)` weights. **Baseline trap**: the published `recon_tas_mean` is referenced to **3–5 ka**, not the pre-industrial (paper: “The reference period of each reconstruction is 3-5 ka”), so `Erbetal2022_6ka_reconstruction.nc` is *not* comparable with `midHolocene − piControl` — used raw it makes every model look far too warm. Step 1 reads `Erbetal2022_6ka_minus_0-1ka_anom.nc` (5.5–6.5 ka mean minus 0–1 ka mean, so the shared 3–5 ka baseline cancels), built with `ncwa`/`ncdiff` from the 3.7 GB Zenodo original; see `recons/README` for the exact commands. Plus **P2F** (`P2F_recommended_sst_midHolocene.csv`, 124 marine cores, unit weights), the recommended annual-SST records from the supplied workbook; the anomaly is the `SST_Anomaly` column (annual SST minus the **ERSSTv6 1870-1899** mean at the core's grid cell, so already an anomaly, needing no re-referencing). It is an SST compilation benchmarked against model `tas`, exactly as the lig127k and midPliocene marine rows already are. Step 1 reads a **verbatim CSV dump** of the workbook rather than the `.xlsx`, because `my-cli-py` has no `openpyxl`; `recons/README` holds the one-off conversion command. Both the workbooks and the dumps are **gitignored on purpose** (`carpet_diagram/recons/*/P2F*`), so a fresh clone has to be given the workbooks and the dump regenerated before the midHolocene and lgm notebooks will run. One row (`GeoB7112-5`, a Chile-margin core) carries a positive `Latitude` that contradicts its own `ERSST Lat`, so Step 1 takes the hemisphere from the ERSST cell and logs it.
  - `recons/lgm/` holds five compilations. **Bartlein** (`Bartlein_mat_21ka.csv`, NOAA `mat_delta_21ka_ALL_grid_2x2.csv`, 98 cells, anomaly column `mat_anm_mean`, unit weights) is the 21 ka half of the same Bartlein et al. (2011) synthesis as the midHolocene row — verified by the fact that the 6 ka file from the same NOAA directory is identical to `Bartlein_mat.csv` on all 635 shared cells. Two Alaskan cells at 67 °N carry the known anomalously-warm Beringian values (+13.4, +8.1 °C); they are kept, not screened. Then two near-global **gridded** products (Step 1 sets `weight = cos(latitude)` for an area-fair RMSE): **Cleator** (`cleator2020_recon.csv`, ~2214 land grid cells; a long metadata preamble precedes a `#`-prefixed header at line 77 read with `skiprows=77, header=None`; the MAT-anomaly column is `MAT`) and **Osman** (`Osman_LGMR_21ka_SAT_anom_climo.nc`, the LGMR ensemble-mean `sat` field, already a 21 ka − (0–1 ka) anomaly; the 96×144 grid is flattened to ~13824 points, dropping fill values). Note the Osman file was regenerated as an anomaly (`ncdiff` of the 21 ka and 0–1 ka climatologies) — the original `..._SAT_climo.nc` was absolute temperature. And **Annan** (`Annan_etal22.SAT.mean.lgm.nc`, Annan et al. 2022), an ensemble-Kalman-filter product blending proxies with a 19-member PMIP prior; the `SAT.mean` field is already an LGM anomaly (area-weighted global mean −4.46 °C vs the paper's −4.5 ± 0.9 °C), 2°×2°, 16200 points, all finite. Its coordinates are named `latitude`/`longitude`, and its longitude axis is **rolled** (181…359 then 1…179, i.e. non-monotonic) — harmless as used, because Step 1 flattens it to points so each value keeps its own longitude, but it would break anything that tried to interpolate the field. **Cleator is a data assimilation product, not a proxy compilation** (its own preamble: assimilates Bartlein et al. 2011 pollen into averaged PMIP LGM output via 3D-Var; the Reading record has no proxy-only file, and the data are on a regular 2° grid, not at sites). It therefore sits *below* the double rule with the other DA products, and its row is **not independent of the lgm Bartlein row**, which is its own proxy input. The fifth compilation is **P2F** (`P2F_recommended_sst_lgm.csv`, 116 marine cores, unit weights) — the LGM half of the same recommended-SST workbook as the midHolocene row, parsed by the same code.
  - `recons/midPliocene-eoi400/` holds two compilations. **Foley-Dowsett** (`FoleyDowsett2019_cs_mp_sst_data_30k_plus_NOAA.csv`, 39 Uk'37 SST sites, 37 usable). Following Haywood et al. (2020) the anomaly is the **NOAA ERSST5** column — the last, unnamed column (`Unnamed: 14`), i.e. SST minus the NOAA ERSST5 modern climatology, *not* the HadISST-based `Anom` column. The `Supplementary/`, `.zip`, and `.pdf` alongside it are the untouched supplement bundle the CSV was extracted from. And **Tierney** (PlioDA; Tierney et al. 2024), a gridded data assimilation product: `plioDAMainResults.nc` has a `time` axis in **Ma** with slices 0 (modern), 3.25 and 4.75; 3.25 Ma is the mid-Piacenzian that `midPliocene-eoi400` targets, so the anomaly is 3.25 Ma − modern, pre-computed by `ncdiff` into `plioDAMainResults.tas_annual_latePlio_changes.nc` (the file Step 1 reads). Needs no re-referencing — both slices come from the same product. The 180×360 grid flattens to 62640 points (3 southernmost + 3 northernmost latitude rows are all fill); `cos(latitude)` weights; global mean warming +5.3 °C.

## Pipeline

1. Parse + normalise the recon files (Step 1, period-specific) into the shared `compilation, reference, site, Proxy, Latitude, Longitude, Anom` schema (plus an optional `weight` column). Sources vary by period: CSVs with region-split/preamble/plain headers (including the two P2F workbook dumps), and four gridded NetCDFs (Osman, Annan, Erb, Tierney) flattened to points. Note `Bartlein` and `P2F` each appear as a compilation name under both `midHolocene` and `lgm` — rows are keyed by (period, compilation), so they stay distinct.
2. Pair each experiment model with its piControl, compute `<experiment> − piControl` anomaly on the native grid (regrid control if grids ever differ); drop models whose CVDP file lacks the tas field.
3. Bilinearly sample each model's anomaly field at every proxy location — model lon is 0–360°, proxy lon −180–180°, so targets are wrapped and the field is made cyclic in longitude first. A `weight` column is carried through (defaults to 1.0 when Step 1 omits it).
4. Compute per-model weight-weighted RMSE (and bias) against each compilation. Scattered-site compilations use unit weights, so that reduces to the ordinary RMSE; the gridded near-global products (Cleator, Osman, Annan, Erb, Tierney) use `cos(latitude)` weights so their RMSE is area-fair rather than pole-heavy. Weights are set per compilation in Step 1, so a period can mix the two (midHolocene and midPliocene now do).

Outputs land in `carpet_diagram/output/`, one set per period: `recon_points_<period>.csv`, `model_anom_at_recon_<period>.csv`, `rmse_long_<period>.csv`, `rmse_summary_<period>.csv`. The **two point-level files are gitignored** — they hold every site's coordinates and anomaly, P2F included — so a fresh clone has only the aggregate tables and must run the notebooks before `carpet_regional.py` (which reads `model_anom_at_recon_*.csv`) will work. `carpet_figure.py` auto-discovers every `rmse_long_*.csv` and renders all periods/compilations as rows. Row order is set explicitly by `ROW_GROUPS`, a list of groups of `(period, compilation)` keys drawn in order and separated in the figure by a **double rule**: group 1 is the proxy compilations (midHolocene Bartlein/Temp12k, lgm Bartlein, lig127k Hoffman/Capron, midPliocene Foley-Dowsett), group 2 the two **P2F** recommended-SST rows (midHolocene, lgm), group 3 the data assimilation products (midHolocene Erb, lgm Cleator, lgm Osman, lgm Annan, midPliocene Tierney). Three groups means **two** double rules. Any `(period, compilation)` present in the data but absent from `ROW_GROUPS` is appended as a further group and logged, never silently dropped — so adding a period still needs no edit here, but placing its row where you want it does. Running `carpet_figure.py` writes **two** figures. The headline one is `carpet_diagram.png`, the reduced version: a column per CMIP6-era PMIP model (the 22 in `CMIP6_PMIP_MODELS`, transcribed from the PMIP data-holdings table at pcmdi.llnl.gov — `EC-Earth3-Veg`/`EC-Earth3-Veg-LR` are listed there but have no CVDP output) plus `UofT-CCSM-4`, giving 21 columns. `carpet_diagram_all_models.png` is the every-model version (36 columns), which also carries the PMIP4 column so its colours are centred the same way; PMIP3 is left off there because all 15 of its members already appear as their own columns. Which summaries a figure gets is the `labels` argument to `make_figure`, and `keep_only` controls whether the individual columns are restricted to `KEEP_MODELS`. Two **ensemble-mean summary columns** are pinned to the right of the reduced figure behind a vertical double rule, in `COMPOSITE_LABELS` order: **PMIP4** (the 20 CMIP6 archive models) and **PMIP3** (the 15 pre-CMIP6 ones). `UofT-CCSM-4` keeps its own column and is deliberately *not* folded into PMIP4, since it is not in the CMIP6 holdings — change the membership test in `composite_membership` to include it. Each summary is the plain **mean of its members' RMSEs** within a row, *not* the RMSE of their ensemble-mean anomaly field (which would be lower, because model errors partly cancel, and would need the notebooks re-run rather than just this script). Membership varies by row — no pre-CMIP6 model ran midPliocene, so PMIP3 is blank on those two rows. Both summary columns sit at the **left-hand edge** in `COMPOSITE_LABELS` order, ahead of the individual models, which are listed **alphabetically** rather than ranked. Colour is **normalised per row against the PMIP4 ensemble mean**: each cell is its signed deviation from that row's PMIP4 value, scaled by the row's largest absolute deviation, on a diverging map — better than PMIP4 (lower RMSE) is **green**, worse is **red**, and **PMIP4 is white by definition**. `better_worse_cmap` builds that from the Greens and Reds sequential ramps, fading the inner end of each to pure white so the centre is actually white rather than the pale seam the two would otherwise meet at. The scale spans every column, summaries included, so nothing clips when PMIP3 is a row's most extreme value. Both figures now carry PMIP4, so both are centred on it; the row-mean fallback (with colourbar and footnote relabelling themselves) only kicks in if `make_figure` is called with no PMIP4 column. This keeps model ranking legible despite very different absolute RMSE between periods; the printed cell numbers remain absolute RMSE in °C.

## Regional variant (`carpet_regional.py`)

`carpet_regional.py` redraws the same carpet diagram over only the reconstruction points inside a region — either a lat/lon box (default **45–66.33 °N, 30 °W–15 °E**, NW Europe / NE Atlantic, northern edge at the Arctic Circle) or a named **IPCC AR6 reference region** via `--ar6 ABBREV`. AR6 polygons come from `regionmask.defined_regions.ar6.all` — the 58 regions of Iturbide et al. (2020, ESSD 12, 2959–2970) behind the AR6 WGI Atlas, the same source `synthesis_figure/hex_figure.py` uses. Use `.all`, **not** `.land`: combined land-and-ocean regions such as **MED** keep their marine part there, and the marine compilations (P2F, Hoffman, Foley-Dowsett) would otherwise be thrown away. AR6 longitudes run −180–180, matching the proxy convention of the point CSVs; a region crossing the dateline would need more care. **MED (Mediterranean)** is AR6 region 19 and is exactly the rectangle **10 °W–40 °E, 30–45 °N** — its 299 polygon vertices are just a densified rectangle, verified equal to `shapely.geometry.box(-10, 30, 40, 45)` — so `--ar6 MED` and `--lat 30 45 --lon -10 40` select identically; the `--ar6` form is preferred because it carries the provenance and generalises to the non-rectangular regions. It needs nothing re-run upstream **provided the notebooks have been run in this working copy** (that file is gitignored, so it is not in a fresh clone): `output/model_anom_at_recon_<period>.csv` carries every proxy point with its coordinates, reconstructed anomaly, weight and each model's sampled anomaly, so the script just re-does **Step 4** (the same weight-weighted RMSE and bias, with `cos(latitude)` weights still doing the area-fair work for the gridded products) over the boxed rows. Tables are written to `output/region_<name>/rmse_long_<period>.csv` — a **subdirectory on purpose**, so the `rmse_long_*.csv` glob behind the headline figure never picks them up — and drawn by `carpet_figure.make_figure(rmse_dir=...)`. Output is `carpet_diagram_<name>.png` plus `carpet_diagram_<name>_all_models.png`.

**The regional figure prints anomalies, not errors.** Its cell numbers are `model_mean`, the weighted mean of that model's own anomaly over the row's boxed points (so LGM cells read about −5 to −15 °C, not an RMSE), and a **reconstruction column** past a vertical double rule at the right-hand edge carries `recon_mean`, the same weighted mean of the reconstruction — the value the models are being read against. **Colour is still the RMSE** relative to the PMIP4 mean, exactly as in the headline figure, so ranking and numbers answer different questions in the same cell: a green cell whose number is far from the reconstruction column just means every other model is further off. The regional tables therefore carry two columns the notebooks' own `rmse_long_*.csv` do not (`model_mean`, `recon_mean`); `recon_mean` is a row property and repeats down each model column. `make_figure` gained `rmse_dir`, `note`, `value_col`, `value_label` and `recon_col` for all this, and defaults to the old behaviour everywhere else — the headline figure still prints RMSE and has no reconstruction column.

```bash
python carpet_regional.py                                    # the default box
python carpet_regional.py --lat 45 66.33 --lon -30 15 --name euroatlantic
python carpet_regional.py --ar6 MED                          # an IPCC AR6 region
```

`--name` defaults to the AR6 abbreviation in lower case, so `--ar6 MED` writes `output/region_med/` and `carpet_diagram_med.png` / `carpet_diagram_med_all_models.png`. Point counts in MED: lgm Annan 140 / Osman 128 / Cleator 97 / Bartlein 13 / P2F 3, midHolocene Erb 55 / Bartlein 46 / Temp12k 17 / P2F 10, midPliocene Tierney 600 / Foley-Dowsett 2, lig127k Hoffman 1 (Capron drops out again). The lig127k row there rests on a **single** point, so it is noise with a star on it.

Longitude bounds use the **proxy convention (−180–180)**, matching the `Longitude` column of the point CSVs, not the models' 0–360. A compilation with no points in the box is dropped and logged (`lig127k` **Capron** goes, being high-latitude — the default box keeps 10 of the 11 rows). Small samples survive but are noisy and carry the existing `*` flag: in the default box lgm Bartlein has 3 points, lgm P2F 4, lig127k Hoffman 4, midPliocene Foley-Dowsett 3 and midHolocene P2F 7, against midHolocene Bartlein 87 / Temp12k 41 / Erb 40, lgm Annan 88 / Osman 77 / Cleator 36, and midPliocene Tierney 315.

## Outstanding — CSIRO-Mk3L-1-2 piControl not yet restored (as of 7 Aug 2026)

CVDP v6.1 lost three piControl runs that the previous release had: **CSIRO-Mk3L-1-2**, **HadCM3** and **KCM1-2-2**. A missing piControl removes the model from *everything*, because both the carpet diagram and the scatterplots work in `<experiment> − piControl` anomalies.

The one that costs real data is **CSIRO-Mk3L-1-2**: its `1pctCO2`, `abrupt4xCO2` and `midHolocene` runs are all still present in `by_expt_rs1`, but they are orphaned and silently dropped. So the current midHolocene carpet row shows **31 models rather than 32**, and the scatterplots are missing that model entirely. (HadCM3 has no experiment runs in either release, and KCM1-2-2 lost its midHolocene run too, so neither strands anything.)

**The fix is one copy** — the file exists, it just wasn't carried into `by_expt_rs1`:

```bash
cp -n /data/oacd/CVDP/v6.1/MMM_rs3/CSIRO-Mk3L-1-2_piControl.cvdp_data.1-1000.nc \
      /data/oacd/CVDP/v6.1/by_expt_rs1/piControl/
```

Verified before proposing it: the file carries `tas_spatialmean_ann` on the same 56×64 grid as the midHolocene run (implying a +0.12 °C anomaly, matching the +0.13 °C seen before it went missing), and it has the `amm_timeseries_mon` / `nino34` / `pdv_pattern_mon` / `nao_pattern_djf` / `pr_spatialmean_ann` fields the scatterplots need. Writing into `/data/oacd/` is outside the repo, so it needs doing by hand.

**After copying, re-run both chains** (nothing else needs editing):

```bash
cd carpet_diagram && jupyter nbconvert --to notebook --execute --inplace carpet_diagram_midHolocene.ipynb && python carpet_figure.py
cd ../scatterplots && ncl -n compile_tidy_numbers.ncl   # ncl-nco env
python scatter_vs_gmt.py                                 # my-cli-py env
```

Also still available in `MMM_rs3` and not yet restored: `FGOALS-g3_1pctCO2.cvdp_data.370-526.nc`. The other v6.1 losses (CESM2-FV2 1pctCO2/abrupt4xCO2, EC-EARTH-2-2 and KCM1-2-2 midHolocene, GISS-E2-R midPliocene) are not anywhere under `/data/oacd/CVDP/v6.1/` and cannot be recovered this way.

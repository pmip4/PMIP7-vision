# PMIP7 Synthesis Figure — IPCC Hexagonal Region Plot

## Concept

Recreate the IPCC Interactive Atlas hexagonal region figure style
(https://interactive-atlas.ipcc.ch/permalink/VcbzRmwn) for PMIP7 palaeoclimate experiments.

### How the figure works

- One **hexagon per IPCC AR6 land region**, arranged in a schematic geographic layout
- Each hexagon is divided into **6 triangular wedges**: 3 scenarios × 2 variables
  - **Teal wedges** = Mean surface temperature change
  - **Orange wedges** = Mean precipitation change
- **Arrows** within each wedge encode:
  - Direction: upward (increasing) or downward (decreasing)
  - Confidence: filled/large (high), outline/medium, small (low)
  - Attribution: red (high attribution) or blue (low attribution)

## Open Questions

1. **Data source** — Where are the AR6 regional mean values coming from?
   - CVDP output from `extract_all_AR6_regions.ncl`?
   - Pre-computed CSV/NetCDF files?

2. **Scenarios** — What are the 3 experiments/scenarios to compare?
   - e.g., LGM, mid-Holocene, Last Interglacial?
   - Or future SSPs for comparison?

3. **Signal encoding** — Replicate full IPCC arrow+confidence design, or simpler?
   - Option A: Exact replication (arrows, fill by direction, confidence encoding)
   - Option B: Simplified (fill color = magnitude, opacity = confidence)

4. **Region subset** — All ~46 AR6 land regions, or a subset?

5. **Language** — Python (matplotlib) recommended for bespoke geometric layout

## Proposed Implementation Plan

1. Define hexagonal grid coordinates for AR6 regions (axial/offset coordinates)
2. For each hexagon, draw 6 wedge patches (matplotlib `Wedge` or `Polygon`)
3. Color wedges by variable (teal/orange) and fill by change direction
4. Overlay arrows scaled/styled by confidence level
5. Add region labels and legend

## File Structure

```
synthesis_figure/
├── README.md               # This file
├── hex_figure.py           # Main figure script
├── ar6_hex_layout.py       # AR6 region hexagonal grid coordinates
├── data/                   # Input data (regional means per scenario)
└── output/                 # Generated figures
```

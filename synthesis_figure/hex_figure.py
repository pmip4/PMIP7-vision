"""
PMIP7 Synthesis Figure — IPCC-style hexagonal AR6 region plot.

Reads AR6 region geometries directly from the cached shapefile using pyshp
(avoids fiona/geopandas incompatibility). Uses regionmask for spatial masking.

Each hexagon = one AR6 land region (46 total, including Land-Ocean).
8 wedges per hexagon: left half = temperature, right half = precipitation,
each half split into 4 wedges ordered top→bottom: LGM, mid-Holocene, LIG, mid-Pliocene.
Wedge fill colour encodes direction and magnitude of change vs. pre-industrial.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colorbar import ColorbarBase
import matplotlib.cm as cm
import xarray as xr
import shapefile
import shapely.geometry as sg
import regionmask

# ── Paths & constants ─────────────────────────────────────────────────────────
DATA_DIR = '/mnt/c/Users/ucfaccb/local_repos/PMIP7-vision/synthesis_figure/synthesis_data_files'
SHP_PATH  = '/tmp/ar6regions/IPCC-WGI-reference-regions-v4.shp'
SHP_ZIP   = os.path.expanduser('~/.cache/regionmask/v0.13.0/IPCC-WGI-reference-regions-v4.zip')
OUT_DIR  = '/mnt/c/Users/ucfaccb/local_repos/PMIP7-vision/synthesis_figure/output'

SCENARIOS = ['lgm', 'midHolocene', 'lig127k', 'midPliocene-eoi400']
SCENARIO_LABELS = ['LGM', 'mid-Holocene', 'LIG', 'mid-Pliocene']

TAS_VMAX  = 6.0   # °C   — symmetric diverging range
PR_VMAX   = 50.0  # %    — symmetric diverging range

TEMP_CMAP   = cm.RdBu_r
PRECIP_CMAP = cm.BrBG

HEX_R = 1.0  # circumradius of each hexagon

# ── Build regionmask from local shapefile (pyshp) ─────────────────────────────
def _ensure_shapefile():
    if not os.path.exists(SHP_PATH):
        import zipfile
        os.makedirs(os.path.dirname(SHP_PATH), exist_ok=True)
        with zipfile.ZipFile(SHP_ZIP) as z:
            z.extractall(os.path.dirname(SHP_PATH))


def build_ar6_regions():
    """Return (regionmask.Regions, list[abbrev]) for all Land/Land-Ocean regions."""
    _ensure_shapefile()
    sf = shapefile.Reader(SHP_PATH)
    fields = [f[0] for f in sf.fields[1:]]
    geoms, names, abbrevs = [], [], []
    for i, rec in enumerate(sf.records()):
        d = dict(zip(fields, rec))
        if 'Land' in d['Type']:
            geoms.append(sg.shape(sf.shape(i).__geo_interface__))
            names.append(d['Name'])
            abbrevs.append(d['Acronym'])
    regions = regionmask.Regions(
        outlines=geoms,
        names=names,
        abbrevs=abbrevs,
        numbers=list(range(len(geoms))),
    )
    return regions, abbrevs


# ── Load data & compute area-weighted regional means ──────────────────────────
def compute_regional_means(regions, abbrevs):
    """
    Returns results[scenario][abbrev] = {'tas': float_degC, 'pr': float_pct}.
    Values are area-weighted means within each AR6 region, after averaging
    over the model ensemble dimension first.
    """
    results = {s: {a: {} for a in abbrevs} for s in SCENARIOS}

    # cos(lat) weights — same grid for all files
    ds0 = xr.open_dataset(f'{DATA_DIR}/tas_annClim_PMIP4_lgm-piControl.nc')
    weights = np.cos(np.deg2rad(ds0['lat'])).rename('weights')

    for scen in SCENARIOS:
        print(f'  Processing {scen}...')

        # Temperature
        ds = xr.open_dataset(f'{DATA_DIR}/tas_annClim_PMIP4_{scen}-piControl.nc')
        tas = ds['tas'].mean(dim='model')
        mask = regions.mask(tas)
        for idx, abbrev in enumerate(abbrevs):
            val = float(tas.where(mask == idx).weighted(weights).mean(('lat', 'lon')))
            results[scen][abbrev]['tas'] = val

        # Precipitation (percentage change)
        ds = xr.open_dataset(f'{DATA_DIR}/pr_annClim_PMIP4_{scen}-piControl_percentage.nc')
        pr = ds['pr'].mean(dim='model')
        mask = regions.mask(pr)
        for idx, abbrev in enumerate(abbrevs):
            val = float(pr.where(mask == idx).weighted(weights).mean(('lat', 'lon')))
            results[scen][abbrev]['pr'] = val

    return results


# ── Hex geometry ──────────────────────────────────────────────────────────────
def hex_center(col, row, r=HEX_R):
    """Offset hex grid (pointy-top) → Cartesian. Odd rows shift right by half-width."""
    x = col * np.sqrt(3) * r + (row % 2) * np.sqrt(3) * r / 2
    y = row * 1.5 * r
    return x, y


def wedge_polygon(cx, cy, a_start_deg, a_end_deg, r=HEX_R, n_pts=30):
    """Filled triangle-fan wedge from centre to arc on rim."""
    angles = np.linspace(np.radians(a_start_deg), np.radians(a_end_deg), n_pts)
    rim = np.column_stack([cx + r * np.cos(angles), cy + r * np.sin(angles)])
    return np.vstack([[cx, cy], rim])


def hex_outline(cx, cy, r=HEX_R):
    """6 vertices of a pointy-top hexagon."""
    angles = np.radians(90 + np.arange(7) * 60)
    return cx + r * np.cos(angles), cy + r * np.sin(angles)


# Wedge definitions: (variable, scenario_index, angle_start°, angle_end°)
# Angles measured anticlockwise from positive-x axis (standard maths convention).
# Left half (90°→270°)  = temperature, 4 wedges top→bottom = LGM…mid-Pliocene.
# Right half (270°→90°) = precipitation, same order.
WEDGES = [
    ('tas', 0,  90, 135),
    ('tas', 1, 135, 180),
    ('tas', 2, 180, 225),
    ('tas', 3, 225, 270),
    ('pr',  0, 270, 315),
    ('pr',  1, 315, 360),
    ('pr',  2,   0,  45),
    ('pr',  3,  45,  90),
]


def val_to_color(val, var):
    norm = TwoSlopeNorm(
        vmin=(-TAS_VMAX if var == 'tas' else -PR_VMAX),
        vcenter=0,
        vmax=( TAS_VMAX if var == 'tas' else  PR_VMAX),
    )
    cmap = TEMP_CMAP if var == 'tas' else PRECIP_CMAP
    return cmap(norm(np.clip(val, norm.vmin, norm.vmax)))


def draw_hexagon(ax, cx, cy, region_data):
    for var, scen_idx, a0, a1 in WEDGES:
        scen = SCENARIOS[scen_idx]
        val  = region_data[scen].get(var, np.nan)
        color = '#cccccc' if np.isnan(val) else val_to_color(val, var)
        verts = wedge_polygon(cx, cy, a0, a1)
        ax.add_patch(MplPolygon(verts, closed=True,
                                facecolor=color, edgecolor='white', linewidth=0.3))
    # Thin vertical divider separating temperature (left) from precipitation (right)
    ax.plot([cx, cx], [cy - HEX_R, cy + HEX_R], color='#888888',
            linewidth=0.5, zorder=4)
    hx, hy = hex_outline(cx, cy)
    ax.plot(hx, hy, 'k-', linewidth=0.7, zorder=5)


# ── AR6 hex grid layout ───────────────────────────────────────────────────────
# Approximate geographic positions using offset (col, row) coordinates.
# Row 0 = south (Antarctica), row 10 = north.
HEX_LAYOUT = {
    # ── Americas (cols 0–4) ───────────────────────────────────────────────────
    # Polar / Greenland
    'GIC': (2, 10),
    # North America
    'NWN': (1, 9),  'NEN': (3, 9),
    'WNA': (0, 8),  'CNA': (2, 8),  'ENA': (4, 8),
    'NCA': (1, 7),  'CAR': (3, 7),
    'SCA': (1, 6),
    # South America
    'NWS': (0, 6),  'NSA': (1, 5),  'NES': (3, 6),
    'SAM': (1, 4),  'SWS': (0, 4),  'SES': (2, 5),
    'SSA': (1, 3),
    # W. Antarctica (below S. America tip)
    'WAN': (2, 2),

    # ── 2-column gap (cols 5–6) separates Americas from rest ─────────────────

    # ── Europe (cols 7–9) ─────────────────────────────────────────────────────
    'NEU': (8, 10), 'WCE': (7, 9),  'EEU': (9, 9),
    'MED': (7, 8),
    # ── Africa (cols 7–11) ───────────────────────────────────────────────────
    'SAH':  (7, 7), 'WAF':  (7, 6), 'CAF':  (9, 7),
    'NEAF': (9, 6), 'SEAF': (10, 6),
    'WSAF': (8, 5), 'ESAF': (10, 5), 'MDG': (11, 5),
    # E. Antarctica (below Africa)
    'EAN':  (9, 2),

    # ── Russia / Central Asia / E Asia (cols 11–17) ──────────────────────────
    'RAR': (11, 10), 'WSB': (13, 10), 'ESB': (15, 10), 'RFE': (17, 10),
    'WCA': (11, 9),  'ECA': (13, 9),  'TIB': (15, 9),  'EAS': (17, 9),
    # ── Middle East / South Asia / SE Asia ───────────────────────────────────
    'ARP': (11, 8), 'SAS': (13, 8),  'SEA': (16, 8),

    # ── Australasia — compact cluster (cols 15–17, rows 5–7) ─────────────────
    'NAU': (15, 7), 'EAU': (16, 7),   # N and E Australia, top row
    'CAU': (15, 6), 'NZ':  (17, 6),   # Central Aus left, NZ right
    'SAU': (15, 5),                    # S Australia below CAU
}


# ── Legend hex ────────────────────────────────────────────────────────────────
def draw_legend_inset(fig, rect=(0.02, 0.02, 0.13, 0.17)):
    """
    Draw a self-contained legend hexagon as a figure inset.
    rect = [left, bottom, width, height] in figure-fraction coordinates.
    Scenario labels sit inside the wedges; variable labels are below the hex.
    """
    ax = fig.add_axes(rect)
    ax.set_aspect('equal')
    ax.axis('off')
    r, cx, cy = 1.0, 0.0, 0.2   # cy shifted up slightly to leave room below

    # Wedge fills
    for var, scen_idx, a0, a1 in WEDGES:
        color = TEMP_CMAP(0.80) if var == 'tas' else PRECIP_CMAP(0.80)
        verts = wedge_polygon(cx, cy, a0, a1, r=r)
        ax.add_patch(MplPolygon(verts, closed=True,
                                facecolor=color, edgecolor='white', linewidth=0.5, alpha=0.7))

    # Vertical temp/precip divider
    ax.plot([cx, cx], [cy - r, cy + r], color='#666666', linewidth=0.8, zorder=4)

    # Hex outline
    hx, hy = hex_outline(cx, cy, r=r)
    ax.plot(hx, hy, 'k-', linewidth=0.9)

    # Scenario labels inside each wedge
    label_r = r * 0.65
    for var, scen_idx, a0, a1 in WEDGES:
        a_mid = np.radians((a0 + a1) / 2)
        lx = cx + label_r * np.cos(a_mid)
        ly = cy + label_r * np.sin(a_mid)
        ax.text(lx, ly, SCENARIO_LABELS[scen_idx],
                ha='center', va='center', fontsize=6.5,
                path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])

    # Variable labels below the hex
    ax.text(cx - r * 0.5, cy - r * 1.25, 'Temperature', ha='center', va='top',
            fontsize=8, fontweight='bold', color=TEMP_CMAP(0.85))
    ax.text(cx + r * 0.5, cy - r * 1.25, 'Precipitation', ha='center', va='top',
            fontsize=8, fontweight='bold', color=PRECIP_CMAP(0.15))

    ax.set_xlim(-r * 1.6, r * 1.6)
    ax.set_ylim(-r * 1.6, r * 1.6)
    ax.set_title('How to read', fontsize=8, pad=3)


# ── Main ──────────────────────────────────────────────────────────────────────
def make_figure():
    print('Building AR6 region masks...')
    regions, abbrevs = build_ar6_regions()
    print(f'  {len(abbrevs)} land regions: {abbrevs}')

    print('Computing regional means...')
    data = compute_regional_means(regions, abbrevs)

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(22, 14))
    ax = fig.add_axes([0.02, 0.22, 0.96, 0.76])
    ax.set_aspect('equal')
    ax.axis('off')

    label_opts = dict(
        ha='center', va='center', fontsize=5.5, fontweight='bold', zorder=6,
        path_effects=[pe.withStroke(linewidth=1.5, foreground='white')],
    )

    missing = []
    for abbrev in abbrevs:
        if abbrev not in HEX_LAYOUT:
            missing.append(abbrev)
            continue
        col, row = HEX_LAYOUT[abbrev]
        cx, cy   = hex_center(col, row)
        draw_hexagon(ax, cx, cy, {s: data[s][abbrev] for s in SCENARIOS})
        ax.text(cx, cy, abbrev, color='k', **label_opts)
    if missing:
        print(f'Warning: no layout position for: {missing}')

    # Auto-scale axes with a little padding
    all_x = [hex_center(*HEX_LAYOUT[a])[0] for a in abbrevs if a in HEX_LAYOUT]
    all_y = [hex_center(*HEX_LAYOUT[a])[1] for a in abbrevs if a in HEX_LAYOUT]
    pad = 1.5 * HEX_R
    ax.set_xlim(min(all_x) - pad, max(all_x) + pad)
    ax.set_ylim(min(all_y) - pad, max(all_y) + pad)

    # ── Legend inset (bottom-left, fully below main axes) ─────────────────
    draw_legend_inset(fig, rect=(0.02, 0.02, 0.13, 0.17))

    # ── Colorbars (centred at bottom) ──────────────────────────────────────
    ax_tas = fig.add_axes([0.30, 0.08, 0.18, 0.020])
    ColorbarBase(ax_tas, cmap=TEMP_CMAP,
                 norm=TwoSlopeNorm(vmin=-TAS_VMAX, vcenter=0, vmax=TAS_VMAX),
                 orientation='horizontal',
                 label='Temperature change (°C)')
    ax_tas.tick_params(labelsize=8)

    ax_pr = fig.add_axes([0.55, 0.08, 0.18, 0.020])
    ColorbarBase(ax_pr, cmap=PRECIP_CMAP,
                 norm=TwoSlopeNorm(vmin=-PR_VMAX, vcenter=0, vmax=PR_VMAX),
                 orientation='horizontal',
                 label='Precipitation change (%)')
    ax_pr.tick_params(labelsize=8)

    ax.set_title('PMIP4 Regional Climate Changes relative to pre-industrial control',
                 fontsize=13, pad=8)

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = f'{OUT_DIR}/synthesis_hexfig.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f'Saved → {out_path}')


if __name__ == '__main__':
    make_figure()

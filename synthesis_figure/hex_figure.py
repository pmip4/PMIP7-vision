"""
PMIP7 Synthesis Figure — IPCC-style hexagonal AR6 region plot.

Styled after the IPCC WGI Interactive Atlas regional-synthesis hexagon figure:
tightly packed hexagons forming continent-shaped clusters.

Each hexagon = one AR6 land region (46 total, including Land-Ocean).
8 wedges per hexagon: left half = temperature, right half = precipitation.
Scenarios are mirrored so each pair of wedges meets at the vertical divider,
ordered top→bottom: mid-Pliocene, LIG, mid-Holocene, LGM.
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
import shapely.geometry as sg
import regionmask

# ── Paths & constants ─────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, 'synthesis_data_files')
OUT_DIR  = os.path.join(SCRIPT_DIR, 'output')

SCENARIOS = ['lgm', 'midHolocene', 'lig127k', 'midPliocene-eoi400']
SCENARIO_LABELS = ['LGM', 'mid-Holocene', 'LIG', 'mid-Pliocene']

TAS_VMAX  = 6.0   # °C   — symmetric diverging range
PR_VMAX   = 50.0  # %    — symmetric diverging range

TEMP_CMAP   = cm.RdBu_r
PRECIP_CMAP = cm.BrBG

HEX_R = 1.0  # circumradius of each hexagon

OCEAN_COLOR  = '#f2f7fb'
MISSING_COLOR = '#d8d8d8'


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
        for idx, abbrev in zip(regions.numbers, abbrevs):
            val = float(tas.where(mask == idx).weighted(weights).mean(('lat', 'lon')))
            results[scen][abbrev]['tas'] = val

        # Precipitation (percentage change)
        ds = xr.open_dataset(f'{DATA_DIR}/pr_annClim_PMIP4_{scen}-piControl_percentage.nc')
        pr = ds['pr'].mean(dim='model')
        mask = regions.mask(pr)
        for idx, abbrev in zip(regions.numbers, abbrevs):
            val = float(pr.where(mask == idx).weighted(weights).mean(('lat', 'lon')))
            results[scen][abbrev]['pr'] = val

    return results


# ── Hex geometry ──────────────────────────────────────────────────────────────
def hex_center(col, row, r=HEX_R):
    """Offset hex grid (pointy-top) → Cartesian. Odd rows shift right by half-width."""
    x = col * np.sqrt(3) * r + (row % 2) * np.sqrt(3) * r / 2
    y = row * 1.5 * r
    return x, y


def hex_polygon(cx, cy, r=HEX_R):
    angles = np.radians(90 + np.arange(6) * 60)
    return sg.Polygon(zip(cx + r * np.cos(angles), cy + r * np.sin(angles)))


def wedge_polygon(cx, cy, a_start_deg, a_end_deg, r=HEX_R, n_pts=30):
    """Angular sector clipped to the hexagon, so wedge rims follow the hex edges."""
    angles = np.linspace(np.radians(a_start_deg), np.radians(a_end_deg), n_pts)
    big = 2.5 * r
    fan = sg.Polygon([(cx, cy),
                      *zip(cx + big * np.cos(angles), cy + big * np.sin(angles))])
    clipped = fan.intersection(hex_polygon(cx, cy, r))
    return np.array(clipped.exterior.coords)


def hex_outline(cx, cy, r=HEX_R):
    """6 vertices of a pointy-top hexagon."""
    angles = np.radians(90 + np.arange(7) * 60)
    return cx + r * np.cos(angles), cy + r * np.sin(angles)


# Wedge definitions: (variable, scenario_index, angle_start°, angle_end°)
# Angles measured anticlockwise from positive-x axis (standard maths convention).
# Left half (90°→270°)  = temperature, right half (270°→90°) = precipitation.
# Scenario pairs mirror about the vertical divider, ordered top→bottom:
# mid-Pliocene, LIG, mid-Holocene, LGM.
WEDGES = [
    ('tas', 3,  90, 135),
    ('tas', 2, 135, 180),
    ('tas', 1, 180, 225),
    ('tas', 0, 225, 270),
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
        color = MISSING_COLOR if np.isnan(val) else val_to_color(val, var)
        verts = wedge_polygon(cx, cy, a0, a1)
        ax.add_patch(MplPolygon(verts, closed=True,
                                facecolor=color, edgecolor='white', linewidth=0.4))
    # Thin vertical divider separating temperature (left) from precipitation (right)
    ax.plot([cx, cx], [cy - HEX_R, cy + HEX_R], color='white',
            linewidth=1.0, zorder=4)
    hx, hy = hex_outline(cx, cy)
    ax.plot(hx, hy, color='#4d4d4d', linewidth=0.9, zorder=5,
            solid_capstyle='round')


# ── AR6 hex grid layout ───────────────────────────────────────────────────────
# Offset (col, row) coordinates; row 0 = south, higher rows = north.
# Hexagons pack tightly into continent-shaped clusters, mimicking the
# IPCC Interactive Atlas regional-synthesis figure.
HEX_LAYOUT = {
    # ── North & Central America ───────────────────────────────────────────────
    'NWN': (1, 10), 'NEN': (2, 10), 'GIC': (3, 10),
    'WNA': (0, 9),  'CNA': (1, 9),  'ENA': (2, 9),
    'NCA': (1, 8),                  'CAR': (3, 8),
    'SCA': (1, 7),

    # ── South America ─────────────────────────────────────────────────────────
    'NWS': (2, 6),  'NSA': (3, 6),
    'SAM': (2, 5),  'NES': (3, 5),
    'SWS': (2, 4),  'SES': (3, 4),
    'SSA': (2, 3),

    # ── Antarctica ────────────────────────────────────────────────────────────
    'WAN': (6, 0),  'EAN': (7, 0),

    # ── Europe ────────────────────────────────────────────────────────────────
    'NEU': (6, 10),
    'WCE': (5, 9),  'EEU': (6, 9),
    'MED': (6, 8),

    # ── Africa ────────────────────────────────────────────────────────────────
    'SAH': (5, 7),
    'WAF': (5, 6),  'CAF': (6, 6),  'NEAF': (7, 6),
    'SEAF': (6, 5),
    'WSAF': (6, 4), 'ESAF': (7, 4), 'MDG': (8, 4),

    # ── Asia ──────────────────────────────────────────────────────────────────
    'RAR': (8, 10), 'WSB': (9, 10), 'ESB': (10, 10), 'RFE': (11, 10),
    'WCA': (8, 9),  'ECA': (9, 9),  'TIB': (10, 9),  'EAS': (11, 9),
    'ARP': (8, 8),  'SAS': (10, 8), 'SEA': (12, 8),

    # ── Australasia ───────────────────────────────────────────────────────────
    'NAU': (12, 5),
    'CAU': (12, 4), 'EAU': (13, 4),
    'SAU': (12, 3),                 'NZ': (14, 3),
}

# Continent labels: name → (x, y) in plot coordinates
SQRT3 = np.sqrt(3)
CONTINENT_LABELS = {
    'NORTH AMERICA': (1.8 * SQRT3, 16.9),
    'SOUTH AMERICA': (2.75 * SQRT3, 2.7),
    'ANTARCTICA':    (7.0 * SQRT3, -1.8),
    'EUROPE':        (6.0 * SQRT3, 16.9),
    'AFRICA':        (6.75 * SQRT3, 4.2),
    'ASIA':          (10.0 * SQRT3, 16.9),
    'AUSTRALASIA':   (13.0 * SQRT3, 2.7),
}


# ── Legend hex ────────────────────────────────────────────────────────────────
def draw_legend_inset(fig, rect=(0.015, 0.015, 0.15, 0.20)):
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
    ax.plot([cx, cx], [cy - r, cy + r], color='white', linewidth=1.0, zorder=4)

    # Hex outline
    hx, hy = hex_outline(cx, cy, r=r)
    ax.plot(hx, hy, color='#4d4d4d', linewidth=1.0)

    # One label per scenario pair, centred on the vertical divider at the
    # height of that pair of wedges (top→bottom: mid-Pliocene … LGM).
    pair_mid_angles = [112.5, 157.5, 202.5, 247.5]  # left-wedge mid-angles, top→bottom
    for scen_idx, a_mid in zip([3, 2, 1, 0], pair_mid_angles):
        ly = cy + r * 0.78 * np.sin(np.radians(a_mid))
        ax.text(cx, ly, SCENARIO_LABELS[scen_idx],
                ha='center', va='center', fontsize=6.5, zorder=6,
                path_effects=[pe.withStroke(linewidth=1.8, foreground='white')])

    # Variable labels below the hex
    ax.text(cx - r * 0.55, cy - r * 1.25, 'Temperature', ha='center', va='top',
            fontsize=8.5, fontweight='bold', color=TEMP_CMAP(0.85))
    ax.text(cx + r * 0.55, cy - r * 1.25, 'Precipitation', ha='center', va='top',
            fontsize=8.5, fontweight='bold', color=PRECIP_CMAP(0.15))

    ax.set_xlim(-r * 1.7, r * 1.7)
    ax.set_ylim(-r * 1.7, r * 1.7)
    ax.set_title('How to read', fontsize=9, pad=3)


# ── Main ──────────────────────────────────────────────────────────────────────
def make_figure():
    print('Loading AR6 land regions from regionmask...')
    regions = regionmask.defined_regions.ar6.land
    abbrevs = list(regions.abbrevs)
    print(f'  {len(abbrevs)} land regions')

    print('Computing regional means...')
    data = compute_regional_means(regions, abbrevs)

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(20, 14))
    fig.patch.set_facecolor(OCEAN_COLOR)
    ax = fig.add_axes([0.01, 0.16, 0.98, 0.75])
    ax.set_aspect('equal')
    ax.axis('off')

    label_opts = dict(
        ha='center', va='center', fontsize=7, fontweight='bold', zorder=6,
        color='#1a1a1a',
        path_effects=[pe.withStroke(linewidth=1.8, foreground='white')],
    )

    missing = []
    for abbrev in abbrevs:
        if abbrev not in HEX_LAYOUT:
            missing.append(abbrev)
            continue
        col, row = HEX_LAYOUT[abbrev]
        cx, cy   = hex_center(col, row)
        draw_hexagon(ax, cx, cy, {s: data[s][abbrev] for s in SCENARIOS})
        ax.text(cx, cy, abbrev, **label_opts)
    if missing:
        print(f'Warning: no layout position for: {missing}')

    # Continent labels
    for name, (lx, ly) in CONTINENT_LABELS.items():
        ax.text(lx, ly, name, ha='center', va='center',
                fontsize=11, fontweight='bold', color='#8aa0b4',
                fontstyle='normal', zorder=2)

    # Auto-scale axes with a little padding
    all_x = [hex_center(*HEX_LAYOUT[a])[0] for a in abbrevs if a in HEX_LAYOUT]
    all_y = [hex_center(*HEX_LAYOUT[a])[1] for a in abbrevs if a in HEX_LAYOUT]
    pad = 2.0 * HEX_R
    ax.set_xlim(min(all_x) - pad, max(all_x) + pad)
    ax.set_ylim(min(all_y) - 2.5 * HEX_R, max(all_y) + 2.5 * HEX_R)

    # ── Legend inset (bottom-left) ─────────────────────────────────────────
    draw_legend_inset(fig)

    # ── Colorbars (centred at bottom) ──────────────────────────────────────
    ax_tas = fig.add_axes([0.30, 0.075, 0.18, 0.018])
    cb_tas = ColorbarBase(ax_tas, cmap=TEMP_CMAP,
                          norm=TwoSlopeNorm(vmin=-TAS_VMAX, vcenter=0, vmax=TAS_VMAX),
                          orientation='horizontal', extend='both',
                          label='Temperature change (°C)')
    ax_tas.tick_params(labelsize=8, length=2)
    cb_tas.outline.set_linewidth(0.5)

    ax_pr = fig.add_axes([0.55, 0.075, 0.18, 0.018])
    cb_pr = ColorbarBase(ax_pr, cmap=PRECIP_CMAP,
                         norm=TwoSlopeNorm(vmin=-PR_VMAX, vcenter=0, vmax=PR_VMAX),
                         orientation='horizontal', extend='both',
                         label='Precipitation change (%)')
    ax_pr.tick_params(labelsize=8, length=2)
    cb_pr.outline.set_linewidth(0.5)

    fig.suptitle('PMIP4 Regional Climate Changes relative to pre-industrial',
                 fontsize=16, fontweight='bold', y=0.97, color='#222222')
    fig.text(0.5, 0.935,
             'Ensemble-mean annual temperature and precipitation change over IPCC AR6 land regions',
             ha='center', fontsize=10.5, color='#555555')

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = f'{OUT_DIR}/synthesis_hexfig.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(f'Saved → {out_path}')


if __name__ == '__main__':
    make_figure()

"""Quick-look emission maps (matplotlib). Publication cartography lives in analysis_R/."""
from __future__ import annotations

import numpy as np


def plot_field(field, lon, lat, *, title, units, outpath, cmap="YlOrRd",
               log=False, mask=None):
    """Render a 2-D field over lon/lat to ``outpath`` (PNG).

    Uses a plain geographic axes with a colorbar; adds coastlines/borders if
    cartopy is available, otherwise a labelled lat/lon plot. ``log`` uses a log
    color scale (good for emissions spanning orders of magnitude).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    field = np.asarray(field, dtype=float)
    if mask is not None:
        field = np.where(mask, field, np.nan)

    norm = None
    if log:
        pos = field[np.isfinite(field) & (field > 0)]
        vmin = pos.min() if pos.size else 1e-12
        norm = LogNorm(vmin=max(vmin, 1e-12), vmax=np.nanmax(field))

    extent = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        fig = plt.figure(figsize=(7, 5))
        ax = plt.axes(projection=ccrs.PlateCarree())
        im = ax.imshow(field, origin="lower", extent=extent, cmap=cmap, norm=norm,
                       transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
        ax.add_feature(cfeature.BORDERS, linewidth=0.4)
        ax.add_feature(cfeature.STATES if hasattr(cfeature, "STATES") else cfeature.BORDERS,
                       linewidth=0.3)
        gl = ax.gridlines(draw_labels=True, linewidth=0.2)
        gl.top_labels = gl.right_labels = False
    except Exception:
        fig, ax = plt.subplots(figsize=(7, 5))
        im = ax.imshow(field, origin="lower", extent=extent, cmap=cmap, norm=norm,
                       aspect="auto")
        ax.set_xlabel("longitude"); ax.set_ylabel("latitude")

    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.03)
    cb.set_label(units)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=140)
    plt.close(fig)
    return outpath

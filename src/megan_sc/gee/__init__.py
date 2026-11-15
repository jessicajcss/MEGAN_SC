"""Google Earth Engine data-acquisition layer for MEGAN-SC.

Pulls the meteorology (ERA5-Land) and leaf-area (MODIS) MEGAN needs, for a
user-defined domain and period, following Earth-Engine best practices for large
imagery: server-side computation, filter-early/clip-late, batch-by-month exports
with explicit crs/scale/maxPixels/tileScale, and idempotent skipping of months
already produced.

These modules are authored to run under the user's EE project (e.g.
``ee-yourname``); they are not executed inside the build sandbox. Call
``init(project)`` once before using the collection builders.
"""
from __future__ import annotations


def initialize(project=None, *, auth_mode="notebook", service_account_key=None,
               service_account_email=None, force=False):
    """Initialise Earth Engine, choosing an auth method that works off-localhost.

    The default browser/localhost token flow (``ee.Authenticate()``) often fails on lab
    machines (no localhost redirect, firewall, headless). This helper avoids that:

    * **Already authenticated once?** It just calls ``ee.Initialize`` — no prompt. EE caches
      a refresh token at ``~/.config/earthengine/credentials`` after the first success, so you
      authenticate *once per machine*, not once per run.
    * ``auth_mode="notebook"`` (default) — ``ee.Authenticate(auth_mode='notebook')`` prints a
      URL and asks you to paste back a verification code. **No localhost redirect needed** —
      this is the fix when the pop-up/token flow won't complete.
    * ``auth_mode="localhost"`` — the original browser flow (only if it works for you).
    * ``auth_mode="gcloud"`` — reuse Google Cloud SDK application-default credentials.
    * ``service_account_key=<path.json>`` — **fully non-interactive** (best for automation and
      reproducible runs): a Google Cloud service-account key. The account must be granted Earth
      Engine access (register it at https://code.earthengine.google.com or give it the
      *Earth Engine Resource Viewer/Writer* role in the project).

    Returns the ``ee`` module, ready to use.
    """
    import ee

    # 1) Service account — no interaction at all.
    if service_account_key:
        email = service_account_email
        if email is None:
            import json
            with open(service_account_key) as f:
                email = json.load(f).get("client_email")
        creds = ee.ServiceAccountCredentials(email, key_file=service_account_key)
        ee.Initialize(creds, project=project)
        return ee

    # 2) Reuse cached credentials if we have them (no prompt).
    if not force:
        try:
            ee.Initialize(project=project)
            return ee
        except Exception:
            pass

    # 3) Interactive auth via the chosen mode, then initialise.
    ee.Authenticate(auth_mode=auth_mode)
    ee.Initialize(project=project)
    return ee


def initialize_from_config(config):
    """Initialise EE from a loaded config's ``earth_engine`` block."""
    ee_cfg = config.get("earth_engine", {})
    return initialize(
        project=ee_cfg.get("project"),
        auth_mode=ee_cfg.get("auth_mode", "notebook"),
        service_account_key=ee_cfg.get("service_account_key"),
        service_account_email=ee_cfg.get("service_account_email"),
    )


def init(project: str):
    """Backward-compatible alias for :func:`initialize` (notebook auth by default)."""
    return initialize(project=project)


def region_from_bounds(bounds):
    """ee.Geometry.Rectangle from (lon_min, lat_min, lon_max, lat_max)."""
    import ee
    lon_min, lat_min, lon_max, lat_max = bounds
    return ee.Geometry.Rectangle([lon_min, lat_min, lon_max, lat_max], proj="EPSG:4326", geodesic=False)

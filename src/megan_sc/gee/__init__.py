"""Google Earth Engine data-acquisition layer for MEGAN-SC.

Pulls the meteorology (ERA5-Land) and leaf-area (MODIS) MEGAN needs, for a
user-defined domain and period, following Earth-Engine best practices for large
imagery: server-side computation, filter-early/clip-late, batch-by-month exports
with explicit crs/scale/maxPixels/tileScale, and idempotent skipping of months
already produced.

These modules run under the user's Earth Engine project (a Google Cloud project
registered for Earth Engine; tutorial, chapter 3). Call ``initialize(project)`` once
before using the collection builders.
"""
from __future__ import annotations


def initialize(project=None, *, auth_mode=None, service_account_key=None,
               service_account_email=None, force=False):
    """Initialise Earth Engine for ``project`` (a Google Cloud project registered for Earth Engine).

    The order follows the Earth Engine guide "Authentication and Initialization"
    (https://developers.google.com/earth-engine/guides/auth):

    * ``service_account_key=<path.json>``: a service-account key, for unattended runs
      (``ee.ServiceAccountCredentials`` then ``ee.Initialize(credentials, project=...)``). The service account needs
      access to Earth Engine in the project (Earth Engine service-account guide).
    * Otherwise the credentials stored by an earlier ``ee.Authenticate()`` or ``earthengine authenticate`` (file
      ``~/.config/earthengine/credentials``) are tried first, so a machine authenticates once.
    * If that fails (or ``force=True``), ``ee.Authenticate(auth_mode=auth_mode)`` runs. ``auth_mode=None`` lets the
      Earth Engine client choose its default for the environment: ``colab`` in Google Colab, ``notebook`` in other
      Jupyter notebooks, ``localhost`` when a browser is available and gcloud is not installed, ``gcloud`` otherwise.
      Pass ``"notebook"`` on a remote machine without a browser.

    Returns the ``ee`` module, ready to use.
    """
    import ee

    if not project or project == "your-ee-project":
        raise ValueError("Earth Engine needs the id of a Google Cloud project registered for Earth Engine "
                         "(tutorial, chapter 3); got %r" % (project,))
    if service_account_key:
        email = service_account_email
        if email is None:
            import json
            with open(service_account_key) as f:
                email = json.load(f).get("client_email")
        creds = ee.ServiceAccountCredentials(email, key_file=service_account_key)
        ee.Initialize(creds, project=project)
        return ee
    if not force:
        try:
            ee.Initialize(project=project)
            return ee
        except Exception:
            pass
    if auth_mode is None:
        ee.Authenticate(force=force)
    else:
        ee.Authenticate(auth_mode=auth_mode, force=force)
    ee.Initialize(project=project)
    return ee


def initialize_from_config(config):
    """Initialise EE from a loaded config's ``earth_engine`` block."""
    ee_cfg = config.get("earth_engine", {})
    return initialize(
        project=ee_cfg.get("project"),
        auth_mode=ee_cfg.get("auth_mode"),
        service_account_key=ee_cfg.get("service_account_key"),
        service_account_email=ee_cfg.get("service_account_email"),
    )


def init(project: str):
    """Backward-compatible alias for :func:`initialize`."""
    return initialize(project=project)


def region_from_bounds(bounds):
    """ee.Geometry.Rectangle from (lon_min, lat_min, lon_max, lat_max)."""
    import ee
    lon_min, lat_min, lon_max, lat_max = bounds
    return ee.Geometry.Rectangle([lon_min, lat_min, lon_max, lat_max], proj="EPSG:4326", geodesic=False)

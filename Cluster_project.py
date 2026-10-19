#!/usr/bin/env python3
"""Identify open-cluster members and estimate their ages.

The script queries Gaia DR3 ``gaiadr3.gaia_source_lite`` for the Pleiades, Praesepe,
and M67, clusters stars in (pmra, pmdec, parallax), writes member catalogues,
plots diagnostic and Hertzsprung-Russell diagrams, and fits Gaia isochrones for
all clusters. Run ``python Cluster_project.py --help`` for usage.

Isochrone fitting can use a supplied CSV grid with the columns ``log_age``,
``feh``, ``abs_g``, and ``bp_rp`` or, when no CSV is supplied, a built-in MIST
Gaia grid generated through the installed ``isochrones`` package.
"""

from __future__ import annotations

import argparse
import json
import logging
import ssl
import tarfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.request import urlretrieve

import numpy as np


LOG = logging.getLogger("gaia_clusters")
GAIA_ARCHIVE_TABLE = "gaiadr3.gaia_source_lite"
QUERY_VERSION = "4"
PARALLAX_OFFSET_MAS = 0.017
MIN_MEMBER_SCORE = 0.05
MAX_OBSERVED_G_MAG = 18.0
REQUIRED_COLUMNS = {
    "source_id",
    "ra",
    "dec",
    "parallax",
    "parallax_error",
    "pmra",
    "pmra_error",
    "pmdec",
    "pmdec_error",
    "phot_g_mean_mag",
    "phot_bp_mean_mag",
    "phot_rp_mean_mag",
}
OPTIONAL_COLUMNS = {
    "ruwe",
    "parallax_over_error",
    "visibility_periods_used",
    "phot_bp_rp_excess_factor",
    "phot_g_mean_flux",
    "phot_g_mean_flux_error",
}


@dataclass(frozen=True)
class ClusterConfig:
    """Fixed starting values; all scientific cuts remain CLI-configurable."""

    name: str
    ra_deg: float
    dec_deg: float
    radius_deg: float
    expected_pmra: float
    expected_pmdec: float
    expected_parallax: float
    av_prior: float
    feh_prior: float
    log_age_min: float
    log_age_max: float
    max_g_mag: float | None = None


CLUSTERS: dict[str, ClusterConfig] = {
    "pleiades": ClusterConfig(
        "Pleiades", 56.75, 24.1167, 2.0, 19.15, -45.72, 7.35, 0.12, 0.00, 7.0, 9.0, 18.0
    ),
    "praesepe": ClusterConfig(
        "Praesepe", 130.05, 19.67, 2.0, -35.70, -12.80, 5.40, 0.08, 0.12, 8.2, 9.6
    ),
    "m67": ClusterConfig(
        "M67", 132.825, 11.800, 2.0, -10.90, -2.90, 1.14, 0.13, 0.03, 9.0, 10.2
    ),
}


def require_runtime() -> dict[str, Any]:
    """Load optional scientific dependencies only when executing the pipeline."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd
        from astroquery.gaia import Gaia
        from astropy.coordinates import SkyCoord
        from astropy.table import Table
        import astropy.units as units
        from scipy.spatial import cKDTree
        from sklearn.cluster import DBSCAN
        from sklearn.preprocessing import RobustScaler
    except ImportError as error:
        raise RuntimeError(
            "Missing a required package. Install: "
            "pip install astroquery astropy numpy pandas scipy scikit-learn matplotlib"
        ) from error
    return {
        "pd": pd,
        "plt": plt,
        "Gaia": Gaia,
        "SkyCoord": SkyCoord,
        "units": units,
        "Table": Table,
        "cKDTree": cKDTree,
        "DBSCAN": DBSCAN,
        "RobustScaler": RobustScaler,
    }


def cache_paths(cache_dir: Path, cluster: ClusterConfig) -> tuple[Path, Path]:
    stem = cluster.name.lower().replace(" ", "_")
    return cache_dir / f"{stem}_raw.ecsv", cache_dir / f"{stem}_query.json"


def query_gaia_archive_cone_search(
    cluster: ClusterConfig,
    columns: list[str],
    runtime: dict[str, Any],
    *,
    allow_insecure_ssl: bool = False,
) -> Any:
    gaia = runtime["Gaia"]
    gaia.ROW_LIMIT = -1
    if allow_insecure_ssl:
        previous_context = ssl._create_default_https_context
        ssl._create_default_https_context = ssl._create_unverified_context
    try:
        job = gaia.cone_search_async(
            runtime["SkyCoord"](cluster.ra_deg, cluster.dec_deg, unit="deg"),
            radius=cluster.radius_deg * runtime["units"].deg,
            table_name=GAIA_ARCHIVE_TABLE,
            columns=columns,
        )
        return job.get_results()
    finally:
        if allow_insecure_ssl:
            ssl._create_default_https_context = previous_context


def fetch_catalogue(
    cluster: ClusterConfig,
    cache_dir: Path,
    refresh: bool,
    offline: bool,
    runtime: dict[str, Any],
) -> tuple[Any, dict[str, Any]]:
    """Return cached data or query Gaia DR3 gaia_source_lite."""
    table_path, metadata_path = cache_paths(cache_dir, cluster)
    table_class = runtime["Table"]
    expected_cluster = asdict(cluster)
    if table_path.exists() and metadata_path.exists() and not refresh:
        cached_metadata = json.loads(metadata_path.read_text())
        if (
            cached_metadata.get("query_version") == QUERY_VERSION
            and cached_metadata.get("catalog") == GAIA_ARCHIVE_TABLE
            and cached_metadata.get("cluster") == expected_cluster
        ):
            LOG.info("Using cached catalogue for %s", cluster.name)
            return table_class.read(table_path, format="ascii.ecsv"), cached_metadata
        if offline:
            raise FileNotFoundError(
                f"Cached data for {cluster.name} do not match the current configuration. "
                "Re-run without --offline to refresh the Gaia Archive catalogue."
            )
        LOG.info("Refreshing cached catalogue for %s because the configuration changed", cluster.name)
    elif offline:
        raise FileNotFoundError(
            f"No cached data for {cluster.name}. Re-run without --offline to query the Gaia Archive."
        )

    columns = [
        "source_id", "ra", "dec", "parallax", "parallax_error", "pmra", "pmra_error",
        "pmdec", "pmdec_error", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag",
    ]
    LOG.info("Querying %s (%s degree cone) via Gaia Archive %s", cluster.name, cluster.radius_deg, GAIA_ARCHIVE_TABLE)
    try:
        table = query_gaia_archive_cone_search(cluster, columns, runtime)
    except Exception as error:
        ssl_error = isinstance(error, ssl.SSLError) or "CERTIFICATE_VERIFY_FAILED" in str(error)
        if not ssl_error:
            raise RuntimeError(f"Gaia Archive query failed for {cluster.name}: {error}") from error
        LOG.warning(
            "Gaia Archive SSL verification failed for %s; retrying with certificate verification disabled.",
            cluster.name,
        )
        try:
            table = query_gaia_archive_cone_search(cluster, columns, runtime, allow_insecure_ssl=True)
        except Exception as retry_error:
            raise RuntimeError(f"Gaia Archive query failed for {cluster.name}: {retry_error}") from retry_error
    if len(table) == 0:
        raise RuntimeError("Gaia Archive returned no rows")

    metadata = {
        "cluster": asdict(cluster),
        "source": "ESA Gaia Archive",
        "catalog": GAIA_ARCHIVE_TABLE,
        "query": {"ra_deg": cluster.ra_deg, "dec_deg": cluster.dec_deg, "radius_deg": cluster.radius_deg},
        "query_version": QUERY_VERSION,
        "selected_columns": columns,
        "row_count": len(table),
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
    }

    cache_dir.mkdir(parents=True, exist_ok=True)
    table.write(table_path, format="ascii.ecsv", overwrite=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return table, metadata


def to_dataframe(table: Any) -> Any:
    dataframe = table.to_pandas()
    for column in dataframe.columns:
        dataframe[column] = dataframe[column].replace([np.inf, -np.inf], np.nan)
    return dataframe


def prepare_catalogue(dataframe: Any, cluster: ClusterConfig) -> Any:
    """Add quality flags and derived astrometry without discarding the raw sample."""
    required_finite = list(REQUIRED_COLUMNS - {"source_id"})
    valid = dataframe[required_finite].notna().all(axis=1)
    valid &= dataframe["parallax_error"] > 0
    valid &= dataframe["pmra_error"] > 0
    valid &= dataframe["pmdec_error"] > 0
    observed_g_for_cut = dataframe["phot_g_mean_mag"] - 0.86 * cluster.av_prior
    valid &= observed_g_for_cut <= MAX_OBSERVED_G_MAG
    if cluster.max_g_mag is not None:
        valid &= dataframe["phot_g_mean_mag"] <= cluster.max_g_mag
    if "parallax_over_error" in dataframe:
        valid &= dataframe["parallax_over_error"] >= 5.0
    else:
        valid &= dataframe["parallax"] / dataframe["parallax_error"] >= 5.0
    if "ruwe" in dataframe:
        valid &= dataframe["ruwe"].fillna(np.inf) < 1.4
    if "visibility_periods_used" in dataframe:
        valid &= dataframe["visibility_periods_used"].fillna(0) >= 8
    if "phot_bp_rp_excess_factor" in dataframe:
        colour = dataframe["phot_bp_mean_mag"] - dataframe["phot_rp_mean_mag"]
        lower = 1.0 + 0.015 * colour**2
        upper = 1.3 + 0.06 * colour**2
        excess = dataframe["phot_bp_rp_excess_factor"]
        valid &= excess.isna() | ((excess > lower) & (excess < upper))

    dataframe = dataframe.copy()
    dataframe["quality_flag"] = valid
    dataframe["bp_rp"] = dataframe["phot_bp_mean_mag"] - dataframe["phot_rp_mean_mag"]
    dataframe["parallax_corrected"] = dataframe["parallax"] + PARALLAX_OFFSET_MAS
    return dataframe


def select_members(
    dataframe: Any,
    cluster: ClusterConfig,
    eps: float,
    min_samples: int,
    runtime: dict[str, Any],
) -> tuple[Any, dict[str, float]]:
    """Find the astrometric component closest to the configured cluster locus."""
    result = dataframe.copy()
    result["cluster_label"] = -1
    result["member_score"] = 0.0
    result["member_flag"] = False
    usable = result["quality_flag"].to_numpy()
    if usable.sum() < min_samples:
        raise RuntimeError(f"Only {usable.sum()} quality sources remain; cannot cluster {cluster.name}.")

    features = result.loc[usable, ["pmra", "pmdec", "parallax_corrected"]].to_numpy()
    scaler = runtime["RobustScaler"]()
    scaled = scaler.fit_transform(features)
    labels = runtime["DBSCAN"](eps=eps, min_samples=min_samples).fit_predict(scaled)
    result.loc[usable, "cluster_label"] = labels
    expected_scaled = scaler.transform(
        [[cluster.expected_pmra, cluster.expected_pmdec, cluster.expected_parallax]]
    )[0]
    candidates = [label for label in np.unique(labels) if label >= 0]
    if not candidates:
        raise RuntimeError(
            f"No DBSCAN component found for {cluster.name}. Try increasing --dbscan-eps."
        )

    def label_cost(label: int) -> float:
        component = scaled[labels == label]
        centre = np.median(component, axis=0)
        expected_distance = float(np.linalg.norm(centre - expected_scaled))
        sky = result.loc[usable].iloc[np.flatnonzero(labels == label)][["ra", "dec"]].to_numpy()
        sky_spread = float(np.nanmedian(np.hypot(sky[:, 0] - cluster.ra_deg, sky[:, 1] - cluster.dec_deg)))
        return expected_distance + 0.10 * sky_spread - 0.001 * len(component)

    chosen_label = min(candidates, key=label_cost)
    member_positions = np.flatnonzero(labels == chosen_label)
    component = scaled[member_positions]
    centre = np.median(component, axis=0)
    scatter = np.maximum(1.4826 * np.median(np.abs(component - centre), axis=0), 0.10)
    distance = np.sqrt(np.sum(((scaled - centre) / scatter) ** 2, axis=1))
    scores = np.exp(-0.5 * distance**2)
    usable_indices = result.index[usable]
    result.loc[usable_indices, "member_score"] = scores
    result.loc[usable_indices, "member_flag"] = (labels == chosen_label) & (scores >= MIN_MEMBER_SCORE)
    parameters = {
        "label": float(chosen_label),
        "pmra": float(np.median(result.loc[result.member_flag, "pmra"])),
        "pmdec": float(np.median(result.loc[result.member_flag, "pmdec"])),
        "parallax": float(np.median(result.loc[result.member_flag, "parallax_corrected"])),
        "n_members": float(result.member_flag.sum()),
    }
    return result, parameters


def estimate_distance(members: Any) -> tuple[float, float, float]:
    parallaxes = members["parallax_corrected"].to_numpy()
    errors = members["parallax_error"].to_numpy()
    weights = 1.0 / np.maximum(errors, 0.02) ** 2
    parallax = float(np.average(parallaxes, weights=weights))
    parallax_error = float(max(1.0 / np.sqrt(weights.sum()), 0.015))
    if parallax <= 0:
        raise RuntimeError("Member parallax is non-positive; no reliable distance modulus is available.")
    distance_pc = 1000.0 / parallax
    distance_modulus = 5.0 * np.log10(distance_pc) - 5.0
    return parallax, parallax_error, float(distance_modulus)


def add_hr_quantities(dataframe: Any, distance_modulus: float, av: float) -> Any:
    result = dataframe.copy()
    # Broad-band Gaia extinction coefficients used as a transparent first-pass model.
    result["abs_g"] = result["phot_g_mean_mag"] - distance_modulus - 0.86 * av
    result["bp_rp_0"] = result["bp_rp"] - 0.42 * av
    return result


def plot_diagnostics(dataframe: Any, cluster: ClusterConfig, output_dir: Path, runtime: dict[str, Any]) -> None:
    plt = runtime["plt"]
    output_dir.mkdir(parents=True, exist_ok=True)
    members = dataframe[dataframe.member_flag]
    figure, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    axes[0, 0].scatter(dataframe.ra, dataframe.dec, s=1, c="0.75", label="all")
    axes[0, 0].scatter(members.ra, members.dec, s=4, c="tab:red", label="members")
    axes[0, 0].set(xlabel="RA (deg)", ylabel="Dec (deg)", title=f"{cluster.name}: sky position")
    axes[0, 0].legend(markerscale=3)
    axes[0, 1].scatter(dataframe.pmra, dataframe.pmdec, s=1, c="0.65")
    axes[0, 1].scatter(members.pmra, members.pmdec, s=4, c="tab:red")
    axes[0, 1].set(xlabel="pmRA (mas/yr)", ylabel="pmDec (mas/yr)", title="Proper motion")
    axes[1, 0].hist(dataframe.parallax_corrected, bins=80, color="0.65")
    axes[1, 0].hist(members.parallax_corrected, bins=40, color="tab:red", alpha=0.8)
    axes[1, 0].set(xlabel="Corrected parallax (mas)", ylabel="Sources", title="Parallax")
    axes[1, 1].scatter(dataframe.pmra, dataframe.parallax_corrected, s=1, c="0.65")
    axes[1, 1].scatter(members.pmra, members.parallax_corrected, s=4, c="tab:red")
    axes[1, 1].set(xlabel="pmRA (mas/yr)", ylabel="Corrected parallax (mas)", title="Astrometric locus")
    figure.savefig(output_dir / f"{cluster.name.lower()}_diagnostics.png", dpi=180)
    plt.close(figure)


def plot_hr_diagram(
    dataframe: Any,
    cluster: ClusterConfig,
    output_dir: Path,
    fit: dict[str, Any] | None,
    runtime: dict[str, Any],
) -> None:
    plt = runtime["plt"]
    members = dataframe[dataframe.member_flag]
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for axis in axes:
        axis.scatter(dataframe.bp_rp_0, dataframe.abs_g, s=1, c="0.80", rasterized=True)
        axis.invert_yaxis()
        axis.set(xlabel=r"$(G_{BP} - G_{RP})_0$", ylabel=r"$M_G$", title=cluster.name)
    points = axes[0].scatter(members.bp_rp_0, members.abs_g, c=members.member_score, s=5, cmap="viridis")
    figure.colorbar(points, ax=axes[0], label="membership score")
    axes[1].scatter(members.bp_rp_0, members.abs_g, s=5, c="tab:red", rasterized=True)
    if fit is not None:
        label = (
            f"best isochrone\nlog age={fit['log_age']:.2f}, [Fe/H]={fit['feh']:.2f}, "
            f"A_V={fit['av']:.2f}, DM={fit['distance_modulus']:.2f}"
        )
        for axis in axes:
            axis.plot(fit["model_bp_rp_0"], fit["model_abs_g"], color="black", lw=1.5, label=label)
        axes[1].legend(loc="best", fontsize=8)
    figure.savefig(output_dir / f"{cluster.name.lower()}_hr_diagram.png", dpi=180)
    plt.close(figure)


def estimate_turnoff(members: Any) -> dict[str, float]:
    """Report a reproducible turn-off proxy; the isochrone fit is the age estimate."""
    candidates = members[(members.bp_rp_0 > -0.2) & (members.bp_rp_0 < 1.3)]
    if len(candidates) < 5:
        candidates = members
    brightest = candidates.nsmallest(max(1, min(5, len(candidates))), "abs_g")
    return {
        "turnoff_abs_g_proxy": float(brightest.abs_g.median()),
        "turnoff_bp_rp_proxy": float(brightest.bp_rp_0.median()),
    }


def ensure_mist_gaia_bolometric_tables() -> None:
    """Seed the MIST Gaia bolometric-correction tarball at the path expected by isochrones."""
    from isochrones.mist.bc import MISTBolometricCorrectionGrid

    grid = MISTBolometricCorrectionGrid(bands=["G", "BP", "RP"])
    tarball = Path(grid.get_tarball_file("UBVRIplus"))
    tarball.parent.mkdir(parents=True, exist_ok=True)
    if tarball.exists():
        try:
            with tarfile.open(tarball) as archive:
                if any(member.name.endswith(".UBVRIplus") for member in archive.getmembers()):
                    return
        except tarfile.TarError:
            pass
        tarball.unlink()
    urlretrieve("https://mist.science/BC_tables/v1/UBVRIplus.txz", tarball)
    with tarfile.open(tarball) as archive:
        if not any(member.name.endswith(".UBVRIplus") for member in archive.getmembers()):
            raise RuntimeError("Downloaded MIST UBVRIplus bolometric-correction archive has an unexpected layout.")


def build_builtin_isochrone_grid(cluster: ClusterConfig, cache_dir: Path, runtime: dict[str, Any]) -> Any:
    """Generate a cluster-focused Gaia isochrone grid with the installed MIST models."""
    pd = runtime["pd"]
    grid_path = cache_dir / f"{cluster.name.lower()}_mist_gaia_grid.csv"
    if grid_path.exists():
        return pd.read_csv(grid_path)
    try:
        from isochrones import get_ichrone
    except ImportError as error:
        raise RuntimeError(
            "Automatic isochrone fitting requires the 'isochrones' package or a CSV passed to --isochrone-grid."
        ) from error

    ensure_mist_gaia_bolometric_tables()
    ichrone = get_ichrone("mist", bands=["G", "BP", "RP"])
    frames = []
    age_grid = np.arange(cluster.log_age_min, cluster.log_age_max + 0.001, 0.1)
    feh_grid = np.arange(cluster.feh_prior - 0.2, cluster.feh_prior + 0.201, 0.1)
    for log_age in age_grid:
        for feh in feh_grid:
            model = ichrone.isochrone(log_age, feh=feh, distance=10.0, AV=0.0).dropna(
                subset=["G_mag", "BP_mag", "RP_mag"]
            )
            if len(model) < 4:
                continue
            frames.append(
                pd.DataFrame(
                    {
                        "log_age": float(log_age),
                        "feh": float(feh),
                        "abs_g": model["G_mag"].to_numpy(),
                        "bp_rp": (model["BP_mag"] - model["RP_mag"]).to_numpy(),
                    }
                )
            )
    if not frames:
        raise RuntimeError(f"No built-in MIST Gaia isochrone grid could be generated for {cluster.name}.")
    grid = pd.concat(frames, ignore_index=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    grid.to_csv(grid_path, index=False)
    return grid


def fit_isochrone_grid(
    members: Any,
    cluster: ClusterConfig,
    distance_modulus: float,
    cache_dir: Path,
    runtime: dict[str, Any],
    grid_path: Path | None = None,
) -> dict[str, Any]:
    """Grid-fit Gaia colour/magnitude data to a supplied or generated model grid."""
    pd = runtime["pd"]
    if grid_path is None:
        grid = build_builtin_isochrone_grid(cluster, cache_dir, runtime)
    else:
        grid = pd.read_csv(grid_path)
    required = {"log_age", "feh", "abs_g", "bp_rp"}
    if missing := required - set(grid.columns):
        raise ValueError(f"Isochrone grid is missing: {', '.join(sorted(missing))}")
    observed = members.sort_values("member_score", ascending=False).head(600)
    observations = observed[["bp_rp", "phot_g_mean_mag", "member_score"]].to_numpy()
    best: dict[str, Any] | None = None
    ages = np.sort(grid.log_age.unique())
    fehs = np.sort(grid.feh.unique())
    distance_grid = np.arange(distance_modulus - 0.35, distance_modulus + 0.351, 0.05)
    av_grid = np.arange(max(0.0, cluster.av_prior - 0.20), cluster.av_prior + 0.201, 0.025)
    for log_age in ages:
        for feh in fehs:
            model = grid[(grid.log_age == log_age) & (grid.feh == feh)].dropna(subset=["abs_g", "bp_rp"])
            if len(model) < 4:
                continue
            model_bp_rp_0 = model.bp_rp.to_numpy()
            model_abs_g = model.abs_g.to_numpy()
            for av in av_grid:
                model_colour = model_bp_rp_0 + 0.42 * av
                for dm in distance_grid:
                    model_mag = model_abs_g + dm + 0.86 * av
                    tree = runtime["cKDTree"](np.column_stack((model_colour / 0.04, model_mag / 0.12)))
                    distances, _ = tree.query(
                        np.column_stack((observations[:, 0] / 0.04, observations[:, 1] / 0.12))
                    )
                    weights = np.clip(observations[:, 2], 0.05, 1.0)
                    objective = float(np.average(distances**2, weights=weights))
                    candidate = {
                        "log_age": float(log_age),
                        "feh": float(feh),
                        "av": float(av),
                        "distance_modulus": float(dm),
                        "objective": objective,
                        "model_bp_rp_0": model_bp_rp_0.tolist(),
                        "model_abs_g": model_abs_g.tolist(),
                    }
                    if best is None or candidate["objective"] < best["objective"]:
                        best = candidate
    if best is None:
        raise RuntimeError("No usable isochrone model sequence was found in the grid.")
    return best


def write_summary(rows: list[dict[str, Any]], output_dir: Path, runtime: dict[str, Any]) -> None:
    dataframe = runtime["pd"].DataFrame(rows)
    dataframe.to_csv(output_dir / "cluster_summary.csv", index=False)
    lines = ["# Gaia Cluster Analysis Summary", "", dataframe.to_markdown(index=False), ""]
    (output_dir / "report.md").write_text("\n".join(lines))


def process_cluster(cluster: ClusterConfig, arguments: argparse.Namespace, runtime: dict[str, Any]) -> dict[str, Any]:
    table, metadata = fetch_catalogue(
        cluster,
        arguments.cache_dir,
        arguments.refresh,
        arguments.offline,
        runtime,
    )
    dataframe = prepare_catalogue(to_dataframe(table), cluster)
    dataframe, membership = select_members(
        dataframe, cluster, arguments.dbscan_eps, arguments.min_samples, runtime
    )
    members = dataframe[dataframe.member_flag].copy()
    if len(members) < 5:
        raise RuntimeError(f"Only {len(members)} members selected for {cluster.name}.")
    parallax, parallax_error, distance_modulus = estimate_distance(members)
    dataframe = add_hr_quantities(dataframe, distance_modulus, cluster.av_prior)
    members = dataframe[dataframe.member_flag].copy()
    fit = None
    try:
        fit = fit_isochrone_grid(
            members,
            cluster,
            distance_modulus,
            arguments.cache_dir,
            runtime,
            grid_path=arguments.isochrone_grid,
        )
    except RuntimeError as error:
        missing_builtin = "requires the 'isochrones' package" in str(error)
        if arguments.isochrone_grid is not None or not missing_builtin:
            raise
        LOG.warning(
            "Skipping built-in isochrone fit for %s because the 'isochrones' package is unavailable.",
            cluster.name,
        )
    if fit is not None:
        dataframe = add_hr_quantities(dataframe, fit["distance_modulus"], fit["av"])
        members = dataframe[dataframe.member_flag].copy()
    plot_diagnostics(dataframe, cluster, arguments.output_dir, runtime)
    plot_hr_diagram(dataframe, cluster, arguments.output_dir, fit, runtime)
    turnoff = estimate_turnoff(members)
    member_path = arguments.output_dir / f"{cluster.name.lower()}_members.csv"
    members.to_csv(member_path, index=False)
    (arguments.output_dir / f"{cluster.name.lower()}_metadata.json").write_text(
        json.dumps({"query": metadata, "membership": membership, "fit": fit, "turnoff": turnoff}, indent=2)
        + "\n"
    )
    return {
        "cluster": cluster.name,
        "n_queried": len(dataframe),
        "n_quality": int(dataframe.quality_flag.sum()),
        "n_members": int(len(members)),
        "parallax_mas": parallax,
        "parallax_error_mas": parallax_error,
        "distance_modulus": fit["distance_modulus"] if fit else distance_modulus,
        "av": fit["av"] if fit else cluster.av_prior,
        "feh": fit["feh"] if fit else cluster.feh_prior,
        "log_age": fit["log_age"] if fit else np.nan,
        "age_myr": 10 ** (fit["log_age"] - 6) if fit else np.nan,
        "fit_objective": fit["objective"] if fit else np.nan,
        **turnoff,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clusters", nargs="+", choices=sorted(CLUSTERS), default=sorted(CLUSTERS))
    parser.add_argument("--cache-dir", type=Path, default=Path("data"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--isochrone-grid",
        type=Path,
        help="Optional CSV with log_age, feh, abs_g, bp_rp columns; otherwise use built-in MIST Gaia isochrones",
    )
    parser.add_argument("--refresh", action="store_true", help="Ignore cached Gaia results")
    parser.add_argument("--offline", action="store_true", help="Require cached Gaia results")
    parser.add_argument("--dbscan-eps", type=float, default=0.55)
    parser.add_argument("--min-samples", type=int, default=12)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    logging.basicConfig(level=logging.DEBUG if arguments.verbose else logging.INFO, format="%(levelname)s: %(message)s")
    if arguments.isochrone_grid is not None and not arguments.isochrone_grid.is_file():
        raise FileNotFoundError(f"Isochrone grid does not exist: {arguments.isochrone_grid}")
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    runtime = require_runtime()
    summaries: list[dict[str, Any]] = []
    for name in arguments.clusters:
        summaries.append(process_cluster(CLUSTERS[name], arguments, runtime))
    if all(np.isfinite(row["log_age"]) for row in summaries):
        summaries.sort(key=lambda row: row["log_age"])
        for rank, row in enumerate(summaries, start=1):
            row["age_rank_youngest_first"] = rank
    write_summary(summaries, arguments.output_dir, runtime)
    LOG.info("Completed %d cluster(s). Results: %s", len(summaries), arguments.output_dir)


if __name__ == "__main__":
    main()
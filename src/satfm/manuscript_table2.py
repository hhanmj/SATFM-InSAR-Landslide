"""Manuscript Table 2 synthetic benchmark reproduction.

This module is fully self-contained. It preserves the 15 m grid, 150 x 150
scene, 30 Monte Carlo realizations, published LOS geometry, and legacy-linear
Eq. (12) edge kernel used for the manuscript Table 2 synthetic experiment.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, uniform_filter
from scipy.sparse.linalg import LinearOperator, cg


GRID_NY = 150
GRID_NX = 150
PIXEL_SIZE_M = 15.0
MC_SEEDS = 30
REPRESENTATIVE_SEED = 20260424
CG_RTOL = 2.0e-4
CG_MAXITER = 180

LOS_VECTORS = np.array(
    [
        [0.4297962784767151, 0.09018063545227051, 0.8984111547470093],
        [-0.6052070260047913, 0.12709592282772064, 0.7858569025993347],
        [-0.6607598066329956, -0.14975716173648834, 0.7355061173439026],
    ],
    dtype=np.float64,
)
SIGMA2_TRACK = np.array(
    [2.1379730992704067, 2.56630558102138, 3.159767145923468],
    dtype=np.float64,
)

LAMBDA_D_BASE = np.array([1.0e-2, 2.4e-2, 5.2e-2], dtype=np.float64)
LAMBDA_T_BASE = np.array([5.0e-3, 1.5e-2, 3.0e-2], dtype=np.float64)
ADAPTIVE_GAMMA = 1.05
LAMBDA_S = 2.4e-2

METHODS = ["WLS", "APFM", "SPFM", "Full SATFM"]


def build_grid(ny: int, nx: int) -> tuple[np.ndarray, np.ndarray]:
    y = np.linspace(-1.0, 1.0, ny)
    x = np.linspace(-1.0, 1.0, nx)
    return np.meshgrid(x, y)


def gaussian2(x: np.ndarray, y: np.ndarray, cx: float, cy: float, sx: float, sy: float) -> np.ndarray:
    return np.exp(-(((x - cx) / sx) ** 2 + ((y - cy) / sy) ** 2))


def sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def rotate_horizontal(east: np.ndarray, north: np.ndarray, angle_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    angle = np.deg2rad(angle_deg)
    ca = np.cos(angle)
    sa = np.sin(angle)
    return ca * east - sa * north, sa * east + ca * north


def downhill_from_dem(dem: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dz_drow, dz_dcol = np.gradient(dem, PIXEL_SIZE_M, PIXEL_SIZE_M)
    dz_dx = dz_dcol
    dz_dy = -dz_drow
    slope = np.sqrt(dz_dx * dz_dx + dz_dy * dz_dy)
    down_e = -dz_dx / np.clip(slope, 1.0e-9, None)
    down_n = -dz_dy / np.clip(slope, 1.0e-9, None)
    slope_deg = np.degrees(np.arctan(slope))
    return down_e, down_n, dz_dx, dz_dy, slope_deg


def terrain_parallel_vertical(east: np.ndarray, north: np.ndarray, dz_dx: np.ndarray, dz_dy: np.ndarray) -> np.ndarray:
    return dz_dx * east + dz_dy * north


def make_blind_scene(ny: int = GRID_NY, nx: int = GRID_NX) -> dict:
    """Synthetic scene used for manuscript Table 2."""

    xx, yy = build_grid(ny, nx)
    g_head = gaussian2(xx, yy, -0.30, -0.55, 0.48, 0.34)
    g_core = gaussian2(xx, yy, 0.12, -0.03, 0.42, 0.42)
    g_toe = gaussian2(xx, yy, 0.36, 0.46, 0.32, 0.34)
    g_left = gaussian2(xx, yy, -0.52, 0.08, 0.30, 0.56)
    g_right = gaussian2(xx, yy, 0.55, 0.10, 0.34, 0.52)
    shear = np.exp(-((xx + 0.18 * np.sin(1.4 * np.pi * yy)) ** 2) / 0.026)

    dem = (
        700.0
        - 285.0 * yy
        + 45.0 * xx
        + 34.0 * g_head
        - 42.0 * g_core
        - 28.0 * g_toe
        + 22.0 * np.sin(1.1 * np.pi * xx) * np.cos(0.65 * np.pi * yy)
        + 12.0 * np.sin(2.3 * np.pi * (xx + 0.18 * yy))
    )
    down_e, down_n, dz_dx, dz_dy, _ = downhill_from_dem(dem)

    speed = 12.0 + 9.0 * sigmoid(3.0 * (yy + 0.58)) + 18.0 * g_core + 10.0 * g_toe + 7.0 * shear - 4.0 * g_left
    speed = np.clip(speed, 6.0, 48.0)
    deflection = (
        30.0 * g_left
        - 26.0 * g_right
        + 16.0 * shear * np.sin(1.15 * np.pi * yy + 0.35 * np.pi * xx)
        + 10.0 * np.sin(np.pi * yy) * np.exp(-1.7 * xx * xx)
    )
    east, north = rotate_horizontal(down_e * speed, down_n * speed, deflection)
    vertical_spf = terrain_parallel_vertical(east, north, dz_dx, dz_dy)
    normal_deform = (
        -5.2 * g_core
        - 3.8 * g_toe
        + 3.2 * g_head
        + 2.4 * shear * np.sin(np.pi * yy)
        - 1.8 * np.sin(1.7 * np.pi * xx) * np.exp(-1.8 * yy * yy)
    )
    vertical = 0.45 * vertical_spf + normal_deform
    return {"truth": {"E": east, "N": north, "U": vertical}, "vertical_spf": vertical_spf, "dem": dem}


def nan_gaussian(arr: np.ndarray, mask: np.ndarray, sigma: float) -> np.ndarray:
    vals = np.where(mask & np.isfinite(arr), arr, 0.0)
    w = gaussian_filter((mask & np.isfinite(arr)).astype(np.float64), sigma=sigma, mode="nearest")
    sm = gaussian_filter(vals, sigma=sigma, mode="nearest")
    return np.where(w > 1.0e-8, sm / w, np.nan)


def robust_scale(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    med = np.nanmedian(x, axis=0)
    q25 = np.nanpercentile(x, 25, axis=0)
    q75 = np.nanpercentile(x, 75, axis=0)
    scale = q75 - q25
    fallback = np.nanstd(x, axis=0)
    scale = np.where(scale > 1.0e-9, scale, np.where(fallback > 1.0e-9, fallback, 1.0))
    return (x - med) / scale, med, scale


def normalize2(e: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    norm = np.sqrt(e * e + n * n)
    ok = norm > 1.0e-12
    ee = np.zeros_like(e, dtype=np.float64)
    nn = np.zeros_like(n, dtype=np.float64)
    ee[ok] = e[ok] / norm[ok]
    nn[ok] = n[ok] / norm[ok]
    return ee, nn, ok


def local_std(arr: np.ndarray, mask: np.ndarray, size: int = 5) -> np.ndarray:
    vals = np.where(mask & np.isfinite(arr), arr, 0.0).astype(np.float64)
    w = uniform_filter((mask & np.isfinite(arr)).astype(np.float64), size=size, mode="nearest")
    mean = uniform_filter(vals, size=size, mode="nearest") / np.clip(w, 1.0e-9, None)
    mean2 = uniform_filter(vals * vals, size=size, mode="nearest") / np.clip(w, 1.0e-9, None)
    out = np.sqrt(np.maximum(mean2 - mean * mean, 0.0))
    out[w <= 1.0e-6] = np.nan
    return out


def mask_aware_gradient_2d(arr: np.ndarray, dx: float, dy: float, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    filled = np.where(mask & np.isfinite(arr), arr, np.nan)
    if not np.all(np.isfinite(filled)):
        replacement = nan_gaussian(filled, mask, sigma=1.0)
        filled = np.where(np.isfinite(filled), filled, replacement)
        filled = np.where(np.isfinite(filled), filled, np.nanmedian(arr[mask]))
    return np.gradient(filled, dy, dx)


def solve_weighted_local(los_rates: np.ndarray, los_vectors: np.ndarray, sigma2: np.ndarray) -> dict:
    los_rates = np.asarray(los_rates, dtype=np.float64)
    los_vectors = np.asarray(los_vectors, dtype=np.float64)
    sigma2 = np.asarray(sigma2, dtype=np.float64)
    _, n_pix = los_rates.shape
    inv_sigma2 = 1.0 / np.clip(sigma2, 1.0e-8, None)
    a_outer = np.einsum("si,sj->sij", los_vectors, los_vectors)
    fisher_obs = np.einsum("sn,sij->nij", inv_sigma2, a_outer)
    rhs = np.einsum("sn,si->ni", inv_sigma2 * los_rates, los_vectors)
    x = np.full((n_pix, 3), np.nan, dtype=np.float64)
    det = np.linalg.det(fisher_obs)
    ok = np.isfinite(det) & (np.abs(det) > 1.0e-14)
    if np.any(ok):
        q = np.linalg.inv(fisher_obs[ok])
        x[ok] = np.einsum("nij,nj->ni", q, rhs[ok])
    return {"x": x, "fisher_obs": fisher_obs, "rhs": rhs, "valid": ok}


def solve_strict_basis(
    los_rates: np.ndarray,
    los_vectors: np.ndarray,
    sigma2: np.ndarray,
    basis: np.ndarray,
    nonnegative_first: bool = False,
) -> np.ndarray:
    g = np.einsum("si,nim->snm", los_vectors, basis)
    inv_sigma2 = 1.0 / np.clip(sigma2, 1.0e-8, None)
    fisher = np.einsum("sn,snm,snk->nmk", inv_sigma2, g, g)
    rhs = np.einsum("sn,snm->nm", inv_sigma2 * los_rates, g)
    coeff = np.full_like(rhs, np.nan, dtype=np.float64)
    det = np.linalg.det(fisher)
    ok = np.isfinite(det) & (np.abs(det) > 1.0e-14)
    if np.any(ok):
        coeff[ok] = np.linalg.solve(fisher[ok], rhs[ok][..., None])[..., 0]
    if nonnegative_first:
        bad = np.isfinite(coeff[:, 0]) & (coeff[:, 0] < 0.0)
        if np.any(bad):
            denom = np.clip(fisher[bad, 1, 1], 1.0e-12, None)
            coeff[bad, 0] = 0.0
            coeff[bad, 1] = rhs[bad, 1] / denom
    return np.einsum("nim,nm->ni", basis, coeff)


def solve_strict_apfm(los_rates: np.ndarray, los_vectors: np.ndarray, sigma2: np.ndarray, geom: dict) -> np.ndarray:
    n_pix = los_rates.shape[1]
    basis = np.zeros((n_pix, 3, 2), dtype=np.float64)
    basis[:, 0, 0] = geom["down_e"].reshape(-1)
    basis[:, 1, 0] = geom["down_n"].reshape(-1)
    basis[:, 2, 1] = 1.0
    return solve_strict_basis(los_rates, los_vectors, sigma2, basis, nonnegative_first=True)


def solve_strict_spfm(los_rates: np.ndarray, los_vectors: np.ndarray, sigma2: np.ndarray, geom: dict) -> np.ndarray:
    n_pix = los_rates.shape[1]
    dz_dx = geom["dz_dx"].reshape(-1)
    dz_dy = geom["dz_dy"].reshape(-1)
    basis = np.zeros((n_pix, 3, 2), dtype=np.float64)
    basis[:, 0, 0] = 1.0
    basis[:, 2, 0] = dz_dx
    basis[:, 1, 1] = 1.0
    basis[:, 2, 1] = dz_dy
    return solve_strict_basis(los_rates, los_vectors, sigma2, basis)


def kmeans_init(x: np.ndarray, k: int, rng: np.random.Generator, n_iter: int = 20) -> np.ndarray:
    n = x.shape[0]
    centers = np.empty((k, x.shape[1]), dtype=np.float64)
    centers[0] = x[rng.integers(0, n)]
    dist2 = np.sum((x - centers[0]) ** 2, axis=1)
    for i in range(1, k):
        probs = dist2 / np.clip(np.sum(dist2), 1.0e-12, None)
        centers[i] = x[rng.choice(n, p=probs)]
        dist2 = np.minimum(dist2, np.sum((x - centers[i]) ** 2, axis=1))
    for _ in range(n_iter):
        d2 = np.sum((x[:, None, :] - centers[None, :, :]) ** 2, axis=2)
        labels = np.argmin(d2, axis=1)
        for i in range(k):
            if np.any(labels == i):
                centers[i] = np.mean(x[labels == i], axis=0)
    return centers


def fit_diag_gmm(x: np.ndarray, k: int = 3, n_iter: int = 60, seed: int = 20260423) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    centers = kmeans_init(x, k, rng)
    cov = np.tile(np.var(x, axis=0) + 1.0e-3, (k, 1))
    weights = np.full(k, 1.0 / k)
    eps = 1.0e-6
    for _ in range(n_iter):
        logp = np.empty((x.shape[0], k), dtype=np.float64)
        for j in range(k):
            v = np.clip(cov[j], eps, None)
            log_det = np.sum(np.log(v))
            mahal = np.sum((x - centers[j]) ** 2 / v, axis=1)
            logp[:, j] = np.log(weights[j] + eps) - 0.5 * (log_det + mahal)
        logp -= np.max(logp, axis=1, keepdims=True)
        resp = np.exp(logp)
        resp /= np.clip(np.sum(resp, axis=1, keepdims=True), eps, None)
        nk = np.sum(resp, axis=0) + eps
        weights = nk / np.sum(nk)
        centers = (resp.T @ x) / nk[:, None]
        for j in range(k):
            cov[j] = (resp[:, j, None] * (x - centers[j]) ** 2).sum(axis=0) / nk[j]
        cov = np.clip(cov, 1.0e-4, None)
    return weights, centers, cov


def predict_diag_gmm_chunks(
    x: np.ndarray,
    weights: np.ndarray,
    centers: np.ndarray,
    cov: np.ndarray,
    chunk_size: int = 250000,
) -> np.ndarray:
    k = len(weights)
    out = np.empty((x.shape[0], k), dtype=np.float32)
    eps = 1.0e-12
    for start in range(0, x.shape[0], chunk_size):
        xx = x[start : start + chunk_size]
        logp = np.empty((xx.shape[0], k), dtype=np.float64)
        for j in range(k):
            v = np.clip(cov[j], 1.0e-6, None)
            logp[:, j] = np.log(weights[j] + eps) - 0.5 * (
                np.sum(np.log(v)) + np.sum((xx - centers[j]) ** 2 / v, axis=1)
            )
        logp -= np.max(logp, axis=1, keepdims=True)
        resp = np.exp(logp)
        resp /= np.clip(np.sum(resp, axis=1, keepdims=True), eps, None)
        out[start : start + chunk_size] = resp.astype(np.float32)
    return out


def build_edges(
    mask: np.ndarray,
    memberships: np.ndarray,
    slope_deg: np.ndarray,
    t0_e: np.ndarray,
    t0_n: np.ndarray,
    c_dir: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    idx = -np.ones(mask.shape, dtype=np.int32)
    rr, cc = np.where(mask)
    idx[rr, cc] = np.arange(rr.size, dtype=np.int32)
    edges_i = []
    edges_j = []
    weights = []
    for dr, dc in [(0, 1), (1, 0)]:
        a = mask.copy()
        b = np.roll(mask, shift=(-dr, -dc), axis=(0, 1))
        if dr == 0:
            a[:, -1] = False
            b[:, -1] = False
        if dc == 0:
            a[-1, :] = False
            b[-1, :] = False
        pair = a & b
        r0, c0 = np.where(pair)
        r1, c1 = r0 + dr, c0 + dc
        ii = idx[r0, c0]
        jj = idx[r1, c1]
        dm2 = np.sum((memberships[ii] - memberships[jj]) ** 2, axis=1)
        ds = np.abs(slope_deg[r0, c0] - slope_deg[r1, c1])
        dot = np.clip(t0_e[r0, c0] * t0_e[r1, c1] + t0_n[r0, c0] * t0_n[r1, c1], -1.0, 1.0)
        da = np.arccos(dot)
        if c_dir is not None:
            c_pair = np.clip(c_dir[ii] * c_dir[jj], 0.0, 1.0)
            da = c_pair * da
        tau_m = 0.55
        tau_s = 8.0
        tau_a = np.deg2rad(35.0)
        w = np.exp(-dm2 / (tau_m * tau_m)) * np.exp(-ds / (tau_s * tau_s)) * np.exp(-da / (tau_a * tau_a))
        ok = np.isfinite(w) & (w > 1.0e-6)
        edges_i.append(ii[ok].astype(np.int32))
        edges_j.append(jj[ok].astype(np.int32))
        weights.append(w[ok].astype(np.float64))
    return np.concatenate(edges_i), np.concatenate(edges_j), np.concatenate(weights)


def solve_adaptive_global(
    fisher_obs: np.ndarray,
    rhs: np.ndarray,
    x0: np.ndarray,
    ndir: np.ndarray,
    nsurf: np.ndarray,
    lambda_d: np.ndarray,
    lambda_t: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    edge_w: np.ndarray,
    lambda_s: float,
    rtol: float = CG_RTOL,
    maxiter: int = CG_MAXITER,
) -> tuple[np.ndarray, int, int]:
    n = rhs.shape[0]
    diag = np.einsum("nii->ni", fisher_obs).copy()
    diag[:, 0] += lambda_d * ndir[:, 0] ** 2 + lambda_t * nsurf[:, 0] ** 2
    diag[:, 1] += lambda_d * ndir[:, 1] ** 2 + lambda_t * nsurf[:, 1] ** 2
    diag[:, 2] += lambda_t * nsurf[:, 2] ** 2
    deg = np.zeros(n, dtype=np.float64)
    np.add.at(deg, edge_i, edge_w)
    np.add.at(deg, edge_j, edge_w)
    diag += lambda_s * deg[:, None]
    inv_diag = 1.0 / np.clip(diag, 1.0e-10, None)

    def matvec(vec: np.ndarray) -> np.ndarray:
        x = vec.reshape(n, 3)
        y = np.einsum("nij,nj->ni", fisher_obs, x)
        hd = ndir[:, 0] * x[:, 0] + ndir[:, 1] * x[:, 1]
        y[:, 0] += lambda_d * hd * ndir[:, 0]
        y[:, 1] += lambda_d * hd * ndir[:, 1]
        ht = np.einsum("ni,ni->n", nsurf, x)
        y += lambda_t[:, None] * ht[:, None] * nsurf
        if edge_w.size:
            diff = x[edge_i] - x[edge_j]
            contrib = lambda_s * edge_w[:, None] * diff
            np.add.at(y, edge_i, contrib)
            np.add.at(y, edge_j, -contrib)
        return y.reshape(-1)

    def psolve(vec: np.ndarray) -> np.ndarray:
        return (vec.reshape(n, 3) * inv_diag).reshape(-1)

    op = LinearOperator((3 * n, 3 * n), matvec=matvec, dtype=np.float64)
    pre = LinearOperator((3 * n, 3 * n), matvec=psolve, dtype=np.float64)
    info_holder = {"iterations": 0}

    def callback(_: np.ndarray) -> None:
        info_holder["iterations"] += 1

    sol, info = cg(
        op,
        rhs.reshape(-1),
        x0=x0.reshape(-1),
        M=pre,
        rtol=rtol,
        atol=0.0,
        maxiter=maxiter,
        callback=callback,
    )
    return sol.reshape(n, 3), int(info), int(info_holder["iterations"])


def build_scene_geometry(scene: dict) -> dict:
    mask = np.ones(scene["dem"].shape, dtype=bool)
    dem_smooth = nan_gaussian(scene["dem"], mask, sigma=1.2)
    dz_drow, dz_dcol = mask_aware_gradient_2d(dem_smooth, PIXEL_SIZE_M, PIXEL_SIZE_M, mask)
    dz_dx = dz_dcol
    dz_dy = -dz_drow
    slope_deg = np.degrees(np.arctan(np.sqrt(dz_dx * dz_dx + dz_dy * dz_dy)))
    slope_norm = np.sqrt(dz_dx * dz_dx + dz_dy * dz_dy)
    down_e = -dz_dx / np.clip(slope_norm, 1.0e-9, None)
    down_n = -dz_dy / np.clip(slope_norm, 1.0e-9, None)
    aspect = np.mod(np.arctan2(-dz_dx, -dz_dy), 2.0 * np.pi)
    surf = np.stack([-dz_dx, -dz_dy, np.ones_like(dz_dx)], axis=-1)
    surf /= np.clip(np.linalg.norm(surf, axis=-1, keepdims=True), 1.0e-12, None)
    dem_large = nan_gaussian(dem_smooth, mask, sigma=3.0)
    rough = local_std(dem_smooth - dem_large, mask, size=7)
    curvature = np.gradient(dz_dx, PIXEL_SIZE_M, axis=1) + np.gradient(dz_dy, PIXEL_SIZE_M, axis=0)
    return {
        "mask": mask,
        "surf": surf.reshape(-1, 3),
        "down_e": down_e,
        "down_n": down_n,
        "dz_dx": dz_dx,
        "dz_dy": dz_dy,
        "slope_deg": slope_deg,
        "aspect": aspect,
        "rough": rough,
        "curvature": curvature,
        "dem_smooth": dem_smooth,
        "dem_large": dem_large,
    }


def simulate_observation(scene: dict, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    truth_stack = np.stack([scene["truth"]["E"], scene["truth"]["N"], scene["truth"]["U"]], axis=0)
    obs_clean = np.einsum("si,ihw->shw", LOS_VECTORS, truth_stack)
    sigma = np.sqrt(SIGMA2_TRACK)[:, None, None]
    noise = rng.normal(0.0, 1.0, size=obs_clean.shape) * sigma
    return obs_clean + noise


def build_direction_prior(x_wls: np.ndarray, geom: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ny, nx = geom["mask"].shape
    raw_e = x_wls[:, 0].reshape(ny, nx)
    raw_n = x_wls[:, 1].reshape(ny, nx)
    hspeed = np.sqrt(raw_e * raw_e + raw_n * raw_n)
    speed_lo, speed_hi = np.percentile(hspeed, [35, 80])
    speed_conf = np.clip((hspeed - speed_lo) / max(speed_hi - speed_lo, 1.0e-6), 0.0, 1.0)
    sm_e = nan_gaussian(raw_e * speed_conf, geom["mask"], sigma=2.0)
    sm_n = nan_gaussian(raw_n * speed_conf, geom["mask"], sigma=2.0)
    sm_w = nan_gaussian(speed_conf, geom["mask"], sigma=2.0)
    sm_e = np.where(sm_w > 1.0e-6, sm_e / np.clip(sm_w, 1.0e-6, None), 0.0)
    sm_n = np.where(sm_w > 1.0e-6, sm_n / np.clip(sm_w, 1.0e-6, None), 0.0)
    sm_e, sm_n, _ = normalize2(sm_e, sm_n)
    slope_conf = np.clip((geom["slope_deg"] - 5.0) / 20.0, 0.0, 1.0)
    down_e = np.sin(geom["aspect"])
    down_n = np.cos(geom["aspect"])
    smooth_resultant = np.sqrt(sm_e * sm_e + sm_n * sm_n)
    raw_e_unit, raw_n_unit, _ = normalize2(raw_e, raw_n)
    rel_raw = speed_conf * np.clip(smooth_resultant, 0.0, 1.0)
    blend_e = rel_raw * raw_e_unit + (1.0 - rel_raw) * (0.72 * sm_e + 0.28 * down_e * slope_conf)
    blend_n = rel_raw * raw_n_unit + (1.0 - rel_raw) * (0.72 * sm_n + 0.28 * down_n * slope_conf)
    t0_e, t0_n, _ = normalize2(blend_e, blend_n)
    t0_rel = np.clip(0.65 * rel_raw + 0.35 * slope_conf, 0.0, 1.0)
    return t0_e, t0_n, t0_rel, speed_conf


def build_data_driven_membership(x_wls: np.ndarray, geom: dict, seed: int, include_xy: bool = True) -> np.ndarray:
    mask = geom["mask"]
    ny, nx = mask.shape
    raw_e = x_wls[:, 0].reshape(ny, nx)
    raw_n = x_wls[:, 1].reshape(ny, nx)
    raw_u = x_wls[:, 2].reshape(ny, nx)
    raw_h = np.sqrt(raw_e * raw_e + raw_n * raw_n)
    speed_lo, speed_hi = np.nanpercentile(raw_h, [35, 85])
    speed_conf = np.clip((raw_h - speed_lo) / max(speed_hi - speed_lo, 1.0e-6), 0.0, 1.0)
    sm_e = nan_gaussian(raw_e * speed_conf, mask, sigma=2.4)
    sm_n = nan_gaussian(raw_n * speed_conf, mask, sigma=2.4)
    sm_u = nan_gaussian(raw_u * speed_conf, mask, sigma=2.4)
    sm_w = nan_gaussian(speed_conf, mask, sigma=2.4)
    sm_e = np.where(sm_w > 1.0e-6, sm_e / np.clip(sm_w, 1.0e-6, None), 0.0)
    sm_n = np.where(sm_w > 1.0e-6, sm_n / np.clip(sm_w, 1.0e-6, None), 0.0)
    sm_u = np.where(sm_w > 1.0e-6, sm_u / np.clip(sm_w, 1.0e-6, None), 0.0)
    sm_h = np.sqrt(sm_e * sm_e + sm_n * sm_n)
    sm_dir = np.mod(np.arctan2(sm_e, sm_n), 2.0 * np.pi)
    aspect = geom["aspect"]
    dir_misfit = np.abs(np.arctan2(np.sin(sm_dir - aspect), np.cos(sm_dir - aspect)))
    ns = geom["surf"].reshape(ny, nx, 3)
    normal_proxy = ns[..., 0] * sm_e + ns[..., 1] * sm_n + ns[..., 2] * sm_u
    ycoord = np.linspace(-1.0, 1.0, ny)[:, None] * np.ones((ny, nx))
    xcoord = np.ones((ny, nx)) * np.linspace(-1.0, 1.0, nx)[None, :]
    dem_norm = (geom["dem_smooth"] - np.nanmedian(geom["dem_smooth"])) / max(
        np.nanpercentile(geom["dem_smooth"], 90) - np.nanpercentile(geom["dem_smooth"], 10),
        1.0e-6,
    )
    feature_items = []
    if include_xy:
        feature_items.extend([("x", xcoord), ("y", ycoord)])
    feature_items.extend(
        [
            ("elevation", dem_norm),
            ("slope", geom["slope_deg"]),
            ("sin_aspect", np.sin(aspect)),
            ("cos_aspect", np.cos(aspect)),
            ("smooth_horizontal_speed", sm_h),
            ("smooth_vertical_speed", sm_u),
            ("direction_misfit_to_aspect", dir_misfit),
            ("terrain_normal_velocity_proxy", normal_proxy),
            ("curvature", geom["curvature"]),
            ("roughness", geom["rough"]),
        ]
    )
    feature_names = [name for name, _ in feature_items]
    features = np.column_stack([arr.reshape(-1) for _, arr in feature_items])
    good = np.isfinite(features).all(axis=1)
    med = np.nanmedian(np.where(np.isfinite(features), features, np.nan), axis=0)
    features = np.where(np.isfinite(features), features, med)
    scaled, center, scale = robust_scale(features)
    fit_idx = np.where(good & mask.reshape(-1))[0]
    weights, centers, cov = fit_diag_gmm(scaled[fit_idx], k=3, n_iter=80, seed=seed)
    membership = predict_diag_gmm_chunks(scaled, weights, centers, cov).astype(np.float64)
    centers_original = centers * scale + center
    score_col = (
        centers_original[:, feature_names.index("smooth_horizontal_speed")]
        + 5.0 * np.abs(centers_original[:, feature_names.index("terrain_normal_velocity_proxy")])
        + 4.0 * centers_original[:, feature_names.index("direction_misfit_to_aspect")]
    )
    order = np.argsort(score_col)
    return membership[:, order]


def solve_one_seed(scene: dict, geom: dict, seed: int) -> tuple[dict[str, np.ndarray], dict]:
    obs_stack = simulate_observation(scene, seed)
    n_pix = obs_stack.shape[1] * obs_stack.shape[2]
    los_rates = obs_stack.reshape(3, n_pix)
    sigma2 = np.repeat(SIGMA2_TRACK[:, None], n_pix, axis=1)
    wls = solve_weighted_local(los_rates, LOS_VECTORS, sigma2)
    apfm = solve_strict_apfm(los_rates, LOS_VECTORS, sigma2, geom)
    spfm = solve_strict_spfm(los_rates, LOS_VECTORS, sigma2, geom)
    membership = build_data_driven_membership(wls["x"], geom, seed, include_xy=True)
    t0_e, t0_n, t0_rel, speed_conf = build_direction_prior(wls["x"], geom)
    ndir = np.column_stack([-t0_n.reshape(-1), t0_e.reshape(-1)])
    dem_resid = np.abs(geom["dem_smooth"] - geom["dem_large"])
    dem_q = float(np.nanpercentile(dem_resid, 80))
    terrain_quality = np.clip(1.0 - dem_resid / max(dem_q, 1.0e-6), 0.0, 1.0)
    slope_conf = np.clip((geom["slope_deg"] - 5.0) / 20.0, 0.0, 1.0)
    c_dir = np.clip(t0_rel.reshape(-1) * speed_conf.reshape(-1), 0.0, 1.0)
    c_tan = np.clip(0.15 + 0.85 * slope_conf * terrain_quality, 0.05, 1.0).reshape(-1)
    lambda_d = ADAPTIVE_GAMMA * c_dir * (membership @ LAMBDA_D_BASE)
    lambda_t = ADAPTIVE_GAMMA * c_tan * (membership @ LAMBDA_T_BASE)
    edge_i, edge_j, edge_w = build_edges(geom["mask"], membership, geom["slope_deg"], t0_e, t0_n, c_dir=c_dir)
    satfm, cg_info, cg_iter = solve_adaptive_global(
        wls["fisher_obs"],
        wls["rhs"],
        spfm,
        ndir,
        geom["surf"],
        lambda_d,
        lambda_t,
        edge_i,
        edge_j,
        edge_w,
        LAMBDA_S,
    )
    estimates = {"WLS": wls["x"], "APFM": apfm, "SPFM": spfm, "Full SATFM": satfm}
    diagnostics = {
        "seed": int(seed),
        "cg_info": int(cg_info),
        "cg_iterations": int(cg_iter),
        "edge_count": int(edge_w.size),
        "edge_kernel": "legacy_linear",
        "edge_direction_confidence_gate": True,
    }
    return estimates, diagnostics


def compute_enu_rmse(truth: np.ndarray, estimate: np.ndarray) -> dict[str, float]:
    err = estimate - truth
    return {
        "East RMSE (mm/yr)": float(np.sqrt(np.nanmean(err[:, 0] ** 2))),
        "North RMSE (mm/yr)": float(np.sqrt(np.nanmean(err[:, 1] ** 2))),
        "Up RMSE (mm/yr)": float(np.sqrt(np.nanmean(err[:, 2] ** 2))),
    }


def _format_pm(values: np.ndarray) -> str:
    return f"{np.nanmean(values):.3f} +/- {np.nanstd(values, ddof=1):.3f}"


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method in METHODS:
        group = detail[detail["Method"] == method]
        row = {"Method": method}
        for col in ["East RMSE (mm/yr)", "North RMSE (mm/yr)", "Up RMSE (mm/yr)"]:
            vals = group[col].to_numpy(dtype=np.float64)
            row[col] = _format_pm(vals)
            row[f"{col}_mean"] = float(np.nanmean(vals))
            row[f"{col}_sd"] = float(np.nanstd(vals, ddof=1))
        rows.append(row)
    return pd.DataFrame(rows)


def run_manuscript_table2(output: str | Path, n_seeds: int = MC_SEEDS) -> tuple[pd.DataFrame, pd.DataFrame]:
    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)
    scene = make_blind_scene(GRID_NY, GRID_NX)
    geom = build_scene_geometry(scene)
    truth = np.column_stack(
        [
            scene["truth"]["E"].reshape(-1),
            scene["truth"]["N"].reshape(-1),
            scene["truth"]["U"].reshape(-1),
        ]
    )
    rows = []
    diagnostics = []
    for i in range(n_seeds):
        seed = REPRESENTATIVE_SEED + i
        estimates, diag = solve_one_seed(scene, geom, seed)
        diagnostics.append(diag)
        for method, estimate in estimates.items():
            row = {"seed": int(seed), "Method": method}
            row.update(compute_enu_rmse(truth, estimate))
            rows.append(row)
    detail = pd.DataFrame(rows)
    summary = summarize(detail)
    summary_display = summary[["Method", "East RMSE (mm/yr)", "North RMSE (mm/yr)", "Up RMSE (mm/yr)"]]
    detail.to_csv(out_dir / "Table_2_maintext_ENU_RMSE_per_seed_mm_yr.csv", index=False)
    summary.to_csv(out_dir / "Table_2_maintext_ENU_RMSE_numeric_mean_sd_mm_yr.csv", index=False)
    summary_display.to_csv(out_dir / "Table_2_maintext_ENU_RMSE_mean_pm_sd_mm_yr.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(out_dir / "Table_2_maintext_cdir_edge_diagnostics.csv", index=False)
    metadata = {
        "grid_ny": GRID_NY,
        "grid_nx": GRID_NX,
        "pixel_size_m": PIXEL_SIZE_M,
        "monte_carlo_realizations": int(n_seeds),
        "seed_start": REPRESENTATIVE_SEED,
        "seed_end": REPRESENTATIVE_SEED + n_seeds - 1,
        "los_vectors": LOS_VECTORS.tolist(),
        "sigma2_track": SIGMA2_TRACK.tolist(),
        "adaptive_gamma": ADAPTIVE_GAMMA,
        "lambda_s": LAMBDA_S,
        "lambda_d_base": LAMBDA_D_BASE.tolist(),
        "lambda_t_base": LAMBDA_T_BASE.tolist(),
        "edge_kernel": "legacy_linear",
    }
    (out_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return detail, summary

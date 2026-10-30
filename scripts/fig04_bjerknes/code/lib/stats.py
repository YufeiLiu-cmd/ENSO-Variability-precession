"""Shared climate-index statistics (anomaly, OLS, 95% CI)."""
from __future__ import annotations

import numpy as np
from scipy import stats


def monthly_anomaly(ts: np.ndarray) -> np.ndarray:
    """Remove 12-month climatology (matches original MATLAB monthlyAnomaly)."""
    ts = np.asarray(ts, dtype=float).ravel()
    n = ts.size
    clim = np.full(12, np.nan)
    for m in range(12):
        clim[m] = np.nanmean(ts[m::12])
    month_idx = np.arange(n) % 12
    return ts - clim[month_idx]


def area_weights(lat_sub: np.ndarray, n_lon: int) -> np.ndarray:
    """cos(lat) area weights, shape (n_lon, n_lat), normalized to sum=1."""
    w_lat = np.cos(np.deg2rad(np.asarray(lat_sub, float).ravel()))
    w2d = np.ones((n_lon, 1)) * w_lat.reshape(1, -1)
    w2d = w2d / np.nansum(w2d)
    return w2d


def regional_mean_llt(field_tll: np.ndarray, lat_idx: np.ndarray, lon_idx: np.ndarray, w2d: np.ndarray) -> np.ndarray:
    """Area-weighted regional mean.

    field_tll: (time, lat, lon)  — Python netCDF order
    w2d: (n_lon_sub, n_lat_sub) matching lon_idx × lat_idx ordering used in MATLAB
         where MATLAB used field(lon, lat, time) * W(lon, lat)
    Returns time series (n_time,).
    """
    sub = field_tll[:, lat_idx, :][:, :, lon_idx]  # time, lat_sub, lon_sub
    # w2d is (n_lon, n_lat); transpose to (lat, lon) for broadcasting
    w_lt = w2d.T  # (n_lat, n_lon)
    # Match MATLAB: squeeze(sum(sum(field.*W,1,'omitnan'),2,'omitnan'))
    # with W already normalized to sum=1 (no re-normalization for NaNs).
    return np.nansum(sub * w_lt[None, :, :], axis=(1, 2))


def ols_slope_through_origin(x: np.ndarray, y: np.ndarray, ci_level: float = 0.95):
    """Zero-intercept OLS: y = μ x  (MATLAB used cov/var form).

    Returns μ, (ci_lo, ci_hi), r, n.
    SE matches original: se = sqrt(sse/(n-1)) / sqrt(sum(x^2))
    """
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    valid = ~(np.isnan(x) | np.isnan(y))
    x = x[valid]
    y = y[valid]
    n = int(x.size)
    if n < 3 or np.var(x, ddof=0) < np.finfo(float).eps:
        return np.nan, (np.nan, np.nan), np.nan, n

    # MATLAB: mu = cov(x,y)/var(x) with default N-1 normalization → same as
    # sum((x-meanx)(y-meany))/sum((x-meanx)^2); but for zero-intercept they used
    # cov/var which is WITH demeaning. Match exactly:
    c = np.cov(x, y, ddof=1)
    mu = c[0, 1] / np.var(x, ddof=1)
    r = float(np.corrcoef(x, y)[0, 1])

    resid = y - mu * x
    sse = np.sum(resid**2)
    se = np.sqrt(sse / (n - 1)) / np.sqrt(np.sum(x**2))
    t_crit = float(stats.t.ppf(1 - (1 - ci_level) / 2, n - 1))
    return float(mu), (mu - t_crit * se, mu + t_crit * se), r, n


def ols_slope_with_intercept(x: np.ndarray, y: np.ndarray, ci_level: float = 0.95):
    """OLS with intercept: y = a + β x  (used for thermocline α / ocean β).

    Matches MATLAB thermoclineFeedbackRegress:
      alpha = cov(T,H)/var(H)
      se = sqrt(sse/(n-1)) / sqrt(sum((H-mean(H))^2))
    which is the standard demeaned slope (intercept absorbed by demeaning).
    """
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    valid = ~(np.isnan(x) | np.isnan(y))
    x = x[valid]
    y = y[valid]
    n = int(x.size)
    if n < 3 or np.var(x, ddof=0) < np.finfo(float).eps:
        return np.nan, (np.nan, np.nan), np.nan, n

    c = np.cov(np.column_stack([y, x]), rowvar=False, ddof=1)
    # cov([T, H]) in MATLAB → C(1,2)/var(H); here y=T, x=H
    beta = c[0, 1] / np.var(x, ddof=1)
    r = float(np.corrcoef(x, y)[0, 1])
    resid = y - beta * x
    sse = np.sum(resid**2)
    se = np.sqrt(sse / (n - 1)) / np.sqrt(np.sum((x - np.mean(x)) ** 2))
    t_crit = float(stats.t.ppf(1 - (1 - ci_level) / 2, n - 1))
    return float(beta), (beta - t_crit * se, beta + t_crit * se), r, n


def align_lag(driver: np.ndarray, response: np.ndarray, lag: int):
    """Align series so response is shifted by +lag relative to driver.

    lag > 0: response(t+lag) vs driver(t)  — driver leads
    lag < 0: response(t) vs driver(t+|lag|) — response leads
    lag = 0: contemporaneous
    """
    d = np.asarray(driver, float).ravel()
    r = np.asarray(response, float).ravel()
    n = d.size
    if abs(lag) >= n:
        return np.array([]), np.array([])
    if lag >= 0:
        return d[: n - lag], r[lag:]
    return d[-lag:], r[: n + lag]


def product_ci_delta(mu, mu_ci, beta, beta_ci, alpha, alpha_ci):
    """Propagate 95% CI for BJ = μ·β·α via relative delta method (asymmetric → avg half)."""
    bj = mu * beta * alpha
    def half(v, ci):
        return 0.5 * (abs(v - ci[0]) + abs(ci[1] - v))

    h_mu = half(mu, mu_ci)
    h_b = half(beta, beta_ci)
    h_a = half(alpha, alpha_ci)
    if not np.isfinite(bj) or abs(bj) < 1e-30:
        return bj, (np.nan, np.nan)
    rel = np.sqrt((h_mu / (abs(mu) + 1e-30)) ** 2 + (h_b / (abs(beta) + 1e-30)) ** 2 + (h_a / (abs(alpha) + 1e-30)) ** 2)
    h = abs(bj) * rel
    return float(bj), (bj - h, bj + h)

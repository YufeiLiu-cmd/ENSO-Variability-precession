"""Time-series preprocessing and confidence intervals (Scheme 1).

Mean state (unchanged):
  3-month running mean -> time mean

95% CI of the mean (climate-appropriate):
  3-month running mean -> calendar-year means (n = 100 for 1200 months)
  -> Student-t CI of the mean of annual means

Rationale: monthly residuals after deseason/detrend remain strongly
autocorrelated; using n≈1200 months in s/√n underestimates uncertainty
(~0.5–1% of the mean). Annual means better match interannual climate
degrees of freedom for mean-state figures. No ad-hoc rescaling.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def running_mean_3month(x: np.ndarray, center: bool = True) -> np.ndarray:
    """Centered 3-month running mean; ends use available neighbors (min_periods=1)."""
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    out = np.full(n, np.nan)
    if center:
        for i in range(n):
            i0 = max(0, i - 1)
            i1 = min(n, i + 2)
            seg = x[i0:i1]
            if np.any(np.isfinite(seg)):
                out[i] = np.nanmean(seg)
    else:
        for i in range(n):
            i0 = max(0, i - 2)
            seg = x[i0 : i + 1]
            if np.any(np.isfinite(seg)):
                out[i] = np.nanmean(seg)
    return out


def remove_seasonal_cycle(x: np.ndarray, n_month: int = 12) -> np.ndarray:
    """Subtract month-of-year climatology (MATLAB-compatible month-first reshape)."""
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n % n_month != 0:
        clim = np.array(
            [np.nanmean(x[m::n_month]) for m in range(n_month)], dtype=float
        )
        return x - clim[np.arange(n) % n_month]
    mat = np.reshape(x, (n_month, n // n_month), order="F")
    clim = np.nanmean(mat, axis=1)
    anom = mat - clim[:, None]
    return anom.ravel(order="F")


def linear_detrend(x: np.ndarray) -> np.ndarray:
    """Remove linear least-squares trend; NaNs ignored via finite mask."""
    x = np.asarray(x, dtype=float).ravel().copy()
    t = np.arange(x.size, dtype=float)
    m = np.isfinite(x)
    if m.sum() < 3:
        return x
    coef = np.polyfit(t[m], x[m], 1)
    x[m] = x[m] - np.polyval(coef, t[m])
    return x


def annual_means_from_monthly(x: np.ndarray, n_month: int = 12) -> np.ndarray:
    """Calendar-year means from a monthly series (complete years only)."""
    x = np.asarray(x, dtype=float).ravel()
    n_year = x.size // n_month
    if n_year < 2:
        raise ValueError(f"Need ≥2 complete years for annual CI (got {n_year})")
    block = x[: n_year * n_month].reshape(n_year, n_month)
    return np.nanmean(block, axis=1)


def scheme1_mean_and_ci(
    x: np.ndarray, ci_level: float = 0.95
) -> tuple[float, float, float, dict]:
    """Return mean_state, ci_low, ci_high, diagnostics.

    Scheme 1:
      smoothed = 3-month running mean(x)
      mean_state = nanmean(smoothed)
      annual = year means of smoothed (n≈100)
      CI = mean_state +/- tcrit * std(annual)/sqrt(N_year)
    """
    x = np.asarray(x, dtype=float).ravel()
    smoothed = running_mean_3month(x)
    mean_state = float(np.nanmean(smoothed))

    annual = annual_means_from_monthly(smoothed)
    valid = annual[np.isfinite(annual)]
    n = valid.size
    if n < 3:
        raise ValueError(f"Not enough annual samples for CI (n={n})")

    sigma = float(np.std(valid, ddof=1))
    sem = sigma / np.sqrt(n)
    alpha = 1.0 - ci_level
    tcrit = float(stats.t.ppf(1.0 - alpha / 2.0, df=n - 1))
    half = tcrit * sem
    diag = {
        "n": n,
        "n_month": int(np.isfinite(smoothed).sum()),
        "sigma_annual": sigma,
        "sem": sem,
        "tcrit": tcrit,
        "half_width": half,
        "smoothed": smoothed,
        "annual": annual,
        "ci_basis": "annual_means_of_3mo_smoothed",
    }
    return mean_state, mean_state - half, mean_state + half, diag


def lon_to_360(lon: np.ndarray) -> np.ndarray:
    """Convert longitude to [0, 360). Works for -180..180 or 0..360."""
    lon = np.asarray(lon, dtype=float)
    return np.mod(lon, 360.0)


def sort_lon_360(lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted lon360 and argsort indices."""
    lon360 = lon_to_360(lon)
    order = np.argsort(lon360)
    return lon360[order], order

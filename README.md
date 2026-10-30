# Precession-forced asymmetric continental heating shapes ENSO variability

Analysis and plotting scripts for the *Climate of the Past* manuscript of the same title. The repository contains the scripts used to calculate and plot the published diagnostics. It does not contain model output or proxy time series.

## Overview

The scripts diagnose ENSO variability, the zonal position of the Indo-Pacific convective region, moist-static-energy EOFs, tropical Pacific mean-state indices, and the three-term Bjerknes loop gain across 24 AWI-ESM precession time slices. Orbital experiments themselves are not run here. Each script reads files that must be obtained separately and writes figures or intermediate results locally.

Set the input location once:

```bash
export ENSO_DATA_ROOT=/path/to/your/data
```

If `ENSO_DATA_ROOT` is unset, MATLAB scripts look in `data/` two levels above the script, and Python scripts look in `data/` at the repository root.

## Directory structure

```text
ENSO-Variability-precession/
├── README.md
├── requirements.txt
├── .gitignore
└── scripts/
    ├── fig01_equatorial_ssta/       Figure 1a, equatorial SSTA standard deviation
    ├── fig01_itcz_centroid/         Figure 1b, precipitation-weighted longitude
    ├── fig01_plot/                  Figure 1 plotting notebook
    ├── fig03_mse/                   column MSE and seasonal/annual EOF
    ├── fig03_plot/                  Figure 3 plotting notebook
    ├── fig04_nino3_amplitude/       Niño3 amplitude used in Figure 4a
    ├── fig04_mean_state/code/       Figure 4a–h mean-state diagnostics and final panel
    ├── fig04_bjerknes/code/         Figure 4i–l feedback regressions
    ├── figS2_low_eccentricity/      Supplementary Figure S2
    ├── figS3_enso_events/           Supplementary Figure S3
    ├── figS4_climatology/           Supplementary Figure S4
    ├── figS5_itcz_wind/             Supplementary Figure S5
    ├── figS6_walker/                Supplementary Figure S6
    ├── figS7_thermocline/           Supplementary Figure S7
    └── figS8_precip_anomaly/        Supplementary Figure S8
```

## Data availability

Download the model output separately. Do not commit it to this repository.

- Model output: [https://doi.org/10.5281/zenodo.13681177](https://doi.org/10.5281/zenodo.13681177)
- Proxy archives: [https://www.pangaea.de](https://www.pangaea.de) and [https://www.ncei.noaa.gov/products/paleoclimatology](https://www.ncei.noaa.gov/products/paleoclimatology)

The Zenodo record is the model-output archive cited by the manuscript. This repository does not redistribute those files, and it does not redistribute proxy series. Site lists compiled for the paper are in Supplementary Tables S1 and S2. Public access to a database does not by itself grant permission to copy individual proxy series into this repository.

Place inputs under `$ENSO_DATA_ROOT` using the names the scripts open:

| Location | Used by | Contents |
|---|---|---|
| `tsw/` | Niño3 amplitude and Supplementary Figure S3 | `000tsw.nc` … `345tsw.nc`, step 15°, variable `tsw` |
| `MSEgeopoth/`, `aps/` | column MSE | `000forMSE.nc` … and `000apsforMSE.nc` … |
| `fig04_mean_state/walker/` | Walker center | `psi24.mat` |
| `fig04_mean_state/ustrw/nc/` | equatorial zonal wind stress and WWBs | wind-stress NetCDF |
| `fig04_mean_state/thetao/nc/` | MLD and stratification | `thetao_000_300m.nc` … `thetao_345_300m.nc` |
| `fig04_mean_state/ipwp/` | Figure 4d | `exp24aprczhixin.mat` |
| `fig04_bjerknes/tsw/nc/`, `ustrw/nc/`, `thetao/nc/`, `prep/` | feedback regressions | experiment files plus `z20_eq_lon_time_exp.nc` and `sst_eq_lon_time_exp.nc` |
| `figS7/` | Supplementary Figure S7 plot | `thermocline_diff.mat` and `cmp_b2r.rgb` |

Other plotting notebooks open additional `.mat` files by basename from `$ENSO_DATA_ROOT`. Colormap `.rgb` files are not included; put any colormap you are allowed to use in `$ENSO_DATA_ROOT/colormaps/`.

## Software requirements

Verified by import or by a syntax check on this machine:

- Python 3, with `numpy` 2.5.1, `scipy` 1.18.0, and `matplotlib` 3.11.0 importable during the syntax check
- `netCDF4`, `pandas`, `statsmodels`, and `colormaps`, imported by notebooks; they were not installed in that interpreter
- `basemap` and `cartopy`, imported by `scripts/figS4_climatology/plot_figS4.ipynb`; not installed in that interpreter

Required by scripts but not available to the syntax check:

- MATLAB, including the Mapping Toolbox calls `m_proj` and `m_contourf` (M_Map). M_Map is not bundled.
- NCL, for `scripts/figS6_walker/walkstudygai.ncl`
- CDO with the ECHAM `afterburner` operator, for `scripts/figS6_walker/cdo_afterburner_uv.sh`. The script calls `module load cdo`, which applies on a module-based HPC system.
- `caleof` from the PCA toolbox written by Guillaume Maze (2004). It is not bundled. Place it in `$ENSO_DATA_ROOT/external/pcatool` before running the EOF scripts. `EOF_hyang.m` calls `caleof(M,3,1)`.

No package versions beyond the three imported above were recorded, because the original analysis environment was not stored with the scripts.

## Usage

From the repository root, after `ENSO_DATA_ROOT` is set and the input files are in place:

```bash
# Figure 4a Niño3 amplitude (MATLAB)
# In MATLAB:
#   cd scripts/fig04_nino3_amplitude
#   compute_nino3_enso

# Figure 4i–l, then Figure 4a–h and the two-column figure
python3 scripts/fig04_bjerknes/code/run_all.py
python3 scripts/fig04_mean_state/code/run_all.py

# Supplementary Figure S3
python3 scripts/figS3_enso_events/enso_event_analysis.py
python3 scripts/figS3_enso_events/compose_fig_S3.py

# Supplementary Figure S7 panel, after thermocline_section.m has written its section
python3 scripts/figS7_thermocline/compose_fig_S7.py
```

Run the MATLAB and NCL scripts from their own directories so relative `data_root` resolution and `EOF_hyang.m` are on the path. Open the `.ipynb` files in the same environment, with `ENSO_DATA_ROOT` set. Figure outputs from the Python drivers are written under `scripts/fig04_mean_state/result/figures/` and `scripts/fig04_bjerknes/result/figures/`.

`scripts/fig04_mean_state/code/09_combined_meanstate_feedbacks_2col.py` reads the `.mat` files produced by the other Figure 4 scripts. Run the Bjerknes driver before the mean-state driver.

## Figure mapping

| Manuscript item | Script |
|---|---|
| Figure 1a | `scripts/fig01_equatorial_ssta/equatorial_ssta_std.m`, plotted by `scripts/fig01_plot/plot_fig1.ipynb` |
| Figure 1b | `scripts/fig01_itcz_centroid/itcz_precipitation_centroid.m`, plotted by the same notebook |
| Figure 2 | No analysis script. The directory that held the schematic contained no code. |
| Figure 3 | `scripts/fig03_mse/msecalculate.m`, then `eof_jja.m`, `eof_djf.m`, and `eof_annual.m`, then `interpolate_eof_fields.m`; plotted by `scripts/fig03_plot/plot_fig3.ipynb` |
| Figure 4a | `scripts/fig04_nino3_amplitude/compute_nino3_enso.m` |
| Figure 4b | `scripts/fig04_mean_state/code/01_walker_edges.py` |
| Figure 4c | `scripts/fig04_mean_state/code/02_ustrw.py` |
| Figure 4d | `scripts/fig04_mean_state/code/03_ipwp_plot.py` (reads the centroid file; it does not recompute it) |
| Figure 4e | `scripts/fig04_mean_state/code/04_mld.py` |
| Figure 4f | `scripts/fig04_mean_state/code/05_stratification.py` |
| Figure 4g | `scripts/fig04_mean_state/code/07_wwb.py` |
| Figure 4h | `scripts/fig04_mean_state/code/08_wwv.py` |
| Figure 4i | `scripts/fig04_bjerknes/code/01_wind_feedback.py` |
| Figure 4j | `scripts/fig04_bjerknes/code/03_ocean_response_beta.py` |
| Figure 4k | `scripts/fig04_bjerknes/code/02_thermocline_feedback.py` |
| Figure 4l | `scripts/fig04_bjerknes/code/04_bjerknes_index.py` |
| Figure 4, assembled | `scripts/fig04_mean_state/code/09_combined_meanstate_feedbacks_2col.py` |
| Figure 5 | Not included. The plotting notebook embeds proxy values taken from published records. |
| Supplementary Figure S1 | Not included. The only mesh-interpolation program found is a third-party FESOM utility and contains paths that should not be published. |
| Supplementary Figure S2 | `scripts/figS2_low_eccentricity/equatorial_ssta_std_ecc002.m` and `plot_figS2.ipynb` |
| Supplementary Figure S3 | `scripts/figS3_enso_events/enso_event_analysis.py` and `compose_fig_S3.py` |
| Supplementary Figure S4 | `scripts/figS4_climatology/climatology_sst_precip.m` and `plot_figS4.ipynb` |
| Supplementary Figure S5 | `scripts/figS5_itcz_wind/plot_figS5.ipynb` |
| Supplementary Figure S6 | `scripts/figS6_walker/cdo_afterburner_uv.sh`, `walkstudygai.ncl`, `walk_threeocean.m`, and `plot_figS6.ipynb` |
| Supplementary Figure S7 | `scripts/figS7_thermocline/thermocline_section.m` and `compose_fig_S7.py` |
| Supplementary Figure S8 | `scripts/figS8_precip_anomaly/precip_anomaly_90_270.m` and `plot_figS8.ipynb` |
| Supplementary Figures S9 and S10 | Not included. They plot proxy series and a proxy–model regression. |

The manuscript states that lead–lag correlations and El Niño composites were supplied to reviewers and are not part of the paper. Those scripts were not copied. Sensitivity figures for Niño3.4 and Niño4 are described as not shown and were not copied.

## Method notes

These points are differences between the archived scripts and the manuscript text. The calculations were not edited to remove them.

- Figure 1a uses 2400 months and does not remove a linear trend. The low-eccentricity Supplementary Figure S2 script uses 1200 months, longitude 120–280°E, and latitude 5°S–5°N, and also does not detrend. The methods text says the final 100 years, Niño3, and linear detrending.
- The precipitation centroid uses a threshold of 4.93 mm day⁻¹ after converting convective precipitation with a factor of 86400. The manuscript states 5 mm day⁻¹. The saved script reads `aprc` from a file whose name is a FESOM temperature file; confirm the precipitation file that was actually used.
- Figure 4a, from `compute_nino3_enso.m`, is the standard deviation of deseasoned, linearly detrended monthly Niño3 SST. It does not apply the three-month running mean, and the area mean is not weighted by the cosine of latitude. Supplementary Figure S3 does use a cosine-weighted mean and a three-month running mean before the standard deviation, event count, and spectrum.
- Figure 4d error bars are `1.96 * std / sqrt(1200)` from the centroid script. Wind stress, mixed-layer depth, and stratification use a Student-t interval of annual means.
- Walker center: mass streamfunction averaged over 200–500 hPa; positive-ψ weighted longitude over 100–260°E.
- Equatorial τx, WWB, WWV, and stratification: 5°S–5°N, 140°E–150°W. Mixed-layer depth uses the same longitudes but 2.5°S–2.5°N and a 0.5 K threshold.
- WWB months are months with deseasoned, detrended τx′ > 0.015 Pa.
- Wind feedback μ is the zero-lag regression of western Pacific τx′ (140–180°E) on Niño3 SST′.
- Ocean response β is the largest regression of eastern Pacific Z20′ (150°W–90°W) on western Pacific τx′ over lags 0–6 months, selecting the most positive slope.
- Thermocline feedback α is the largest regression of Niño3 SST′ on Z20′ over 140–260°E for lags 1–6 months, with Z20 leading.
- BJ = μβα, with a relative delta-method confidence interval. This is the three-term product described in the Supplementary Information, not the full Jin et al. (2006) index.
- The EOF driver sets `MON=1` before calling `caleof`. For these calls the third dimension is the 24 precession phases, not a monthly time axis.
- Column MSE uses Cp = 1005 J kg⁻¹ K⁻¹, g = 9.8 m s⁻², and Lv = 2.26×10⁶ J kg⁻¹, integrated with layer pressure thicknesses and divided by surface pressure, then converted to kJ kg⁻¹.
- Supplementary Figures S6 and S7 read files named `105` and `285`. The manuscript labels the corresponding states 90° and 270°.

## Citation

Liu, Y., Yang, H., Dang, H., Lohmann, G., Lian, T., Lu, Z., Shi, X., Liu, J., Wang, H., Jiang, F., Deng, X., Hu, X., and Chen, D.: Precession-forced asymmetric continental heating shapes ENSO variability, *Climate of the Past*.

Corresponding author: Hu Yang (yanghu@sml-zhuhai.cn).

No software license was found with these scripts. None has been added. Until the authors choose a license, copyright remains with the authors and the code is not licensed for reuse by others.

## License

No license is specified in the original scripts. This repository does not add one.

# Algorithmic Reference & Scientific Foundations

The International Reference Ionosphere (IRI) is the internationally recognized standard empirical model of the Earth's ionosphere, sponsored by the Committee on Space Research (COSPAR) and the International Union of Radio Science (URSI). This document outlines the scientific, empirical, and mathematical models implemented in the pure Rust codebase of `iri2020-rust`.

---

## 1. Core Simulation Engine
The main execution flow is controlled by the `irisub.rs` module, translating the legacy `IRI_SUB` control routine.

- **Physical Quantities**: Computes electron density ($N_e$), electron temperature ($T_e$), ion temperature ($T_i$), neutral temperature ($T_n$), and relative percentage abundances of positive ions ($O^+$, $H^+$, $He^+$, $O_2^+$, $NO^+$).
- **Execution Settings**: Evaluates 50 user-toggleable control flags (`jf`) that specify model options, index overrides, and peak density/height selections.
- **Reference**: [Bilitza et al., 2022 (Journal of Space Weather and Space Climate)](https://doi.org/10.1051/swsc/2022009)

## 2. Total Electron Content (TEC) Integration
The module `iritec.rs` computes the integrated electron column density along a vertical path from a lower altitude boundary ($h_{start}$) to an upper altitude boundary ($h_{end}$):

$$\text{vTEC} = \int_{h_{start}}^{h_{end}} N_e(h) \, dh$$

- **Numerical Quadrature**: Employs a multi-segment trapezoidal integration method with 1000 sub-steps per segment.
- **Segment Matching**: Features precise segment-end state recovery. In compliance with original empirical constraints, the integrated topside ($TEC_{top}$) and bottomside ($TEC_{bot}$) bounds match the F2-peak options to preserve output consistency.
- **Reference**: [Bilitza, 2001 (Radio Science)](https://doi.org/10.1029/2000RS002432)

## 3. Geomagnetic Field Model (IGRF-13)
The module `igrf.rs` and static data in `igrf_coeff.rs` implement the 13th generation of the **International Geomagnetic Reference Field (IGRF)**.

- **Mathematical Form**: Evaluates the main geomagnetic field potential ($V$) as a spherical harmonic expansion of degree $N=13$:

$$V(r, \theta, \phi) = a \sum_{n=1}^{N} \sum_{m=0}^{n} \left(\frac{a}{r}\right)^{n+1} \left[ g_n^m \cos(m\phi) + h_n^m \sin(m\phi) \right] P_n^m(\cos\theta)$$

  where $a$ is the Earth's mean radius, $r$ is the radial distance, $\theta$ is the geocentric co-latitude, $\phi$ is the east longitude, $g_n^m$ and $h_n^m$ are Gauss coefficients, and $P_n^m$ represents Schmidt quasi-normalized Legendre functions.

- **Temporal Extrapolation**: Interpolates the Gauss coefficients linearly between five-year epochs and applies secular variation models for predictive intervals beyond the latest epoch.
- **Reference**: [Alken et al., 2021 (Earth, Planets and Space)](https://doi.org/10.1186/s40623-020-01288-x)

## 4. Neutral Atmosphere Model (CIRA-86 / NRLMSISE-00)
The module `cira.rs` and coefficients in `cira_coeff.rs` implement the COSPAR International Reference Atmosphere model.

- **Scientific Foundation**: NRLMSISE-00 models the neutral temperature profile and the number densities of major/minor thermospheric species ($He$, $O$, $N_2$, $O_2$, $Ar$, $H$, $N$, and anomalous $O$).
- **Numerical Formulation**: Utilizes cubic spline interpolations and spherical harmonic coefficients mapped globally to calculate the diurnal, semidiurnal, and seasonal variations of neutral densities.
- **Reference**: [Picone et al., 2002 (Journal of Geophysical Research)](https://doi.org/10.1029/2002JA009430)

## 5. ROCSAT-1 Vertical Plasma Drift
Equatorial vertical ion drift velocities in the F-region are calculated in `rocdrift.rs` using the 64,900 coefficient empirical tensor stored in `rocdrift_coeff.rs`.

- **Physical Model**: Represents empirical F-region vertical drifts measured by the ROCSAT-1 satellite.
- **Method**: Evaluates drifts based on local time, longitude, and solar activity ($F_{10.7}$ index) using multidimensional bilinear interpolations.
- **Reference**: [Fejer et al., 2008 (Journal of Geophysical Research)](https://doi.org/10.1029/2008JA013177)

## 6. F2-Layer Bottomside Thickness ($B_0$, $B_1$)
The thickness parameter $B_0$ (which scales the profile width below the F2 peak) and the shape parameter $B_1$ (which controls the profile slope) are computed in `b0_b1_model.rs`.

- **Scientific Approach**: Offers the Altadill spherical harmonics model option (derived from global ionosonde records) alongside traditional Bilitz/Gulyaeva table options.
- **Reference**: [Altadill et al., 2009 (Journal of Atmospheric and Solar-Terrestrial Physics)](https://doi.org/10.1016/j.jastp.2008.09.043)

## 7. Electron Density Profile Construction
The vertical electron density profile $N_e(h)$ is synthesized in `xe_profile.rs` by dividing the ionosphere into structural layers:

- **Profile Layers**: Bottomside D-region, E-region, E-F valley, F1-region, and F2-region.
- **Formulation**: Utilizes Rawer's Epstein LAY-functions to construct smooth, continuous analytical profiles:

$$\text{LAY}(h; h_i, W_i) = \ln\left(1 + e^{\frac{h-h_i}{W_i}}\right) - \ln\left(1 + e^{\frac{h_m-h_i}{W_i}}\right) - \frac{h-h_m}{W_i \left(1 + e^{\frac{h_m-h_i}{W_i}}\right)}$$

  where $h_i$ is the transition height and $W_i$ is the transition width parameter.

- **Reference**: [Rawer, 1988 (Advances in Space Research)](https://doi.org/10.1016/0273-1177(88)90234-5)

## 8. D-Region Electron Density Models
Implements lower ionosphere models (altitudes 60 to 90 km) in `d_region.rs` and `iridreg.rs`.

- **FIRI Model**: The Faraday International Reference Ionosphere (stored as a massive 160,380 float lookup block in `firi_data.rs`) represents Friedrich and Torkar's empirical lower-ionosphere profiles, evaluated via 5D interpolation.
- **Danilov Model**: Implements Danilov's D-region model representing variations in ion composition and electron density under quiet and disturbed conditions.
- **References**: [Friedrich & Torkar, 2001 (Journal of Geophysical Research)](https://doi.org/10.1029/2001JA900070); [Danilov et al., 1995 (Advances in Space Research)](https://doi.org/10.1016/0273-1177(94)00071-G)

## 9. Ion Composition Model
The positive ion density profiles ($O^+$, $H^+$, $He^+$, $O_2^+$, $NO^+$) are computed in `ioncom.rs` with coefficients defined in `calion_coeff.rs`.

- **Physical Formulation**: Solves the altitude transitions where atomic oxygen ions dominate the F2 region, transitioning to molecular species ($O_2^+$, $NO^+$) in the E region, and light ions ($H^+$, $He^+$) in the topside plasmasphere.
- **Reference**: [Bilitza, 1990 (Advances in Space Research)](https://doi.org/10.1016/0273-1177(90)90184-R)

## 10. Coordinate Systems & Field-Line Tracing (GEO-CGM)
The module `iriflip.rs` implements geomagnetic coordinate transformations and magnetic field line tracing using spatial arrays in `cormag_data.rs`.

- **Coordinate Conversions**: Converts between Geocentric coordinates and Corrected Geomagnetic (CGM) coordinates.
- **Numerical Integrator**: Employs a 4th-order Runge-Kutta integration scheme with adaptive step scaling to trace magnetic field lines from geographic points to their conjugate points or F-region magnetic footprints.
- **Reference**: [Gustafsson et al., 1992 (Journal of Atmospheric and Terrestrial Physics)](https://doi.org/10.1016/0021-9169(92)90109-7)

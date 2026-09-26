# Step 6M — Final Project Synthesis

No additional FDTD simulation is required. This step freezes the validated results from Steps 1–6L into a portfolio-ready summary.

## Preferred full-device result

Using the 25-nm subpixel representation:

- focal position: **25.6177 µm**
- nominal focus: **28.0 µm**
- focal shift: **-2.3823 µm**
- lateral FWHM: **1.3752 µm**
- post-lens transmission: **81.32%**
- incident-power-normalized 3×FWHM focusing efficiency: **64.83%**
- 3×FWHM capture fraction: **87.49%**

The 33.333→25 nm subpixel change is only **1.11%** in eta_3×FWHM and **0.27%** in post-lens transmission.

## Spectral result

The broadband robustness study uses the validated 33.333-nm subpixel model for practical cost. Its 1550-nm spectral reference is eta_3×FWHM = **0.638485**, giving an 80%-retention threshold of **0.510788**.

Dense threshold localization gives:

- lower edge: **1427.50 nm**
- first upper downward edge: **2390.51 nm**
- continuous 1550-connected 80%-retention bandwidth: **963.01 nm**
- recovery edge: **2434.56 nm**
- below-threshold notch width: **44.05 nm**

Because the response recovers after the notch, the later pass region must be reported separately rather than merged into the continuous 1550-centered band.

## What can be claimed

A defensible portfolio statement is:

> Developed a custom C++ 2D TMz FDTD framework for discrete metalens design, including CPML, complex Poynting-flux monitors, multi-frequency Fourier analysis, OpenMP acceleration, and subpixel dielectric averaging. Designed and numerically validated an eight-state 1550-nm metalens library and a 21-cell finite metalens. The preferred 25-nm subpixel simulation achieved 64.8% incident-power-normalized focusing efficiency within a 3×FWHM focal window, with a 1.38-µm lateral FWHM. A broadband study of the idealized nondispersive model identified a continuous 80%-retention band from about 1.43 to 2.39 µm.

## Important scope

Do not present these values as an experimental or full-3D device prediction. The model is:

- 2D TMz, per-unit-out-of-plane length
- nondispersive and lossless n = 2
- no substrate
- no material absorption
- no fabrication roughness/tolerance
- no full 3D polarization or coupling study

The very broad spectral interval is therefore a numerical property of this idealized model.

## Project status

**Portfolio numerical project complete.**

Further work such as real-material dispersion, substrate inclusion, tolerance analysis, or 3D FDTD should be treated as a new project extension rather than required verification of this version.

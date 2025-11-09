# MET5510 Project
This repository contains the files, code, and plots for Grace and Matthew's MET5510 project under the guidance of Dr. Cai and Dr. Sun.

## Project Overview
The program simulates various atmospheric wave patterns using the Quasi-Geostrophic Cyclogenesis, Advection, and Instability (QG-CAI) modeling technique for mid-latitude systems. The program computes potential vorticity (PV) gradients, mean zonal wind fields, and related dynamics for the Eady-type and Hoskins-West models. The scripts are designed for wave propagation, stability analysis, and quasi-geostrophic dynamics studies.

## Instructions
### Calculate Model Data
1. `cfg.m`: Sets zonal wave # (m0) and the eMode # (n_mode) for different initial conditions
2. `main.m`: Processes the QG-CAI model's main data, matrices, and eigen vector arrays
3. `calc.m`: Calculates derived values from main data file and iterates with eMode #
4. Generated data files are saved in `data/data_wave-<m0>.mat`
5. Calculated data files are saved in `data/calc_wave-<m0>_emode-<n_mode>.mat`

### Create Initial Plots
1. `plt.m`: Plots all figuresfrom the calculated data from plotting functions in `plots`
2. `plt_one('list')`: Command can be run in terminal to get a list of individual plots to generate
3. `eady_plot.m` to generate initial plots using the Eady data
4. Generated plots are saved in `output/figures/*.png`

### Config Options
1. `hwme_config.m` contains all config, variable, and constants needed for Hoskins-West Modified Eady-type calculations
2. `rossby_config.m` contains all config, variable, and constants needed for Rossby Wave calculations
3. `eady_config.m` contains all the config, variable, and constants needed for Eady calculations
4. Config files are stored in `config/*.m`

## File Structure
### Functions are stored in /functions

#### Tool Scripts
- `plt_one.m`: Can use a command like `plt_one('temp',1)` to get the temperature plot and values at surface

#### Ubar Calculations and other Model Dependant Scripts for QG omega equation
- `ubar.m`: Generates the mean zonal wind field $\bar{u}(y, z)$  in the quasi-geostrophic framework.
- `PV2intgrad.m`: Determines the interior PV gradient for the Eady-type model, including the meridional second derivative of the mean flow.
- `F1matrix.m`, `F2matrix.m`, `F3matrix.m`: Computes the forcing terms F1, F2, and F3 for vertical motion calculations.
- `Gmatrix.m`: Constructs the sparse elliptic operator matrix G (LW×LW) with the result of applying the diagnostic elliptic operator.

#### Utility Scripts
- `bpvy.m`: PVY computes the background potential vorticity gradient (dPV/dy)for Rossby, Eady, or HWME models.
- `eigen.m`: Solves the eigenvalue problem for the quasi-geostrophic model, computing eigenvectors and eigenvalues for matrices.
- `jk2l.m`: Converts 2D indices $(j, k)$ to a single linear index for grid operations, facilitating array manipulation.
- `jk2lw.m`: Transforms 2D indices $(j, k)$ to a linear index with weighting, used for specific grid indexing tasks.
- `l2jk.m`: Converts a linear index back to 2D indices $(j, k)$, supporting grid coordinate recovery.
- `lw2jk.m`: Converts a weighted linear index back to 2D indices $(j, k)$, aiding in weighted grid operations.
- `matrices.m`: Constructs matrices $B$, $C$, and $D$ for PV inversion and advection in the eigenvalue problem, central to stability analysis.
- `stream2pv.m`: Computes potential vorticity from the streamfunction, a key step in quasi-geostrophic dynamics.
- `stream2xPVadv.m`: Calculates the zonal advection of PV, contributing to advection terms in the model.
- `stream2yPVadv.m`: Computes the meridional advection of PV, essential for capturing cross-flow dynamics.
- `w2ellipse.m`: Transforms a vector field to an elliptical representation, used for visualization or analysis.
- `w2wfield.m`: Converts a vector to a 3D field, supporting perturbation field construction.
- `w2vec.m`: Solves the elliptic equation G w = (f0/N²) * (F1 + F2 + F3) for the interior vertical velocity vector.
- `XV2field.m`: Transforms an eigenvector $XV$ to a 3D field, enabling perturbation visualization.
- `XV2streamxtime.m`: Converts an eigenvector to a streamfunction over longitude and time, aiding time-dependent analysis.
- `XV2ugxtime.m`: Computes the time evolution of the zonal wind (ug) along longitude for a Hovmoller diagram at a specific latitude and height.
- `XV2XVx.m`: Computes the x-derivative of an eigenvector, supporting zonal gradient calculations.
- `XV2XVy.m`: Computes the y-derivative of an eigenvector, used for meridional gradient analysis.
- `XV2XVz.m`: Computes the z-derivative of an eigenvector, essential for vertical structure analysis.
- `XVx2field.m`: Converts the x-derivative of an eigenvector to a 3D field, supporting zonal perturbation plots.
- `XVy2field.m`: Converts the y-derivative of an eigenvector to a 3D field, aiding meridional perturbation visualization.
- `XVz2field.m`: Converts the z-derivative of an eigenvector to a 3D field, supporting vertical perturbation analysis.


### Figures are stored in /figures
- `plot_dpvdym_boundaries.m`: Plots d(PVbar)/dy at boundaries with beta.
- `plot_dpvdym_int.m`: Plots d(PVbar)/dy interior contour.
- `plot_evec_amp.m`: Plots eigenvector amplitude contour.
- `plot_gph_hovmoller.m`: Plots Hovmoller diagram for Geopotential Height across latitude-time in days (x-time).
- `plot_gph.m`: Plots geopotential height contour at set hlevels across longitude-latitude (x-y).
- `plot_gph_xsec.m`: Plots geopotential height perturbation at mid-latitude across longitude-height (x-z).
- `plot_meridional_hflux.m`: Plots zonally-averaged meridional eddy heat flux at latitude-height (y-z).
- `plot_meridional_wind.m`: Plots meridional wind perturbation at specified vertical level across longitude-latitude (x-y).
- `plot_meridional_xsec.m`: Plots meridional wind vertical cross-section at mid-latitude across longitude-height (x-z).
- `plot_meridional_flux.m`: Plots zonally-averaged meridional eddy momentum flux <v'u'> across latitude-height (y-z).
- `plot_pvfield.m`: Plots potential vorticity perturbations at given hlevels across longitude-latitude (x-y).

- `plot_hovmoller.m`: Plots Hovmoller diagram
- `plot_zonal_wind.m`: Plots zonal wind contour
- `plot_meridional_wind.m`: Plots meridional wind contour
- `plot_temperature.m`: Plots temperature contour
- `plot_gph_top.m`: Plots geopotential height at top boundary
- `plot_pvfield.m`: Plots potential vorticity contour
- `plot_vg_cross_section.m`: Plots meridional wind vertical cross-section
- `plot_ug_hovmoller.m`: Plots Hovmoller diagram for zonal wind
- `plot_ubar_contour.m`: Plots Ubar contour
- `plot_dpvdym_int.m`: Plots d(PVbar)/dy interior contour
- `plot_dpvdym_boundaries.m`: Plots d(PVbar)/dy at boundaries with beta

#### Example Plots
**Hoskins-West Modified Background Flow**
![hwm_background_flow.png](https://i.imgur.com/Eoe0sqx.png)

**Hoskins-West Eady-type Background Flow**
![hwm_background_flow.png](https://i.imgur.com/JfLbQt6.png)

**Rossby Wave Eigenvector Amplitude**
![row_eigenvector_amplitude.png](https://i.imgur.com/58LJSsK.png)

**Rossby Wave Geopotential Height Hovmoller**
![row_hovmoller_geopotential.png](https://i.imgur.com/KZIcWUE.png)

**Rossby Wave Geopotential Height Hovmoller**
![row_meridional_cross_section_geopotential.png](https://i.imgur.com/fsyTOZN.png)

**Rossby Wave Vertical Motion Diagnosis**
![row_vertical_motion_diagnosis.png](https://i.imgur.com/NriIuM2.png)

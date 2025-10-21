# MET5510 Project
This repository contains the files, code, and plots for Grace and Matthew's MET5510 project under the guidance of Dr. Cai and Dr. Sun.

## Project Overview
The program simulates various atmospheric wave patterns using different modeling techniques for mid-latitude systems. Currently, simulations can be run for the **Hoskins-West Modified Eady-tpe** model, the **Rossby Wave** mode, and the **Eady Wave** model. The program computes potential vorticity (PV) gradients, mean zonal wind fields, and related dynamics for the Eady-type and Hoskins-West models. The scripts are designed for wave propagation, stability analysis, and quasi-geostrophic dynamics studies.

### Available Models
- **HWME** - Hoskins-West Eady-type Model
- **Rossby** - Rossby Wave Model
- **Eady** - Eady Wave Model


## Instructions
### Calculate Model Data
1. `hwme_main.m` to process the Hoskins-West Modified Eady-type calculations
2. `rossby_main.m` to process the Rossby Wave calculations
3. `eady_main.m` to process the Eady calculations
4. Generated data files are saved in `output/data/*.mat`

### Create Initial Plots
1. `hwme_plot.m` to generate initial plots using the Hoskins-West Modified data
2. `rossby_plot.m` to generate initial plots using the Rossby Wave data
3. `eady_plot.m` to generate initial plots using the Eady data
4. Generated plots are saved in `output/figures/*.png`

### Config Options
1. `hwme_config.m` contains all config, variable, and constants needed for Hoskins-West Modified Eady-type calculations
2. `rossby_config.m` contains all config, variable, and constants needed for Rossby Wave calculations
3. `eady_config.m` contains all the config, variable, and constants needed for Eady calculations
4. Config files are stored in `config/*.m`

## File Structure
### Functions are stored in /functions

#### Ubar Calculations and other Model Dependant Scripts
- `hwme_ubar.m`: Generates the mean zonal wind field $\bar{u}(y, z)$ for the Hoskins-West Modified Eady-type model in the quasi-geostrophic framework.
- `rossby_ubar.m`: Initializes for the Rossby wave model in the quasi-geostrophic framework, set to a constant value.
- `eady_ubar.m`: Initializes  for the Eady wave model in the quasi-geostrophic framework, based on a linear vertical profile derived from PV gradient.
- `hwme_PV2intgrad.m`: Determines the interior PV gradient for the Eady-type model, including the meridional second derivative of the mean flow.
- `eady_F123.m`: Computes the forcing terms F1, F2, and F3 for vertical motion calculations in the Eady wave model within the quasi-geostrophic framework.

#### Utility Scripts
- `bpvy.m`: PVY computes the background potential vorticity gradient (dPV/dy)for Rossby, Eady, or HWME models.
- `jk2l.m`: Converts 2D indices $(j, k)$ to a single linear index for grid operations, facilitating array manipulation.
- `jk2lw.m`: Transforms 2D indices $(j, k)$ to a linear index with weighting, used for specific grid indexing tasks.
- `l2jk.m`: Converts a linear index back to 2D indices $(j, k)$, supporting grid coordinate recovery.
- `lw2jk.m`: Converts a weighted linear index back to 2D indices $(j, k)$, aiding in weighted grid operations.
- `matricesBCD.m`: Constructs matrices $B$, $C$, and $D$ for PV inversion and advection in the eigenvalue problem, central to stability analysis.
- `stream2pv.m`: Computes potential vorticity from the streamfunction, a key step in quasi-geostrophic dynamics.
- `stream2xPVadv.m`: Calculates the zonal advection of PV, contributing to advection terms in the model.
- `stream2yPVadv.m`: Computes the meridional advection of PV, essential for capturing cross-flow dynamics.
- `w2ellipse.m`: Transforms a vector field to an elliptical representation, used for visualization or analysis.
- `w2wfield.m`: Converts a vector to a 3D field, supporting perturbation field construction.
- `XV2field.m`: Transforms an eigenvector $XV$ to a 3D field, enabling perturbation visualization.
- `XV2streamxtime.m`: Converts an eigenvector to a streamfunction over longitude and time, aiding time-dependent analysis.
- `XV2XVx.m`: Computes the x-derivative of an eigenvector, supporting zonal gradient calculations.
- `XV2XVy.m`: Computes the y-derivative of an eigenvector, used for meridional gradient analysis.
- `XV2XVz.m`: Computes the z-derivative of an eigenvector, essential for vertical structure analysis.
- `XVx2field.m`: Converts the x-derivative of an eigenvector to a 3D field, supporting zonal perturbation plots.
- `XVy2field.m`: Converts the y-derivative of an eigenvector to a 3D field, aiding meridional perturbation visualization.
- `XVz2field.m`: Converts the z-derivative of an eigenvector to a 3D field, supporting vertical perturbation analysis.

### Figures are stored in /figures
- `plot_evec_amp.m`: Plots eigenvector amplitude contour
- `plot_gph.m`: Plots geopotential height contour
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

### Data files are in /data
- `output/data/*.mat` with model and wave number prefixes


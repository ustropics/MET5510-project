# MET5510 Project
This repository contains the files, code, and plots for Grace and Matthew's MET5510 project under the guidance of Dr. Cai and Dr. Sun.

## Project Overview
The program simulates various atmospheric wave patterns using different modeling techniques for mid-latitude systems. Currently, simulations can be run for Rossby waves,the modified Hoskins-West model, and the Hoskins-West Eady-type model. The program computes potential vorticity (PV) gradients, mean zonal wind fields, and related dynamics for the Eady-type and Hoskins-West models. The scripts are designed for wave propagation, stability analysis, and quasi-geostrophic dynamics studies.

### Available Models
- **ROS** - Rossby Wave Model
- **HWM** - Hoskins-West Modified Model
- **HWE** - Hoskins-West Eady-type Model

## Instructions
### Calculate Model Data
1. `row_main.m` to process the Rossby Wave Model calculations
2. `hwm_main.m` to process the Hoskins-West Modified calculations
3. `hwe_main.m` to process the Hoskins-West Eady-type calculations
4. Generated data files are saved in `output/data/*.mat`

### Create Initial Plots
1. `row_plot.m` to generate initial plots using the generated Rossby Wave data
2. `hwm_plot.m` to generate initial plots using the generated Hoskins-West Modified data
3. `hwe_plot.m` to generate initial plots using the generated Hoskins-West Eady-type data
4. Generated plots are saved in `output/plots/*.png`

### Config Options
1. `row_config.m` contains all config, variable, and constants needed for Rossby Wave calculations
2. `hwm_config.m` contains all config, variable, and constants needed for Hoskins-West Modified calculations
3. `hwe_config.m` contains all config, variable, and constants needed for Hoskins-West Eady-type data
4. Config files are stored in `config/*.m`

## File Structure
### Functions are stored in /functions

#### Hoskins-West Eady-Type Model Scripts
- `hwe_BPVyCalc.m`: Computes the meridional PV gradient at interior grid points for the Eady-type model, using a $\cos^4$-modulated mean zonal wind $\bar{u}(y, z)$. Essential for wave propagation studies.
- `hwe_PV2bndgrad.m`: Calculates PV gradients at the surface ($z=0$) and tropopause ($z=HH$) for the Eady-type model, incorporating vertical shear with $\cos^4$ meridional variation. Critical for boundary condition enforcement.
- `hwe_PV2intgrad.m`: Determines the interior PV gradient for the Eady-type model, including the meridional second derivative of the mean flow. Supports stability and dynamics analysis.
- `hwe_ubar.m`: Generates the mean zonal wind field $\bar{u}(y, z)$ for the Eady-type model, using linear shear with $\cos^4$ meridional modulation. Foundational for background flow setup.

#### Hoskins-West Modified Model Scripts
- `hwm_BPVyCalc.m`: Computes the meridional PV gradient for the original Hoskins-West model, using sinusoidal modulation and vertical shear terms. Included for reference.
- `hwm_PV2bndgrad.m`: Calculates PV gradients at boundaries for the Hoskins-West model, incorporating complex shear terms. Retained for historical context.
- `hwm_PV2intgrad.m`: Computes the interior PV gradient for the Hoskins-West model, using planetary $\beta$ and vertical shear effects. Kept for comparison.
- `hwm_ubar.m`: Generates the mean zonal wind field for the Hoskins-West model with sinh-based modulation. Preserved for legacy support.

#### Rossby Wave Model Scripts
- `row_F123.m`: Generates row vectors $F1$, $F2$, and $F3$ for matrix assembly, supporting the linearized PV equation.

#### Utility Scripts
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

### Plots are stored in /plots
- `output/plots/*.png` with model prefixes

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


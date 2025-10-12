# MET5510 Project
This repository contains the files and code for Grace and Matthew's MET5510 project under the guidance of Dr. Cai and Dr. Sun.

## Project Overview
The program simulates various atmospheric wave patterns using different modeling techniques for mid-latitude systems. Currently, simulations can be run for Rossby waves and the modified Hoskins-West model. The program computes potential vorticity (PV) gradients, mean zonal wind fields, and related dynamics for the Eady-type and Hoskins-West models. The scripts are designed for wave propagation, stability analysis, and quasi-geostrophic dynamics studies. The codebase includes both the modern Eady-type model with meridional modulation and legacy Hoskins-West model scripts for reference.

## Instructions
### Calculate Model Data
1. Run `row_main.m` to process the Rossby Wave Model calculations
2. Run `hwm_main.m` to process the Hoskins-West Modified calculations
3. Run `hwe_main.m` to process the Hoskins-West Eady-type calculations
4. Generated data files are saved in `output/data/*.mat`

### Create Initial Plots
1. Run `row_plot.m` to generate initial plots using the generated Rossby Wave data
2. Run `hwm_plot.m` to generate initial plots using the generated Hoskins-West Modified data
3. Run `hwe_plot.m` to generate initial plots using the generated Hoskins-West Eady-type data
4. Generated plots are saved in `output/plots/*.png`

### Config Options
1. Available config options for the the 

## File Structure
### Functions are stored in /functions
#### Eady-Type Model Scripts
- `hwe_BPVyCalc.m`: Computes the meridional PV gradient at interior grid points for the Eady-type model, using a $\cos^4$-modulated mean zonal wind $\bar{u}(y, z)$. Essential for wave propagation studies.
- `hwe_PV2bndgrad.m`: Calculates PV gradients at the surface ($z=0$) and tropopause ($z=HH$) for the Eady-type model, incorporating vertical shear with $\cos^4$ meridional variation. Critical for boundary condition enforcement.
- `hwe_PV2intgrad.m`: Determines the interior PV gradient for the Eady-type model, including the meridional second derivative of the mean flow. Supports stability and dynamics analysis.
- `hwe_ubar.m`: Generates the mean zonal wind field $\bar{u}(y, z)$ for the Eady-type model, using linear shear with $\cos^4$ meridional modulation. Foundational for background flow setup.

#### Hoskins-West Model Scripts (Legacy)
- `hwm_BPVyCalc.m`: Computes the meridional PV gradient for the original Hoskins-West model, using sinusoidal modulation and vertical shear terms. Included for reference.
- `hwm_PV2bndgrad.m`: Calculates PV gradients at boundaries for the Hoskins-West model, incorporating complex shear terms. Retained for historical context.
- `hwm_PV2intgrad.m`: Computes the interior PV gradient for the Hoskins-West model, using planetary $\beta$ and vertical shear effects. Kept for comparison.
- `hwm_ubar.m`: Generates the mean zonal wind field for the Hoskins-West model with sinh-based modulation. Preserved for legacy support.

#### Utility Scripts
- `jk2l.m`: Converts 2D indices $(j, k)$ to a single linear index for grid operations, facilitating array manipulation.
- `jk2lw.m`: Transforms 2D indices $(j, k)$ to a linear index with weighting, used for specific grid indexing tasks.
- `l2jk.m`: Converts a linear index back to 2D indices $(j, k)$, supporting grid coordinate recovery.
- `lw2jk.m`: Converts a weighted linear index back to 2D indices $(j, k)$, aiding in weighted grid operations.
- `matricesBCD.m`: Constructs matrices $B$, $C$, and $D$ for PV inversion and advection in the eigenvalue problem, central to stability analysis.
- `row_F123.m`: Generates row vectors $F1$, $F2$, and $F3$ for matrix assembly, supporting the linearized PV equation.
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
### Data files are in /data



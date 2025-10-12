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
- `hwe_BPVyCalc.m` - Computes the meridional potential vorticity (PV) gradient for the Eady-type model at interior grid points, using the new $  \bar{u}(y, z)  $ formulation with $  \cos^4  $ modulation, essential for wave propagation studies.
- `hwe_PV2bndgrad.m`
`hwe_PV2intgrad.m`
`hwe_ubar.m`
`hwm_BPVyCalc.m`
`hwm_PV2bndgrad.m`
`hwm_PV2intgrad.m`
`hwm_ubar.m`
`jk2l.m`
`jk2lw.m`
`l2jk.m`
`lw2jk.m`
`matricesBCD.m`
`row_F123.m`
`stream2pv.m`
`stream2xPVadv.m`
`stream2yPVadv.m`
`w2ellipse.m`
`w2wfield.m`
`XV2field.m`
`XV2streamxtime.m`
`XV2XVx.m`
`XV2XVy.m`
`XV2XVz.m`
`XVx2field.m`
`XVy2field.m`
`XVz2field.m`

### Plots are stored in /plots
### Data files are in /data



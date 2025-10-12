# MET5510 Project
This directory contains the files and code for Grace and Matthew's MET5510 project under the guidance of Dr. Cai and Dr. Sun.

## Project Overview
The program simulates various atmospheric wave patterns using different modeling techniques for mid-latitude systems. Currently, simulations can be run for Rossby waves and the modified Hoskins-West model.

## Instructions
### Calculate Model Data
1. Run `row_main.m` to process the Rossby Wave Model calculations
2. Run `hwm_main.m` to process the Hoskins-West Modified calculations
3. Run `hwe_main.m` to process the Hoskins-West Eady-type calculations
4. Generated data files are saved in `output/data/*.mat`

### Create Initial Plots
1. Run `row_plot.m` to generate initial plots using the generated Rossby data
2. Run `hwm_plot.m` to generate initial plots using the generated Hoskins-West Modified data
3. Run `hwe_plot.m` to generate initial plots using the generated Hoskins-West Eady-type data
4. Additional Rossby plots are available using `rossby_vmotion.m`
5. Generated plots are saved in `output/plots/*.png`

### Config Options
1. Available config options for the the 

## File Structure
### Functions are stored in /functions
### Plots are stored in /plots
### Data files are in /data



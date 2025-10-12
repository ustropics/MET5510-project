# MET5510 Project

This directory contains the files and code for Grace and Matthew's MET5510 project under the guidance of Dr. Cai and Dr. Sun.

## Project Overview
The program simulates various atmospheric wave patterns using different modeling techniques for mid-latitude systems. Currently, simulations can be run for 

## Instructions
### Generate Data Files
1. Run `rossby_main.m` to generate an initial rossby_wave_#.mat file
2. Run `hoskin_main.m` to generate an initial hoskin_wave_#.mat file
3. Generated data files are saved in `data/*.mat`

### Create Initial Plots
1. Run `rossby_plot.m` to generate initial plots using the generated Rossby data
2. Run `hoskin_plot.m` to generate initial plots using the generated Hoskin data
3. Additional Rossby plots are available using `rossby_vmotion.m`
4. Generated plots are saved in `plots/*.png`

## File Structure
### Functions are stored in /functions
### Plots are stored in /plots
### Data files are in /data



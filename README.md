[![Build Status](https://travis-ci.com/cmhoove14/optimal_control.svg?branch=master)](https://travis-ci.com/cmhoove14/optimal_control)

## Team members:  
+ Chris Hoover, cmhoove14

## Project:  
Using optimal control theory to identify cost-optimal interventions for control and elimination of schistosomiasis in different settings. Building on a previous project for ESPM 288 investigating optimal control of a simple infectious disease system, work for which can be found in the **Archive** folder. 

### data folder  
Contains raw data used for this project such as epi data (that's non-sensitive) and large `.RData` files that are outputs of more intensive simulations  

### analysis folder  
Contains two `.Rmd`s and their outputs. `manuscript.Rmd` is text of the manuscruipt that's being generated related to this work. `notebook.Rmd` is a daily tracker of my thoughts and actions related to the work I do on this project

### `R` folder  
Contains scripts that are used in the analysis folder to generate figures, test hypotheses, test functions, generate new questions, etc.  

## Common files:  
+ `README.md`: This file
+ `gitignore`: Files ignored by Git commands

## Infrastructure for testing:  
+ `.travis.yml`: Configuration file for automatically running [continuous integration](https://travis-ci.com) checks to verify reproducibility of all `.Rmd` notebooks in the repo.  If all `.Rmd` notebooks can render successfully, the "Build Status" badge above will be green (`build success`), otherwise it will be red (`build failure`).  
+ `DESCRIPTION` : Metadata file for the repository, based on the R package standard. Lists any additional R packages/libraries needed for any of the `.Rmd` files to run.
+ `render_rmds.R` : R script that is run to execute the above described tests, rendering all `.Rmd` notebooks.

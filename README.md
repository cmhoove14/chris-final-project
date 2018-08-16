[![Build Status](https://travis-ci.com/cmhoove14/optimal_control.svg?branch=master)](https://travis-ci.com/cmhoove14/optimal_control)

## Team members:  
+ Chris Hoover, cmhoove14

## Project:  
Stochastic dynamic programming (SDP) to investigate optimal control of a simple infectious disease system. Work is in the *Project* folder and project proposal is in the *Proposal* folder.

### Common files:  
+ `README.md`: This file
+ `gitignore`: Files ignored by Git commands

### Infrastructure for testing:  
+ `.travis.yml`: Configuration file for automatically running [continuous integration](https://travis-ci.com) checks to verify reproducibility of all `.Rmd` notebooks in the repo.  If all `.Rmd` notebooks can render successfully, the "Build Status" badge above will be green (`build success`), otherwise it will be red (`build failure`).  
+ `DESCRIPTION` : Metadata file for the repository, based on the R package standard. Lists any additional R packages/libraries needed for any of the `.Rmd` files to run.
+ `render_rmds.R` : R script that is run to execute the above described tests, rendering all `.Rmd` notebooks.

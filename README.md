# PHSX815 Spring 2021 Project 4

## Computational Methods in Physics

This repository contains several types of programs:

- `PHypoSim.x` : generates and exports Poisson distributed data samples
                 with user defined output file, background rate mean and sigma,
                 varied signal rate, experiments per signal rate, and exposure [C++]
- `PHypoTest.x` : performs analysis on Poisson data sets to determine and
                  visualize the relationship between LLR hypothesis test significance
                  and signal rate [C++]

### Requirements

In order to compile (by typing `make`) and run the C++ examples, you
need the ROOT package installed (for visualization):
- [ROOT](https://root.cern/) (C++)

### Usage

All of the executables can be called from the
command line with the `-h` or `--help` flag, which will print the options

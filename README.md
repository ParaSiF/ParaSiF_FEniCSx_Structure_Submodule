# ParaSiF FEniCSx Solver

## Overview

ParaSiF FEniCSx (v0.9.0) solver

## Folder Structure

* src

It contains the source code on the FEniCSx solvers for ParaSiF coupled framework.

* single_benchmark

It contains the benchmark test case without MUI related coupling.

* coupled_benchmark

It contains the benchmark case on the FEniCSx linear elasticity solver, based on test case 10.2. Three-dimensional cantilever of Slone et al. 2003:

Slone, A. K., C. Bailey, and M. Cross. "Dynamic solid mechanics using finite volume methods." Applied mathematical modelling 27.2 (2003): 69-87.

There are two sub-folders:

1. dummyOF/
It contains a simple C++ script as a dummy fluid solver. It's function is to pass the node forces to the structure domain through MUI library.

2. structureDomain/
It contains the input files of FEniCSx solver, includes the main file (structureDomainRun.py), BCs, SubDomains and control parameters. A new subfolder "structureResults" will be generated to collect the results from the FEniCSx solver.

## Install

The dependencies of this ParaSiF FEniCSx are (correct version of) FEniCSx_v0.9.0 and MUI_v2.0.

* Step One: Install FEniCSx (v0.9.0).

Following FEniCSx homepage (https://fenicsproject.org/) for the installation procedure.

* Step Two: Install MUI Python wrapper by (it may take 5-10 min to compile, depends on the performance of the machine):

```
cd ParaSiF/coupling/MUI/wrappers/Python

make USE_RBF=1 mui4py_mod
```

After these two steps, the installation has been finished.

## Run the Benchmark case

* Step One: go to the coupled_benchmark (or single_benchmark) folder.


* Step Two: Run the case.

Go to the root of this case and execute the run script by

```
./Allrun
```
To clean-up the results from previous run, execute the clean script before run

```
./Allclean
```

* Step Three: Check results.

For the benchmark case, once the simulation finished, there will be a PNG file (result_compare.png) generated in the root of this benchmark case. Open it to check the results compared with published results from Slone et al. 2003.

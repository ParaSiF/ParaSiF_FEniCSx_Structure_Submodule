# ParaSiF FEniCSx Structure Solver

This repository contains the **FEniCSx Structure Solver** integrated with the ParaSiF Parallel Partitioned Simulation Framework.
It is maintained as a **submodule** of the main ParaSiF repository: [ParaSiF Main Repository](https://github.com/ParaSiF/ParaSiF).

---

## Overview

The FEniCSx Structure Solver in ParaSiF allows the simulation of *structural domains* in **multi-physics partitioned simulations**.
It is designed to interface with other solvers (e.g., fluid solvers) via the **[MUI coupling library](https://mxui.github.io/)**.

Key features:

- Linear and Hyper- elasticity solvers.
- Supports parallel execution for high-performance simulations.
- Modular design: can be replaced or updated independently of other solvers.
- Compatible with precompiled or user-installed FEniCSx versions.

---

## Compatible Codebase

This solver has been tested and is compatible with **[FEniCSx v0.9.0](https://github.com/FEniCS/dolfinx/releases/tag/v0.9.0.post1)**.

> Users are recommended to use this version to ensure full compatibility with ParaSiF FEniCSx solvers.

---

## Location in the Main ParaSiF Repository

`ParaSiF/src/structure/FEniCSx/`

---

## Repository Structure

```
ParaSiF/src/structure/FEniCSx/
├── doc/                  # Documentation folder
├── src/                  # ParaSiF-specific source code folder
│ └── structureFSISolver/ # FEniCSx solver source code
└── test/                 # test folder
  ├── single_benchmark/    # Benchmark test case without MUI coupling
  └── coupled_benchmark/   # Benchmark case with FEniCSx solver coupled with dummy solver via MUI
```

---

## Installation

**Note:** This solver is a submodule of ParaSiF. Follow the main ParaSiF repository instructions to initialise submodules and install global dependencies.

### Steps

1. **Obtain and install the codebase**
   - Initialise the submodule from the main ParaSiF repository (or clone this repository).
   - Ensure the correct version of FEniCSx is installed by following instructions from the [FEniCSx homepage](https://fenicsproject.org/) (The Spack approach is recommended).

2. **Install MUI Python wrapper**

```bash
cd ParaSiF/coupling_lib/MUI/wrappers/Python
```
- Install MUI Python wrapper by following instructions in the [MUI reposirory](https://github.com/MxUI/MUI)

> There is no need to compile or install the ParaSiF FEniCSx src/structureFSISolver. These Python-based FEniCSx codes are compiled automatically at runtime when executing the solver scripts.

## Running Tests and Example Cases

Benchmark cases are located in the test/ folder:

### Steps

1. Navigate to the desired benchmark folder:

```bash
cd test/XXX
```

2. Run the simulation:

```bash
./Allrun.sh
```

3. (Optional) Clean up previous results before rerunning:

```bash
./Allclean.sh
```
4. Check results:

During runtime, a `structureResults/` folder will be automatically generated in `structureDomain/` to store the solver output (displacements, stresses, checkpoint data, etc.).

For the benchmark case, a `result_compare.png` file will be generated in the benchmark folder. This allows comparison of simulation results with published results from Slone et al. 2003.

For integrated example cases with other solvers, see the example/ folder in the main ParaSiF repository.

## Contributing

ParaSiF, including this submodule, is an **open-source project**, and contributions from the community are warmly welcomed.

There are many ways you can help improve this submodule, including:

- Adding new features, libs or solvers
- Improving documentation, tests and examples
- Fixing bugs or refining existing functionality
- Sharing feedback and suggestions for enhancements

Your contributions, whether large or small, are highly valued and help make ParaSiF a stronger resource for the research community.

For detailed guidance on contributing, please see the [CONTRIBUTING.md](https://github.com/ParaSiF/ParaSiF/blob/main/CONTRIBUTING.md) in the main ParaSiF repository.

## License

Copyright (C) 2021–2025 The ParaSiF Development Team.  
Licensed under the **GNU General Public License v3 (GPL-3.0)**.

## Contact

For questions or contributions, please contact the ParaSiF Development Team

## References
- Slone, A. K., Bailey, C., & Cross, M. (2003). Dynamic solid mechanics using finite volume methods. Applied Mathematical Modelling, 27(2), 69-87.

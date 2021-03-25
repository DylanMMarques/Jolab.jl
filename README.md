# Jolab
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://dylanmmarques.github.io/Jolab.jl/)
[![dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dylanmmarques.github.io/Jolab.jl/dev/)
[![Build status](https://github.com/DylanMMarques/Jolab.jl/workflows/CI/badge.svg)](https://github.com/DylanMMarques/Jolab.jl/actions)
[![Codecov](https://codecov.io/gh/DylanMMarques/Jolab.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/DylanMMarques/Jolab.jl)
[![Coveralls](https://coveralls.io/repos/github/DylanMMarques/Jolab.jl/badge.svg?branch=master)](https://coveralls.io/github/DylanMMarques/Jolab.jl?branch=master)

Jolab is a free and open-source Julia package to simulate light propagation in optical systems. The package implements rigorous, full-wave optical models to simulate light propagation on setups. The package is built to be easily usable by any user independent of its background in theoretical optics.

The package was designed to allow an user to simulate optical setups based on individual optical components, i.e., the package provides functions that calculate how light propagates on an optical component and the user must propagate the light by the many optical components. The simulation of an optical system is made by simulating how light propagates in each individual optical component.

The package is currently at a very early stage in its development. Due to that, it is very likely the code to have bugs and inconsistencies. Therefore, we advise any user to be careful in its utilization. In the near future, we will focus the development in creating a stable version of the package with scientific rigor.

In the future, we want to develop the package to allow users to simulate a variety of optical systems as microscopes, Optical Coherence Tomography setups (OCT), interferometers, ellipsometers, etc. We want the package to be a tool to help scientists and engineers to design and understand their optical setups.

[![Video introduction of Jolab](https://raw.githubusercontent.com/DylanMMarques/Jolab.jl/master/docs/src/assets/presentation_spie.png)](https://youtu.be/P-S6lMDnE7Y)

---

## Instalation
Jolab is a [Julia](https://julialang.org/) package. To use it, an user only needs to install [Julia](https://julialang.org/) and run:
```julia
] add Jolab
```

---

## Supporting and Citing
Jolab is a community project to create an environment of optical models to simulate systems. If you want to support the project please star the repository, as it might help to secure funding in the future. If you use Jolab in your research please cite:

[Dylan Marques, James A. Guggenheim, and Peter R. T. Munro "Jolab a free and open-source software to simulate light propagation in optical systems", Proc. SPIE 11649, Three-Dimensional and Multidimensional Microscopy: Image Acquisition and Processing XXVIII, 1164914 (5 March 2021); https://doi.org/10.1117/12.2577631](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/11649/2577631/Jolab-a-free-and-open-source-software-to-simulate-light/10.1117/12.2577631.short?SSO=1)

Cite as well the models that you are using. You have more information about the citation of each model in the documentation.

---

The package is developed at the Department of Medical Physics and Biomedical Engineering of UCL (University College London) by Dylan M. Marques, James A. Guggenheim, and Peter R. T. Munro.

Please feel free to get in touch to report bugs, comments and feedback (dylan.marques.17@ucl.ac.uk)

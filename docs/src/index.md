# Jolab

[![Build Status](https://travis-ci.com/DylanMMarques/Jolab.jl.svg?branch=master)](https://travis-ci.com/DylanMMarques/Jolab.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/DylanMMarques/Jolab.jl?svg=true)](https://ci.appveyor.com/project/DylanMMarques/Jolab-jl)
[![Codecov](https://codecov.io/gh/DylanMMarques/Jolab.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/DylanMMarques/Jolab.jl)

Jolab is a free, open-source Julia package to simulate light propagation in optical setups. The package implements rigorous, full-wave optical theories to simulate light propagation on setups. The package is built to be easily usable by any user independent of its background in theoretical optics.

The package was designed to allow an user to simulate optical setups based on individual optical components, i.e., the package provides functions that calculate how light propagates on an optical component and the user must propagate the light by the many optical components. The simulation of an optical setup is made by simulating how light propagates in each individual optical component.
A few examples are provided in the folder example to show the principle of the toolbox.

The package is currently at a very early stage in its development. Due to that, it is very likely the code to have bugs and inconsistencies. Therefore, we advise any user to be careful in its utilization. In the near future, we will focus the development in creating a stable version of the package with scientific rigor.

In the future, we want to develop the package to allow users to simulate a variety of optical setups as microscopes, Optical Coherence Tomography setups (OCT), interferometers, ellipsometers, etc. We want the package to be a tool to help scientists and engineers to design and understand their optical setups.

## Instalation
To install Jolab, a user need to run:
```julia
] add https://github.com/DylanMMarques/Jolab.jl
```

The package is developed at the Department of Medical Physics and Biomedical Engineering of UCL (University College London) by Dylan M. Marques, James A. Guggenheim, and Peter R. T. Munro.

Please feel free to get in touch to report bugs, comments and feedback (dylan.marques.17@ucl.ac.uk)

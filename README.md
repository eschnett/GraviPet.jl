![GraviPet logo](figures/GraviPet-light-background.jpg)

# GraviPet.jl

GraviPet is the **G**eneral **R**elativistic **A**strophysics
**V**isualization, **I**nitialization, and **P**ostprocessing
**E**fficient **T**oolkit.

* [![Documentation](https://img.shields.io/badge/Docs-Dev-blue.svg)](https://eschnett.github.io/GraviPet.jl/dev/)
* [![GitHub
  CI](https://github.com/eschnett/GraviPet.jl/workflows/CI/badge.svg)](https://github.com/eschnett/GraviPet.jl/actions)
* [![codecov](https://codecov.io/gh/eschnett/GraviPet.jl/graph/badge.svg?token=VGMG5U8M41)](https://codecov.io/gh/eschnett/GraviPet.jl)

GraviPet provides ways to represent functions, for example solutions
or initial conditions for PDEs (partial differential equations). Often
such functions are discretized, i.e. they are represented in terms
of a finite number of basis functions. Common choices are sampling
function at grid points (finite differencing) or averaging functions
over grid cells (finite volumes). Many other choices exist.

The main design idea behind GraviPet is to offer *composable
abstractions*. That is, instead of multi-threaded finite differencing
grid function, there exist basic (serial) grid functions, as well as
adapters that render any other kind of grid function multi-threaded.

## Acknowledgements

The GraviPet logo was created by [Grabriela
Secara](https://perimeterinstitute.ca/people/gabriela-secara) at the
[Perimeter Institute for Theoretical
Physics](https://perimeterinstitute.ca/).

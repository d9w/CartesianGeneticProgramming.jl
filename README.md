# CartesianGeneticProgramming.jl

[![Build Status](https://travis-ci.org/d9w/CartesianGeneticProgramming.jl.svg?branch=master)](https://travis-ci.org/d9w/CartesianGeneticProgramming.jl) [![Coverage Status](https://coveralls.io/repos/d9w/CartesianGeneticProgramming.jl/badge.svg?branch=master)](https://coveralls.io/r/d9w/CartesianGeneticProgramming.jl?branch=master) [![codecov](https://codecov.io/gh/d9w/CartesianGeneticProgramming.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/d9w/CartesianGeneticProgramming.jl)

[Cartesian Genetic Programming](http://www.cartesiangp.com/) in
[Julia](http://julialang.org/). Based on the
[Cambrian.jl](https://github.com/d9w/Cambrian.jl) framework.

## Installation

CartesianGeneticProgramming.jl can be installed through the Julia package manager

```julia
pkg> add CartesianGeneticProgramming.jl
```

## Examples

Examples of usage are provided in `scripts`. To run any of these examples, do

```bash
$ julia scripts/gym.jl
```

This will create and populate `gens` and `logs` folders, which store CGP
individuals and evolutionary logs, respectively.

Note that these examples have only been tested on Linux and have additional
dependencies which must be installed, like `PyCall` or
`ArcadeLearningEnvironment`.

## Tests

`CGP.jl` comes with tests which can offer examples of detailed usage of the
different genetic operators and CGP extensions. To run all tests, use

```julia
pkg> test CartesianGeneticProgramming
```

## Configuration

Configuration is handled by `src/config.jl`. Configurable options are read in
through YAML files (an example being `cfg/test.yaml`). All possible node
functions are defined in `src/functions.jl` and are populated in the
configuration files.

## Mixed-Type CGP

This repository is now intended only for floating point inputs and therefore no
longer implements Mixed-Type CGP. The commit reference for the results in [1] is
[096e843](https://github.com/d9w/CartesianGeneticProgramming.jl/commit/096e84316c0c265e17ea75de7b5fcaa6427fa5dc),
however this requires Julia 0.6 and may no longer function. MTCGP is currently
being rewritten at https://github.com/d9w/MTCGP.jl. This is done to separate
heavy image processing package dependencies which are not necessary for base
CGP.

[1] Wilson, Dennis G., et al. "Evolving simple programs for playing Atari
games." Proceedings of the Genetic and Evolutionary Computation Conference. ACM,
2018.

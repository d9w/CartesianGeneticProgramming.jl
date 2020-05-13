# CartesianGeneticProgramming.jl

[![Build Status](https://travis-ci.org/d9w/CartesianGeneticProgramming.jl.svg?branch=master)](https://travis-ci.org/d9w/CartesianGeneticProgramming.jl) [![Coverage Status](https://coveralls.io/repos/d9w/CartesianGeneticProgramming.jl/badge.svg?branch=master)](https://coveralls.io/r/d9w/CartesianGeneticProgramming.jl?branch=master) [![codecov](https://codecov.io/gh/d9w/CartesianGeneticProgramming.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/d9w/CartesianGeneticProgramming.jl)

[Cartesian Genetic Programming](http://www.cartesiangp.co.uk/) for
[julia](http://julialang.org/).

Development is still ongoing, so be warned that major changes might happen.
Stable resources in other languages can be found on
the [CGP site](https://www.cartesiangp.com/).

## Installation

CartesianGeneticProgramming.jl can be installed through the Julia package manager

```julia
pkg> add https://github.com/d9w/CartesianGeneticProgramming.jl
```

This package depends on [Cambrian.jl](https://github.com/d9w/Cambrian.jl) which
is not yet in the General repository and therefore must also be installed from
github.

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

```julilang
pkg> test CartesianGeneticProgramming
```

## Configuration

Configuration is handled by `src/config.jl`. Configurable options are read in
through YAML files (an example being `cfg/test.yaml`). All possible node
functions are defined in `src/functions.jl` and are populated in the
configuration files.

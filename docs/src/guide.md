# Package Guide

CatmapInterface implements the functionality of a subset of the Python package [CatMAP](https://catmap.readthedocs.io) that is used in the Python package [CatINT](https://catint.readthedocs.io).
Whereas [CatINT](https://catint.readthedocs.io) uses an iterative approach for combining the Poisson-Nernst-Planck transport model with a microkinetic model of the surface reactions at the electrode, the output of CatmapInterface can be directly plugged into the functionality of [LiquidElectrolytes](https://j-fu.github.io/LiquidElectrolytes.jl) to solve the coupled system.

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add https://github.com/smaasz/CatmapInterface.jl
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("https://github.com/smaasz/CatmapInterface.jl")
```

Finally, you need to include it to the current module
```julia
julia> using CatmapInterface
```

## Parsing a CatMAP Input File

!!! note
    CatmapInterface only supports a subset of the complete CatMAP functionality.
    See [`parse_catmap_input`](@ref) for details.

Assuming that the microkinetic model is described in a CatMAP input file stored locally at `input_file_path` (see e.g. [CO2-Reduction on Ag](https://github.com/sringe/CatINT/blob/master/examples/02_CO2R_Au_CatMAP/catmap_CO2R_template.mkm)), the input can be parsed by
```julia
julia> catmap_params = parse_catmap_input(input_file_path)
```
which returns a [`CatmapParams`](@ref) struct that stores the needed information.

## Creating the Microkinetic Model

CatmapInterface assembles the microkinetic model as a [ReactionSystem](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/#Catalyst.ReactionSystem) implemented in [Catalyst.jl](https://docs.sciml.ai/Catalyst/):
```julia
julia> microkinetic_rn = create_reaction_network(catmap_params)
```
It includes all specified species but neglects the solvent, water, and the ficitious species of electrons.
```julia
julia> using Pkg; Pkg.add("Catalyst"); using Catalyst
julia> species(microkinetic_rn)
```
The reaction rates are computed symbolically in the parameters of the model
```julia
julia> parameters(microkinetic_rn)
```

## Use the Microkinetic Model in LiquidElectrolytes

To use the rate equations from the microkinetic model in the Poisson-Nernst-Planck transport model implemented in LiquidElectrolytes, an in-place function definition (see [Interface of DiffEq.jl](https://docs.sciml.ai/DiffEqDocs/stable/basics/problem/#Problem-Interface)) of the (micro-)kinetic model ODE can be generated
```julia
julia> f_microkinetics! = CatmapInterface.generate_function(microkinetic_rn, dvs, ps)
```
where `dvs` and `ps` define the order the species and parameters, respectively. See [`CatmapInterface.generate_function`](@ref) for details.
# CatmapInterface

| **Documentation**                 | **Build Status**        |
|:---------------------------------:|:-----------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://smaasz.github.io/CatmapInterface.jl/dev/) | [![](https://github.com/smaasz/CatmapInterface.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/smaasz/CatmapInterface.jl/actions/workflows/CI.yml?query=branch%3Amain) |

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> registry add https://github.com/j-fu/PackageNursery
pkg> add CatmapInterface
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.Registry.add("https://github.com/j-fu/PackageNursery"); Pkg.add("CatmapInterface")
```

## Guide Outline

```@contents
Pages = [
    "guide.md"
]
```

## Library Outline

```@contents
Pages = [
    "public.md",
    "internal.md"
]
Depth = 1
```

## Index

```@index
Pages = ["public.md"]
```
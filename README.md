# Schrodinger1d

This is a Julia package for solving 1D schrodinger equation.

## Installation
This is an unregistered package, and can be installed in the following way:

```julia
julia> pkg"add git@github.com:Genshin1024/schrodinger1d.git"
```

## Usage

```julia
julia> potential(x) = 0.5*x^2


julia> energy, psi = schrodinger1d.solver(potential)

```


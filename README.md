# Schrodinger1d

This is a toy Julia package for solving 1D schrodinger equation.

## Installation
This is an unregistered package, and can be installed in the following way:

```julia
julia> pkg"add git@github.com:Genshin1024/schrodinger1d.git"
```

## Usage

```julia
julia> using Schroinger1d, Printf,Plots

julia> potential(x) = 0.5*x^2

julia> energy, psi = schrodinger1d.solver(potential);

julia> title_str =  @sprintf "ground energy = %6.2f" energy
"ground energy =   0.50"

julia> plot(psi, title=title_str)

```

![psi with energy=0.5](https://picutre999.oss-cn-hangzhou.aliyuncs.com/img/E0.5psi.jpg)

using Printf
using Plots
using schrodinger1d

# harmonic oscillator
potential(x) = 0.5*x^2

energy, psi = schrodinger1d.solver(potential)

title_str =  @sprintf "ground energy = %6.2f" energy
plot(psi, title=title_str)
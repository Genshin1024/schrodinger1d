using Schrodinger1d
using Printf
using Plots
using Test

# harmonic oscillator 
@testset "Schrodinger1d.jl" begin
    include("poltpsi.jl")
    @test plotpsi()
end

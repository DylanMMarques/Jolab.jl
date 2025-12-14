# import Pkg
# Pkg.activate(@__DIR__)
using Jolab, Test, Enzyme, StaticArrays, FiniteDiff
import FiniteDiff: finite_difference_derivative
import Jolab: Forward, Backward

include(joinpath(@__DIR__, "auxialiary_functions.jl"))
@testset "All" begin
    @testset "General beam" begin include("beam.jl") end
    @testset "Reference frames" begin include("reference_frames.jl") end
    @testset "Lens" begin include("lens.jl") end
    @testset "Mirror" begin include("mirror.jl") end
    @testset "Dielectric stack" begin include("dielectric_stack.jl") end
    @testset "Step index waveguides" begin include("circularindexfibre.jl") end
    @testset "Single mode fibre" begin include("single_mode_fibre.jl") end
    @testset "Fourier Transform" begin include("fourier.jl") end
end

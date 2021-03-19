using Jolab, Test

@testset "Fourier transform" begin include("fouriertransform.jl") end
@testset "Multilayer Structure" begin include("multilayerstructure.jl") end
@testset "Rough Multilayer Structure" begin include("roughinterface.jl") end
@testset "Axicon" begin include("axicon.jl") end
@testset "Mirror" begin include("mirror.jl") end
@testset "Lens" begin include("lens.jl") end
@testset "Recursive Algorithms" begin include("lightinteraction_recursive.jl") end

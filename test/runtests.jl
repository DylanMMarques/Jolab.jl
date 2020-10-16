using Jolab, Test

@testset "Jolab.jl" begin
    @test include("fouriertransform.jl")
    @test include("multilayerstructure.jl")
    @test include("roughinterface.jl")
    @test include("axicon.jl")
end

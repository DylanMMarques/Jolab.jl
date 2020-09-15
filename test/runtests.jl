using Jolab, Test

@testset "Jolab.jl" begin
    @test include("multilayerstructure.jl")
    @test include("roughinterface.jl")
    @test include("fouriertransform.jl")
end

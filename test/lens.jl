using Jolab, Test

focal_len = 1e-3
lens = Lens(focal_len, .5, (Medium(1), Medium(1.0)), ReferenceFrame((0,0, focal_len), (0,0,0)))

@test_throws ArgumentError Lens(1E-3, .5, (Medium(1 + im), Medium(1.0)), ReferenceFrame((0,0,0), (0,0,0)))

nsx = range(-1, 1, length=10)
λ = 1500E-9
grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, λ)
field = Jolab.AngularSpectrum(Jolab.Forward, grid_ang, rand(ComplexF64, 10, 10), Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, λ)
field = Jolab.AngularSpectrum(Jolab.Forward, grid_ang, rand(ComplexF64, 10, 10), Medium(1), ReferenceFrame((0,0,1), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, λ)
field = Jolab.AngularSpectrum(Jolab.Forward, grid_ang, rand(ComplexF64, 10, 10), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
(rfield, tfield) = light_interaction(lens, field)
@test_opt light_interaction(lens, field)
@test tfield isa Jolab.MeshedSpatialBeam
@test rfield isa Jolab.MeshedAngularSpectrum


x = range(-1E-3, 1E-3, length = 10)
λ = 1500E-9
grid = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ)
field = Jolab.SpatialBeam(Jolab.Forward, grid, rand(ComplexF64, 10, 10), Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)
grid = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ)
field = Jolab.SpatialBeam(Jolab.Forward, grid, rand(ComplexF64, 10, 10), Medium(1), ReferenceFrame((0,0,1), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

grid = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ)
field = Jolab.SpatialBeam(Jolab.Backward, grid, rand(ComplexF64, 10, 10), Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)
grid = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ)
field = Jolab.SpatialBeam(Jolab.Backward, grid, rand(ComplexF64, 10, 10), Medium(1), ReferenceFrame((0,0,1), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

grid = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ)
field = Jolab.SpatialBeam(Jolab.Backward, grid, rand(ComplexF64, 10, 10), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

grid = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ)
field = Jolab.SpatialBeam(Jolab.Forward, grid, rand(ComplexF64, 10, 10), Medium(1), ReferenceFrame((0,0,2focal_len), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

x = LinRange(-50E-6, 50E-6, 100)
λ = 1500E-9
grid_gauss = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ)
field = Jolab.SpatialBeam(Jolab.Forward, grid_gauss, Jolab.Gaussian(10E-6), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test_opt light_interaction(lens, field)
(bfield, ffield) = light_interaction(lens, field)
@test bfield isa Jolab.MeshedSpatialBeam
@test ffield isa Jolab.MeshedAngularSpectrum
@test iszero(intensity(bfield))
@test isapprox(intensity(ffield), intensity(field); rtol = 1E-5)

x = LinRange(-100E-6, 100E-6, 100)
grid_gauss = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ)
field = Jolab.SpatialBeam(Jolab.Backward, grid_gauss, Jolab.Gaussian(10E-6), Medium(1), ReferenceFrame((0,0,2focal_len), (0,0,0)))
(bfield, ffield) = light_interaction(lens, field)
@test bfield isa Jolab.MeshedAngularSpectrum
@test ffield isa Jolab.MeshedSpatialBeam
@test iszero(intensity(ffield))
@test isapprox(intensity(bfield), intensity(field); rtol = 1E-5)

sca = ScatteringMatrix(lens, field)
@test all(light_interaction(sca, field) .≈ light_interaction(lens, field))

x = LinRange(-.5, .5, 1000)
grid_ang_gauss = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, x, x, λ)
field = Jolab.AngularSpectrum(Jolab.Backward, grid_ang_gauss, Jolab.Gaussian(10E-6), Medium(1), ReferenceFrame((0,0,2focal_len), (0,0,0)))
(bfield, ffield) = light_interaction(lens, field)
@test bfield isa Jolab.MeshedSpatialBeam
@test ffield isa Jolab.MeshedAngularSpectrum
@test iszero(intensity(ffield))
@test isapprox(intensity(bfield), intensity(field); rtol = 1E-2) # Not sure if correct

@test_opt ScatteringMatrix(lens, field)
sca = ScatteringMatrix(lens, field)
@test all(light_interaction(sca, field) .≈ light_interaction(lens, field))

stack = DielectricStack(Medium.([1.0, 1.5, 1.0]), [500E-9], ReferenceFrame((0,0,0.002), (0,0,0)))
x = range(-1, 1, length = 512) * lens.focal_length
e = zeros(ComplexF64, length(x), length(x), 1, 3)
e[:,:,:,1] .= 1
begin
    lens = Lens(focal_len, 1.0, (Medium(1), Medium(1.0)), ReferenceFrame((0,0, focal_len), (0,0,0)))
    grid_xyz = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, 1550E-9)
    beam = Jolab.MeshedBeam{Float64, Jolab.Forward, Jolab.X_Y_λ, Jolab.PolarizationXYZ}(grid_xyz, e, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
    
    @test_opt light_interaction(lens, beam)
    (beam_r, beam_f) = light_interaction(lens, beam)
    
    @test_opt light_interaction(stack, beam_f)
    (_, beam_end) = light_interaction(stack, beam_f)

    @test_opt Jolab.change_polarization_basis(beam_end, Jolab.PolarizationXYZ)
    beam_xyz = Jolab.change_polarization_basis(beam_end, Jolab.PolarizationXYZ)
    space_xyz = light_interaction(Jolab.FourierFFT(), beam_xyz)[2]
    sum(beam_xyz.e, dims = (1,2))
end

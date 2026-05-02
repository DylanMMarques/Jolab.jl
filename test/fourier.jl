using Jolab, Test

# test_opt fails because of FFT use

fft_op = Jolab.FourierFFT()

## Test Fourier Transform against analytical Gaussian beam transform
x = range(-200E-6, 200E-6, length = 1024)
λ = 1550E-9
grid_spatial = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ)
beam_space_analytical = Jolab.SpatialBeam(Float64, Forward, grid_spatial, Jolab.Gaussian(50E-6), Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))
(field_b, ang_spe) = light_interaction(fft_op, beam_space_analytical)

(nsx, nsy, λ_ang) = Jolab.get_ranges(ang_spe.mesh)
grid_angular = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsy, λ_ang)
ang_spe_analitical = Jolab.AngularSpectrum(Float64, Forward, grid_angular, Jolab.Gaussian(50E-6), Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))

@test ang_spe ≈ ang_spe_analitical
@test field_b isa Jolab.MeshedSpatialBeam
@test iszero(intensity(field_b))

(field_b, beam_space) = light_interaction(fft_op, ang_spe_analitical)

(x, y, λ_beam) = Jolab.get_ranges(beam_space.mesh)
grid_spatial_beam = Jolab.CartesianGrid(Jolab.X_Y_λ, x, y, λ_beam)
beam_space_analytical = Jolab.SpatialBeam(Float64, Forward, grid_spatial_beam, Jolab.Gaussian(50E-6), Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))

@test beam_space ≈ beam_space_analytical
@test field_b isa Jolab.MeshedAngularSpectrum
@test iszero(intensity(field_b))


## Test shift property of Fourier Transform
(x, y, λ_shift) = Jolab.get_ranges(beam_space.mesh)
grid_spatial_shift = Jolab.CartesianGrid(Jolab.X_Y_λ, x .+ x[512 + 124], y .+ y[512 + 124], λ_shift)
space_t = Jolab.SpatialBeam(Float64, Forward, grid_spatial_shift, Jolab.Gaussian(50E-6), Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))

(field_b, ang_spe) = light_interaction(fft_op, space_t)
(_, beam_space_trans) = light_interaction(fft_op, ang_spe)

@test findmax(abs, beam_space_trans.e)[2] == (findmax(abs, space_t.e)[2] + CartesianIndex(123, 123, 0, 0))

e = rand(ComplexF64, length(x), length(x), 1, 3)
grid_spatial_e = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, 1550E-9)
beam = Jolab.MeshedBeam{Float64, Forward, Jolab.X_Y_λ, Jolab.PolarizationXYZ}(grid_spatial_e, e, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
(field_b, ang_spe) = light_interaction(fft_op, beam)

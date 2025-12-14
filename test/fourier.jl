using Jolab, Test

fft_op = Jolab.FourierFFT()

## Test Fourier Transform against analytical Gaussian beam transform
x = range(-200E-6, 200E-6, length = 1024)
beam_space_analytical = MonochromaticSpatialBeam_gaussian(Float64, Forward, x, x, 50E-6, 1550E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))
(field_b, ang_spe) = light_interaction(fft_op, beam_space_analytical)

(nsx, nsy, λ) = Jolab.get_ranges(ang_spe.mesh)
ang_spe_analitical = MonochromaticAngularSpectrum_gaussian(Float64, Forward, nsx, nsy, 50E-6, 1550E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))

@test ang_spe ≈ ang_spe_analitical
@test field_b isa Jolab.MeshedSpatialBeam
@test iszero(intensity(field_b))

(field_b, beam_space) = light_interaction(fft_op, ang_spe_analitical)

(x, y, λ) = Jolab.get_ranges(beam_space.mesh)
beam_space_analytical = MonochromaticSpatialBeam_gaussian(Float64, Forward, x, y, 50E-6, 1550E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))

@test beam_space ≈ beam_space_analytical
@test field_b isa Jolab.MeshedAngularSpectrum
@test iszero(intensity(field_b))


## Test shift property of Fourier Transform
(x, y, λ) = Jolab.get_ranges(beam_space.mesh)
space_t = MonochromaticSpatialBeam_gaussian(Float64, Forward, x .+ x[512 + 124], y .+ y[512 + 124], 50E-6, 1550E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))

(field_b, ang_spe) = light_interaction(fft_op, space_t)
(_, beam_space_trans) = light_interaction(fft_op, ang_spe)

@test findmax(abs, beam_space_trans.e)[2] == (findmax(abs, space_t.e)[2] + CartesianIndex(123, 123, 0))

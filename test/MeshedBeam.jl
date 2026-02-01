using Jolab, Test

ns = range(-.1, 0.1, length = 100)
angspe = MonochromaticAngularSpectrum_gaussian(Float64, Forward, ns, ns, 50E-6, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(angspe) ≈ 1

x = range(-100E-6, 100E-6, length = 100)
beam = MonochromaticSpatialBeam_gaussian(Float64, Forward, x, x, 10E-6, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam) ≈ 1

x = range(-100E-6, 100E-6, length = 100)
beam = MonochromaticSpatialBeam_gaussian(Forward, x, x, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam) ≈ 1

nsx = range(-0.5, 0.5, length = 100)
beam = MonochromaticAngularSpectrum_gaussian(Forward, nsx, nsx, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam) ≈ 1

nsr = range(0, 0.5, length = 10000)
beam = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam) ≈ 1 rtol = 1E-3

r = range(0, 40E-6, length = 1000)
beam = MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam) ≈ 1 rtol = 1E-2

nsr = range(0, 0.5, length = 100)
nsy = range(0, stop = eps(), length = 2)
beam = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
beam_cart = MonochromaticAngularSpectrum_gaussian(Forward, nsr, nsy, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test beam_cart.e[:,1] ≈ beam.e

r = range(0, 40E-6, length = 100)
y = range(0, stop = eps(), length = 2)
beam = MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
beam_cart = MonochromaticSpatialBeam_gaussian(Forward, r, y, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test beam_cart.e[:,1] ≈ beam.e

view_beam = @view beam_cart[:,:,1:1]
@test view_beam.e ≈ beam_cart.e
@test view_beam.e ≈ beam_cart.e
@test view_beam.mesh ≈ beam_cart.mesh


x = range(-100E-6, 100E-6, length = 100)
e = zeros(ComplexF64, length(x), length(x), 3)
e[:,:,1] .= 1
beam = Jolab.MonochromaticSpatialBeamVectorial(Forward, x, x, e, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test size(beam.e) == (100, 100, 3)

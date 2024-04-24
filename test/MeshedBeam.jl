using Jolab

ns = range(-.1, 0.1, length = 100)
angspe = MonochromaticAngularSpectrum(Float64, Forward, ns, ns, rand(length(ns), length(ns)), 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
intensity(angspe)

x = range(-100E-6, 100E-6, length = 100)
beam = MonochromaticSpatialBeam(Float64, Forward, x, x, rand(length(x), length(x)), 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@time intensity(beam)

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
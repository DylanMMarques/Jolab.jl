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

x = range(-0.5, 0.5, length = 100)
beam = MonochromaticAngularSpectrum_gaussian(Forward, x, x, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam) ≈ 1
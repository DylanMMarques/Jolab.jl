using Jolab, Test

aperture = Jolab.SpatialLightModulator_aperture(Float64, 1E-1, Medium.((1,1)), ReferenceFrame((0,0,0), (0,0,0)))
disk = Jolab.SpatialLightModulator_diskaperture(Float64, 0.0, 1E-1, Medium.((1,1)), ReferenceFrame((0,0,0), (0,0,0)))

angspe = MonochromaticAngularSpectrum_gaussian(Float64, Forward, range(-0.1, 0.1, length=64), range(-0.1, 0.1, length=64), 10E-2, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test_throws MethodError light_interaction(disk, angspe)
beam = MonochromaticSpatialBeam_gaussian(Float64, Forward, range(-0.1, 0.1, length=64), range(-0.1, 0.1, length=64), 10E-2, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))

(_, t1) = light_interaction(aperture, beam)
(_, t2) = light_interaction(disk, beam)
@test t1 ≈ t2

@test_opt light_interaction(disk, beam)
@test_opt light_interaction(aperture, beam)


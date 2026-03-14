using Jolab, Test

mfd = 10E-6
frames = (ReferenceFrame((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), ReferenceFrame((0,0,1E-2), (0,0,0.0)))
media = Medium.((1,1))
fibre = SingleModeFibre(mfd, 1E-2, media, frames)
Jolab.findmodes!(fibre, 1500E-9)

x = LinRange(-50E-6, 50E-6, 1000)
spa =  MonochromaticSpatialBeam_gaussian(Forward, x, x, mfd, 1500E-9, Medium(1.0), frames[1])
(fieldr, fieldt) = light_interaction(fibre, spa)
@test intensity(fieldt) ≈ 1 rtol = 1E-3
field = MonochromaticSpatialBeam(Float64, Forward, fibre, x, x, 1500e-9)
@test intensity(field) ≈ 1 rtol = 1E-3

r = LinRange(0, 50E-6, 100000)
spa =  MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, mfd, 1500E-9, Medium(1.0), frames[1])
(fieldr, fieldt) = light_interaction(fibre, spa)
@test iszero(intensity(fieldr)) 
@test intensity(fieldt) ≈ 1 rtol = 1E-3
field = MonochromaticSpatialBeamRadialSymmetric(Float64, Forward, fibre, r, 1500e-9)
@test intensity(field) ≈ 1 rtol = 1E-3

nsx = LinRange(-0.5, 0.5, 1000)
field = MonochromaticAngularSpectrum_gaussian(Forward, nsx, nsx, mfd, 1500E-9, Medium(1.0), frames[1])
(fieldr, fieldt) = light_interaction(fibre, field)
@test iszero(intensity(fieldr))
@test intensity(fieldt) ≈ 1 rtol = 1E-3
field = MonochromaticAngularSpectrum(Float64, Forward, fibre, nsx, nsx, 1500e-9)
@test_opt MonochromaticAngularSpectrum(Float64, Forward, fibre, nsx, nsx, 1500e-9)
@test intensity(field) ≈ 1 rtol = 1E-3

nsr = LinRange(0, 0.5, 100000)
field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, mfd, 1500E-9, Medium(1.0), frames[1])
@test_opt light_interaction(fibre, field)
(fieldr, fieldt) = light_interaction(fibre, field)
@time light_interaction(fibre, field)
@test intensity(fieldt) ≈ 1 rtol = 1E-3
field = MonochromaticAngularSpectrumRadialSymmetric(Float64, Forward, fibre, nsr, 1500e-9)
@test intensity(field) ≈ 1 rtol = 1E-3

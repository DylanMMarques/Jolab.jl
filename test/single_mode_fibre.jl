using Jolab, Test

mfd = 10E-6
frames = (ReferenceFrame((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), ReferenceFrame((0,0,1E-2), (0,0,0.0)))
media = Medium.((1,1))
fibre = SingleModeFibre(mfd, 1E-2, media, frames)
Jolab.findmodes!(fibre, 1500E-9)

x = LinRange(-50E-6, 50E-6, 1000)
spa =  MonochromaticSpatialBeam_gaussian(Forward, x, x, mfd, 1500E-9, Medium(1.0), frames[1])
(a, b) = light_interaction(fibre, spa)
@test intensity(b) ≈ 1 rtol = 1E-3
field = MonochromaticSpatialBeam(Float64, Forward, fibre, x, x, 1500e-9)
@test intensity(field) ≈ 1 rtol = 1E-3

r = LinRange(0, 50E-6, 100000)#
spa =  MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, mfd, 1500E-9, Medium(1.0), frames[1])
(a, b) = light_interaction(fibre, spa)
@test intensity(b) ≈ 1 rtol = 1E-3
field = MonochromaticSpatialBeamRadialSymmetric(Float64, Forward, fibre, r, 1500e-9)
@test intensity(field) ≈ 1 rtol = 1E-3

nsx = LinRange(-0.5, 0.5, 1000)
field = MonochromaticAngularSpectrum_gaussian(Forward, nsx, nsx, mfd, 1500E-9, Medium(1.0), frames[1])
(a, b) = light_interaction(fibre, field)
@test intensity(b) ≈ 1 rtol = 1E-3
field = MonochromaticAngularSpectrum(Float64, Forward, fibre, nsx, nsx, 1500e-9)
@test intensity(field) ≈ 1 rtol = 1E-3

nsr = LinRange(0, 0.5, 100000)
field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, mfd, 1500E-9, Medium(1.0), frames[1])
(a, b) = light_interaction(fibre, field)
@test intensity(b) ≈ 1 rtol = 1E-3
field = MonochromaticAngularSpectrumRadialSymmetric(Float64, Forward, fibre, nsr, 1500e-9)
@test intensity(field) ≈ 1 rtol = 1E-3

sr = LinRange(0, 0.5, 10000)
λ = 1500E-9
mfd = 10E-6
re(e, sr, mfd) = e * mfd * exp(-mfd^2 * (2π/ λ)^2 * sr^2 / 16) * r * 2π 
sum(2π / λ *  re.(1, sr, mfd) .* step(sr)).^2

mfd = 20e-6
re_spa(e, r, mfd) = exp(-4π * r^2 / mfd^2) * r * 2π
re_spa2(e, r, mfd) = exp(-4π * r^2 / mfd^2)

@test 8 / (mfd^2) * sum(abs2.(re_spa2.(1, r, mfd)) .* r .* 2π .* step(r)) ≈ 1 rtol = 1E-2
mfd = 20E-6
function detection_fun(r, e, mfd)
    8 / (mfd^2) * abs2(sum(e .* exp.(-4π .* r.^2 ./ mfd^2) .* r .* 2π .* step(r)))
end
r = LinRange(0, 4mfd, 10000)#
spa =  MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, mfd, λ, Medium(1.0), frames[1])
intensity(spa)
detection_fun(r, spa.e, mfd)


function detection_fun_ang(nsr, e, mfd, λ)
    k = 2π/λ
    mfd^2 / 4 * k^2 / 2π * abs2(sum(e .* exp.(-nsr.^2 .* k.^2 .* mfd^2 ./ 16) .* step(nsr) .* 2π .* nsr))
    # cons = mfd * √(1 / 32 / π^3) * k^2;
    # √(64π^6 * real(angsperef.n)) * cons 
end
sr = LinRange(0, 0.5, 10000)
λ = 1500E-9
mfd = 20E-6
field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, mfd, λ, Medium(1.0), frames[1])
intensity(field)
detection_fun_ang(nsr, field.e, mfd, λ) 

k = 2π/λ
mfd = 20E-6
mfd^2 / 4 * k^2 / 2π * sum(abs2.(exp.(-nsr.^2 .* k.^2 .* mfd^2 ./ 16)) .* step(nsr) .* 2π .* nsr)


nsr = LinRange(0,1,10000)[2:end]

f(i) = mfd^2 / 4 * k^2 * sum(abs2.(exp.(-i[1].^2 .* k.^2 .* mfd^2 ./ 16)) .* i) 
Jolab.hcubature(f, [0], [1])
function coupling_test(mfd, λ; rtol= 1E-5)
    r = LinRange(0, 3mfd, 10000)#
    fibre = SingleModeFibre(mfd, 1E-2, media, frames)
    Jolab.findmodes!(fibre, λ)
    mode = Jolab.modes(fibre, λ)[1]
    spa =  MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, mfd, λ, Medium(1.0), frames[1])
    
    bool_1 = isapprox(intensity(mode, Jolab.R_θ_λ, spa.mesh), 1, rtol = rtol)
    (tmp, f) = light_interaction(fibre, spa)

    bool_2 = isapprox(detection_fun(r, spa.e, mfd), intensity(f), rtol = rtol)
    
    nsr = LinRange(0, 0.5, 10000)
    angspe = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, mfd, λ, Medium(1.0), frames[1])

    bool_3 = isapprox(intensity(mode, Jolab.NSR_NSθ_λ, angspe.mesh), 1, rtol = rtol)
    intensity(mode, Jolab.NSR_NSθ_λ, angspe.mesh)
end
# coupling_test(10E-6, 1500E-9, rtol = 1E-3)
# @test all(coupling_test(10E-6, 1500E-9, rtol = 1E-3))
# @test all(coupling_test(20E-6, 400E-9, rtol = 1E-3))
# @test all(coupling_test(30E-6, 600E-9, rtol = 1E-3))
# @test all(coupling_test(40E-6, 100E-9, rtol = 1E-3))

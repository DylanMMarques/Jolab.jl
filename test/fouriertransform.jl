using Jolab, Test

sx = range(-.3, .3, length = 51)
x = range(-15E-6, 15E-6, length = 51)

ref = ReferenceFrame(0.,0,0,0,0)
λ = 1550E-9
dir = 1
angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,λ, 1., dir, ref)
spaceref = FieldSpace_gaussian(x, x, 10E-6,λ, 1., dir, ref)
fourier = Jolab.FourierTransform(x, x, sx, sx)

coef = Jolab.coefficient_general(fourier, spaceref)

fourier2 = Jolab.FourierTransform(x, x, sx, sx, ref)
coef2 = Jolab.coefficient_general(fourier2, spaceref)
@test isapprox(coef.t₁₂, coef2.t₁₂, rtol = 1E-2)
@test isapprox(coef.t₂₁, coef2.t₂₁, rtol = 1E-2)
# stop

#Checking space to angular spectrum
(tmp, angspe) = lightinteraction(coef, spaceref)
@test isapprox(angsperef.e_SXY, angspe.e_SXY, rtol = .5E-2)

angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,λ, 1., -dir, ref)
(space, tmp) = lightinteraction(coef, angsperef)
@test isapprox(spaceref.e_SXY, space.e_SXY, rtol = .5E-2)

spaceref = FieldSpace_gaussian(x, x, 10E-6,λ, 1., -dir, ref)
coef = Jolab.coefficient_general(fourier, spaceref)

fourier2 = Jolab.FourierTransform(x, x, sx, sx, ref)
coef2 = Jolab.coefficient_general(fourier2, spaceref)
@test isapprox(coef.t₁₂, coef2.t₁₂, rtol = 1E-5)
@test isapprox(coef.t₂₁, coef2.t₂₁, rtol = 1E-5)

(angspe, tmp) = lightinteraction(coef, spaceref)
@test isapprox(angsperef.e_SXY, angspe.e_SXY, rtol = .5E-2)

angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,λ, 1., dir, ref)
(tmp, space) = lightinteraction(coef, angsperef)
@test isapprox(spaceref.e_SXY, space.e_SXY, rtol = .5E-2)

# Checking angular spectrum to space
dir = 1
angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,λ, 1., dir, ref)
spaceref = FieldSpace_gaussian(x, x, 10E-6,λ, 1., dir, ref)
fourier = Jolab.FourierTransform(x, x, sx, sx)

(tmp, space) = lightinteraction(fourier, angsperef)
@test isapprox(spaceref.e_SXY, space.e_SXY, rtol = .5E-2)

spaceref = FieldSpace_gaussian(x, x, 10E-6,λ, 1., -dir, ref)
(angspe, tmp) = lightinteraction(fourier, spaceref)
@test isapprox(angsperef.e_SXY, angspe.e_SXY, rtol = .5E-2)

dir = -1
angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,λ, 1., dir, ref)
spaceref = FieldSpace_gaussian(x, x, 10E-6,λ, 1., dir, ref)
fourier = Jolab.FourierTransform(x, x, sx, sx)

(space, tmp) = lightinteraction(fourier, angsperef)
@test isapprox(spaceref.e_SXY, space.e_SXY, rtol = .5E-2)

spaceref = FieldSpace_gaussian(x, x, 10E-6,λ, 1., -dir, ref)
(tmp, angspe) = lightinteraction(fourier, spaceref)
@test isapprox(angsperef.e_SXY, angspe.e_SXY, rtol = .5E-2)

# Checking field propagation
sx = range(-.3, .3, length = 50)
sy = range(-.3, .3, length = 50)
x = range(0E-6, 20E-6, length = 10)
y = range(0E-6, 1E-12, length = 2)

z = 1E-5
ref = ReferenceFrame(0.,0,0,0,0)
ref2 = ReferenceFrame(0.,0,z,0,0)
dir = 1
angsperef = FieldAngularSpectrum_gaussian(sx, sy, 10E-6, λ, 1., dir, ref)
fourier = Jolab.FourierTransform(x, y, sx, sy, ref2)

(tmp, space) = lightinteraction(fourier, angsperef)

# Guassian beam propagation for reference
zr(ω0) = π * (ω0)^2 / λ
ω(z, ω0) = ω0 * √(1 + (z / zr(ω0))^2)
ϕ(z,ω0) = atan(z / zr(ω0))
R(z,ω0) = z * (1 + (zr(ω0) / z)^2)
gauss(ω0, z, r) = ω0 / ω(z,ω0) * exp(-r^2 / (ω(z,ω0))^2) * exp(im * (2π/λ * z + 2π/λ * r^2 / 2 / R(z,ω0) - ϕ(z,ω0)))

refbeam = gauss.(5E-6, z, x)
refbeam .*= maximum(abs.(space.e_SXY[1,:,1])) ./ maximum(abs.(refbeam))
@test isapprox(refbeam, space.e_SXY[1,:,1], rtol = .5E-2)

dir = -1
angsperef = FieldAngularSpectrum_gaussian(sx, sy, 10E-6, λ, 1., dir, ref)
fourier = Jolab.FourierTransform(x, y, sx, sy, ref2)

(space, tmp) = lightinteraction(fourier, angsperef)
refbeam = gauss.(5E-6, -z, x)
refbeam .*= maximum(abs.(space.e_SXY[1,:,1])) ./ maximum(abs.(refbeam))
@test isapprox(refbeam, space.e_SXY[1,:,1], rtol = .5E-2)
return true

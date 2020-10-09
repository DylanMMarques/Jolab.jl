using Jolab
using Test

#if false
sx = range(-.3, .3, length = 50)
x = range(-15E-6, 15E-6, length = 50)

ref = ReferenceFrame(0.,0,0,0,0)
λ = 1550E-9
dir = 1
angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,λ, 1., dir, ref)
spaceref = FieldSpace_gaussian(x, x, 10E-6,λ, 1., dir, ref)
fourier = Jolab.FourierTransform(x, x, sx, sx)
coef = Jolab.coefficient_general(fourier, spaceref)

#Checking space to angular spectrum
(tmp, angspe) = lightinteraction(coef, spaceref)
@test isapprox(angsperef.e_SXY, angspe.e_SXY, rtol = .5E-2)

angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,λ, 1., -dir, ref)
(space, tmp) = lightinteraction(coef, angsperef)
@test isapprox(spaceref.e_SXY, space.e_SXY, rtol = .5E-2)

spaceref = FieldSpace_gaussian(x, x, 10E-6,λ, 1., -dir, ref)
coef = Jolab.coefficient_general(fourier, spaceref)
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
return true

#end
if false
sx = range(-.3, .3, length = 100)
x = range(-200E-6, 200E-6, length = 100)

ref = ReferenceFrame(0.,0,0,0,0)
ref2 = ReferenceFrame(0.,0,1E-3,0,0)
λ = 1550E-9
dir = 1
angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,λ, 1., dir, ref)
spaceref = FieldSpace_gaussian(x, x, 10E-6,λ, 1., dir, ref)
fourier = Jolab.FourierTransform(x, x, sx, sx, ref2)
coef = Jolab.coefficient_general(fourier, angsperef)

@time (tmp, space) = lightinteraction(coef, angsperef)
end

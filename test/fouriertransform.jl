using Jolab
using FFTW
using Test

sx = range(-1, 1, length = 51)
ref = ReferenceFrame(0.,0,0,0,0)
x = fftshift(FFTW.fftfreq(length(sx), 1 / (sx[2] - sx[1]))) * angsperef.λ;

angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,1550E-9, 1., 1, ref)
spaceref = FieldSpace_gaussian(x, x, 10E-6,1550E-9, 1., 1, ref)

fft = Jolab.FourierTransform()
coef = Jolab.coefficient_general(fft, angsperef)
(tmp, space) = lightinteraction(coef, angsperef)
@test isapprox(space.e_SXY, spaceref.e_SXY, rtol = 1E-5)

coef = Jolab.coefficient_general(fft, spaceref)
(tmp, angspe) = lightinteraction(coef, spaceref)
@test isapprox(angspe.e_SXY, angsperef.e_SXY, rtol = 1E-5)


angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,1550E-9, 1., -1, ref)
spaceref = FieldSpace_gaussian(x, x, 10E-6,1550E-9, 1., -1, ref)

fft = Jolab.FourierTransform()
coef = Jolab.coefficient_general(fft, angsperef)
(space, tmp) = lightinteraction(coef, angsperef)
@test isapprox(space.e_SXY, spaceref.e_SXY, rtol = 1E-5)

coef = Jolab.coefficient_general(fft, spaceref)
(angspe, tmp) = lightinteraction(coef, spaceref)
@test isapprox(angspe.e_SXY, angsperef.e_SXY, rtol = 1E-5)

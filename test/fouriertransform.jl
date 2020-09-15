using Jolab
using FFTW
using Test

sx = range(-1, 1, length = 60)
ref = ReferenceFrame(0.,0,0,0,0)
λ = 1550E-9
x = fftshift(FFTW.fftfreq(length(sx), 1 / (sx[2] - sx[1]))) * λ;
sx = fftshift(FFTW.fftfreq(length(x), 1 / (x[2] - x[1]))) * λ;
dir = 1
angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,λ, 1., dir, ref)
spaceref = FieldSpace_gaussian(x, x, 10E-6,1550E-9, 1., dir, ref)
fft = Jolab.FourierTransform()
#analytical evaluation of the fourrier transform
coef = Jolab.coefficient_general(fft, angsperef)
(tmp, space) = lightinteraction(coef, angsperef)
#fft algorithm
coef_fft = Jolab.coefficient_specific(fft, angsperef)
(tmp, space_fft) = lightinteraction(coef_fft, angsperef)
coef_fft = Jolab.coefficient_specific(fft, space_fft)
(tmp, angspe_fft) = lightinteraction(coef_fft, space_fft)
#Test if 2 fft gives the initial field
@test isapprox(angsperef.e_SXY, angspe_fft.e_SXY, rtol = 1E-5)
#Test if fourier transform gives the same as the analytical version
@test isapprox(space.e_SXY, spaceref.e_SXY, rtol = 1E-5)
@test isapprox(space.e_SXY, space_fft.e_SXY, rtol = 1E-5)
coef = Jolab.coefficient_general(fft, spaceref)
(tmp, angspe) = lightinteraction(coef, spaceref)
#Test if fourier transform gives the same result as the analytical of a gausian beam
@test isapprox(angspe.e_SXY, angsperef.e_SXY, rtol = 1E-5)

dir = -1
angsperef = FieldAngularSpectrum_gaussian(sx, sx, 10E-6,λ, 1., dir, ref)
spaceref = FieldSpace_gaussian(x, x, 10E-6,1550E-9, 1., dir, ref)
fft = Jolab.FourierTransform()
#analytical evaluation of the fourrier transform
coef = Jolab.coefficient_general(fft, angsperef)
(space, tmp) = lightinteraction(coef, angsperef)
#fft algorithm
coef_fft = Jolab.coefficient_specific(fft, angsperef)
(space_fft, tmp) = lightinteraction(coef_fft, angsperef)
coef_fft = Jolab.coefficient_specific(fft, space_fft)
(angspe_fft, tmp) = lightinteraction(coef_fft, space_fft)
#Test if 2 fft gives the initial field
@test isapprox(angsperef.e_SXY, angspe_fft.e_SXY, rtol = 1E-5)
#Test if fourier transform gives the same as the analytical version
@test isapprox(space.e_SXY, spaceref.e_SXY, rtol = 1E-5)
@test isapprox(space.e_SXY, space_fft.e_SXY, rtol = 1E-5)
coef = Jolab.coefficient_general(fft, spaceref)
(angspe, tmp) = lightinteraction(coef, spaceref)
#Test if fourier transform gives the same result as the analytical of a gausian beam
@test isapprox(angspe.e_SXY, angsperef.e_SXY, rtol = 1E-5)
return true

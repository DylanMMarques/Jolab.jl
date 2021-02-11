using Jolab, Test

sx = range(-.3, .3, length = 51)
x = range(-15E-6, 15E-6, length = 51)

ref = ReferenceFrame(0.,0,0,0,0)
λ0 = 1550E-9
dir = 1
angspe = FieldAngularSpectrum_gaussian(sx, sx, 10E-6, λ0, 1., dir, ref)

t = range(-.5E-9, .5E-9, length = 20)
iT = @. exp(-(t - .0E-9)^2 / .1E-9^2)
λ = λ0 .+ range(-5E-11, 5E-11, length = 41)
@time poly = Jolab.FieldPolychromatic_monochromaticpulse(angspe, t, iT, λ)

mls = MultilayerStructure([1., 2, 3, 1], [100E-9, 150E-9], ReferenceFrame(0.,0,10E-9,0.,0.))
@time (lfield, rfield) = lightinteraction(mls, poly)

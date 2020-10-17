using Jolab

nsx = range(-0.006, 0.006, length = 50)
x = range(-3E-3, 3E-3, length = 51)

mirror = Mirror(.99, 1, 1., ReferenceFrame(0,0,0))

for dir in [-1, 1]
    field = FieldSpace_gaussian(x, x, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,0.))
    lightinteraction(mirror, field)

    field = FieldAngularSpectrum_gaussian(nsx, nsx, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,0.))
    lightinteraction(mirror, field)
end
return true

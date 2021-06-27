using Jolab, Test

nsx = range(-0.006, 0.006, length = 51)
x = range(-3E-3, 3E-3, length = 51)

mirror = Mirror(.5, 1, 3., ReferenceFrame(0,0,0))

field = FieldSpaceScalar_gaussian(x, x, 2.5E-3, 1550E-9, 1, 1, ReferenceFrame(0,0,0.))
(field1, field2) = lightinteraction(mirror, field)
@test isapprox(intensity(field1) + intensity(field2), 1, atol = 1E-4)

field = FieldAngularSpectrumScalar_gaussian(nsx, nsx, 2.5E-3, 1550E-9, 1, 1, ReferenceFrame(0,0,0.))
(field1, field2) = lightinteraction(mirror, field)
@test isapprox(intensity(field1) + intensity(field2), 1, atol = 1E-4)

field = FieldSpaceScalar_gaussian(x, x, 2.5E-3, 1550E-9, 3, -1, ReferenceFrame(0,0,0.))
(field1, field2) = lightinteraction(mirror, field)
@test isapprox(intensity(field1) + intensity(field2), 1, atol = 1E-4)

field = FieldAngularSpectrumScalar_gaussian(nsx, nsx, 2.5E-3, 1550E-9, 3, -1, ReferenceFrame(0,0,0.))
(field1, field2) = lightinteraction(mirror, field)
@test isapprox(intensity(field1) + intensity(field2), 1, atol = 1E-4)

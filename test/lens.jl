using Jolab

nsx = range(-0.006, 0.006, length = 50)
x = range(-3E-3, 3E-3, length = 51)

mirror = Mirror(.99, 1, 1., ReferenceFrame(0,0,0))
lens = Lens(10E-3,1,ReferenceFrame(0,0,10E-3))

dir = 1
field = FieldSpaceScalar_gaussian(x, x, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,0.))
(fieldr, fieldt) = lightinteraction(lens, field)

field = FieldAngularSpectrumScalar_gaussian(nsx, nsx, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,0.))
(fieldr, fieldt) = lightinteraction(lens, field)

dir = -1
field = FieldSpaceScalar_gaussian(x, x, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,20E-3))
(fieldr, fieldt) = lightinteraction(lens, field)

field = FieldAngularSpectrumScalar_gaussian(nsx, nsx, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,20E-3))
(fieldr, fieldt) = lightinteraction(lens, field)

return true

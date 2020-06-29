using Jolab

z(x,y) = 10E-9
sx = range(-.1, .1, length = 10)
sy = range(-.1, .1, length = 10)
rmls = RoughInterface(1, 2, z, ReferenceFrame(0,0,10E-9,0,0))

field = FieldAngularSpectrum_gaussian(sx, sx, 50E-6, 1500E-9, 1., 1, ReferenceFrame(0.,0,0,0,0))

Jolab.coefficient_general(rmls, field)
# lightinteraction(rmls, field)

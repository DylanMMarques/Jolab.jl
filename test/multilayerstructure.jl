using Jolab, Test

mls = MultilayerStructure([1., 2, 3, 1], [100E-9, 150E-9], ReferenceFrame(0.,0,10E-9,0.,0.))
sx = range(-0.1, 0.1, length = 10)
field = FieldAngularSpectrum_gaussian(sx, sx, 50E-6, 1500E-9, 1., 1, ReferenceFrame(0.,0,0,0,0))
(fieldr, fieldt) = lightinteraction(mls, field)
return true

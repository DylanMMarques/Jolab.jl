using Jolab

x = range(-100E-6, 100E-6, length = 10)

mirror = Mirror(.99, 1, 1., ReferenceFrame(0,0,0))
mmf = CircularStepIndexFibre(50E-6,.1,1.44,1,ReferenceFrame(0,0,0.), 1, 1.)

dir = 1
field = FieldSpace_gaussian(x, x, 50E-6, 1550E-9, 1., dir, ReferenceFrame(50E-6,0,0.))
(fieldr, fieldt) = lightinteraction(mmf, field)

scat = Jolab.coefficient_general(mmf, field)
(fieldr, fieldt) = lightinteraction(scat, field)

dir = -1
field = FieldSpace_gaussian(x, x, 50E-6, 1550E-9, 1., dir, ReferenceFrame(50E-6,0,1))
(fieldr, fieldt) = lightinteraction(mmf, field)

scat = Jolab.coefficient_general(mmf, field)
(fieldr, fieldt) = lightinteraction(scat, field)

# dir = -1
# field = FieldSpace_gaussian(x, x, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,1E-2))
# lightinteraction(mmf, field)
return true

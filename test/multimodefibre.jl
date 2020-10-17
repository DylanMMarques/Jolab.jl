using Jolab

x = range(-3E-3, 3E-3, length = 51)

mirror = Mirror(.99, 1, 1., ReferenceFrame(0,0,0))
mmf = CircularStepIndexFibre(50E-6,.2,1.5,1,ReferenceFrame(0,0,0.), 1, 1E-2)

dir = 1
field = FieldSpace_gaussian(x, x, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,0.))
lightinteraction(mmf, field)

# dir = -1
# field = FieldSpace_gaussian(x, x, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,1E-2))
# lightinteraction(mmf, field)
return true

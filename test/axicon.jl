using Jolab
axicon = Axicon(0.003, ReferenceFrame(0.,0.,0))

# Create the gaussian field incident upon the axicon
x = range(-3E-3, 3E-3, length = 51)
nsx = range(-0.006, 0.006, length = 50)

for dir in [-1, 1]
    dir = 1
    field = FieldSpace_gaussian(x, x, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,0.))

    # calculate the transmission matrix of the axicon
    coef = Jolab.coefficient_general(axicon, field)
    bessel = lightinteraction(axicon, field)

    axiconfourier = AxiconFourier(x, x, nsx, nsx, 0.003, ReferenceFrame(0.,0,0))

    bessel = lightinteraction(axiconfourier, field)
    coef = Jolab.coefficient_general(axiconfourier, field)

    b = lightinteraction(axiconfourier, bessel[dir == 1 ? 2 : 1])
    coef = Jolab.coefficient_general(axiconfourier, bessel[dir==1 ? 2 : 1])
end
return true

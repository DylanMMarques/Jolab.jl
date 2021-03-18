using Jolab, Test

# Create the gaussian field incident upon the axicon
x = range(-3E-3, 3E-3, length = 51)
nsx = range(-0.006, 0.006, length = 50)

axicon = Axicon(0.003, ReferenceFrame(0.,0.,0))
axiconfourier = AxiconFourier(x, x, nsx, nsx, 0.003, ReferenceFrame(0.,0,0))

for dir in [-1, 1]
    field_l = FieldSpaceScalar_gaussian(x, x, 2.5E-3, 1550E-9, 1, dir, ReferenceFrame(0,0,0.))

    # calculate the transmission matrix of the axicon
    coef_l = Jolab.coefficient_general(axicon, field_l)
    bessel_l = lightinteraction(axicon, field_l)


    bessel_l = lightinteraction(axiconfourier, field_l)
    coef_l = Jolab.coefficient_general(axiconfourier, field_l)

    b_l = lightinteraction(axiconfourier, bessel_l[dir == 1 ? 2 : 1])
    coef_l = Jolab.coefficient_general(axiconfourier, bessel_l[dir==1 ? 2 : 1])
end
return true

using Jolab, Test

z(x,y) = 1E-9
sx = range(-.1, .1, length = 10)
sy = range(-.1, .1, length = 10)
rmls = RoughInterface(1, 2., z, ReferenceFrame(0.,0,0,0,0))
mls = MultilayerStructure([1., 1, 2], [z(0,0)], ReferenceFrame(0.,0,0,0,0))
field = FieldAngularSpectrum_gaussian(sx, sx, 50E-6, 1500E-9, 1., 1, ReferenceFrame(0.,0,10E-9,0,0))

(lfield1, rfield1) = lightinteraction(rmls, field)
(lfield2, rfield2) = lightinteraction(mls, field)

@test isapprox(lfield1.e_SXY[1,5,5], lfield2.e_SXY[1,5,5], rtol = 1E-3)

(coef) = Jolab.coefficient_general(rmls, field)
(lfield3, rfield3) = lightinteraction(coef, field)
@test isapprox(lfield1.e_SXY[1,5,5], lfield3.e_SXY[1,5,5], rtol = 1E-5)

changereferenceframe!(rfield2, ReferenceFrame(0.,0,0,0,0))
@test isapprox(rfield1.e_SXY[1,5,5], rfield2.e_SXY[1,5,5], rtol = 1E-3)


return true

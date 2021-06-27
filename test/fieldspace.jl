using Jolab, Test


x = range(-100e-6, 100e-6, length = 101)
space = FieldSpaceScalar_gaussian(x, x, 50e-6, 1500e-9, 1, 1, ReferenceFrame(0,0,0,0,0))
@test isapprox(intensity(space), 1, atol = 1e-5)

space = FieldSpaceScalar_uniform(x, x, 1500e-9, 1, 1, ReferenceFrame(0,0,0,0,0))
@test isapprox(intensity(space), 1, atol = 1e-5)

x = range(0, 100e-6, length = 1001)
space = Jolab.FieldSpaceScalarRadialSymmetric_uniform(x, 1500e-9, 1, 1, ReferenceFrame(0,0,0,0,0))
@test isapprox(intensity(space), 1, atol = 1e-5)

space = Jolab.FieldSpaceScalarRadialSymmetric_gaussian(x, 50E-6, 1500e-9, 1, 1, ReferenceFrame(0,0,0,0,0))
@test isapprox(intensity(space), 1, atol = 1e-5)

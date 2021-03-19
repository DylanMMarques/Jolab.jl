using Jolab, Test, LinearAlgebra

z_f(x,y) = 1E-9
sx = range(-.1, .1, length = 10)
sy = range(-.1, .1, length = 10)
rmls = Jolab.RoughMultilayerStructure([1, 2.],zeros(0), [z_f], ReferenceFrame(0.,0,0,0,0))
mls = MultilayerStructure([1., 1, 2], [z_f(0,0)], ReferenceFrame(0.,0,0,0,0))
field = FieldAngularSpectrumScalar_gaussian(sx, sx, 50E-6, 1500E-9, 1., 1, ReferenceFrame(0.,0,10E-9,0,0))

(lfield1, rfield1) = lightinteraction(rmls, field)
(lfield2, rfield2) = lightinteraction(mls, field)
changereferenceframe!(lfield2, lfield1.ref)
changereferenceframe!(rfield2, rfield1.ref)

@test isapprox(lfield1.e_SXY[55], lfield2.e_SXY[55], rtol = 1E-4)
@test isapprox(rfield1.e_SXY[55], rfield2.e_SXY[55], rtol = 1E-4)

(coef) = Jolab.coefficient_general(rmls, field)
(lfield3, rfield3) = lightinteraction(coef, field)
@test isapprox(lfield1.e_SXY[55], lfield3.e_SXY[55], rtol = 1E-10)
@test isapprox(rfield1.e_SXY[55], rfield3.e_SXY[55], rtol = 1E-10)

rmls = Jolab.RoughMultilayerStructure([2, 1.],zeros(0), [z_f], ReferenceFrame(0.,0,0,0,0))
mls = MultilayerStructure([2., 2, 1], [z_f(0,0)], ReferenceFrame(0.,0,0,0,0))
field = FieldAngularSpectrumScalar_gaussian(sx, sx, 50E-6, 1500E-9, 1., -1, ReferenceFrame(0.,0,10E-9,0,0))

(lfield1, rfield1) = lightinteraction(rmls, field)
(lfield2, rfield2) = lightinteraction(mls, field)
changereferenceframe!(lfield2, lfield1.ref)
changereferenceframe!(rfield2, rfield1.ref)

@test isapprox(rfield1.e_SXY[55], rfield2.e_SXY[55], rtol = 1E-4)
@test isapprox(lfield1.e_SXY[55], lfield2.e_SXY[55], rtol = 1E-4)

(coef) = Jolab.coefficient_general(rmls, field)
coef2 = Jolab.coefficient_general(mls, field)
(lfield3, rfield3) = lightinteraction(coef, field)
@test isapprox(lfield1.e_SXY[55], lfield3.e_SXY[55], rtol = 1E-10)
@test isapprox(rfield1.e_SXY[55], rfield3.e_SXY[55], rtol = 1E-10)

## Compairison agaisnt arctan perturbation on a filter
data = [0.00408221220653889,
    0.0011402066408047376,
    0.000684671687514954,
    0.0026399636502819055,
    0.006679540081315143,
    0.01230492617854357]

ref1 = ReferenceFrame(0,0,0,0,0.)
b = 1500E-9 / 2
ref2 = ReferenceFrame(0,0,b/1.5,0,0.)
ref3 = ReferenceFrame(0,0,ref2.z + b,0,0.)
ref4 = ReferenceFrame(0,0,ref3.z + b/1.5,0,0.)
sx = range(-.7, .7, length = 64+1)
sx = sx[2:end]
λ = range(1450E-9, 1550E-9, length = 6)
zp(x,y) = 5E-9 * (atan(-1E7 * x) / π + .5)
zm(x,y) = -5E-9 * (atan(-1E7 * x) / π + .5)

RoughInterface(n1,n2,z, ref) = RoughMultilayerStructure([n1, n2], zeros(0), [z], ref)

rmls = [RoughInterface(1, 1.5, zp, ref1), RoughInterface(1.5, 1, zm, ref2), RoughInterface(1,1.5,zp,ref3), RoughInterface(1.5,1,zm,ref4)]
int = zeros(length(λ))

for i in 1:length(λ)
    field = FieldAngularSpectrumScalar_gaussian(sx, sx, 10E-6, λ[i], 1, 1, ref1)
    (fieldr, fieldt) = lightinteraction_recursivegridded(rmls, field, rtol = 1E-9)
    int[i] = intensity(fieldr)
end

@test all(isapprox.(int, data, rtol = 1E-8))
return true

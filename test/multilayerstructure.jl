using Jolab, Test

## Test if the scattering matrix approach and lightinteraction do the same
mls = MultilayerStructure([1., 2, 3, 1], [100E-9, 150E-9], ReferenceFrame(0.,0,10E-9,0.,0.))
sx = range(-0.1, 0.1, length = 11)
field = FieldAngularSpectrumScalar_gaussian(sx, sx, 50E-6, 1500E-9, 1., 1, ReferenceFrame(0.,0,0,0,0))
(fieldr, fieldt) = lightinteraction(mls, field)

coef = Jolab.coefficient_general(mls, field)
(fieldr2, fieldt2) = lightinteraction(coef, field)
changereferenceframe!(fieldr2, fieldr.ref)
@test all(isapprox.(fieldt2.e_SXY, fieldt.e_SXY, rtol = 1E-12))
@test all(isapprox.(fieldr2.e_SXY, fieldr.e_SXY, rtol = 1E-12))

field = FieldAngularSpectrumScalar_gaussian(sx, sx, 50E-6, 1500E-9, 1., -1, ReferenceFrame(0.,0,0,0,0))
(fieldr, fieldt) = lightinteraction(mls, field)

coef = Jolab.coefficient_general(mls, field)
(fieldr2, fieldt2) = lightinteraction(coef, field)
changereferenceframe!(fieldr2, fieldr.ref)
changereferenceframe!(fieldt2, fieldt.ref)
@test all(isapprox.(fieldt2.e_SXY, fieldt.e_SXY, rtol = 1E-12))
@test all(isapprox.(fieldr2.e_SXY, fieldr.e_SXY, rtol = 1E-12))

h = [100E-9, 200E-9, 300E-9]
n = [2, 3, 4+im, 5, 6]

## Test plane wave incident fresnel coefficient (s waves)
sx = range(.1/√2, .1/√2, length = 2)
mls = MultilayerStructure(n, h, ReferenceFrame(0.,0,0,0.,0.))
field = FieldAngularSpectrumScalar_gaussian(sx, sx, 50E-6, 1550E-9, n[1], 1, ReferenceFrame(0.,0,0,0,0))

coef = Jolab.coefficient_general(mls, field)
@test isapprox(coef.r₁₂[1,1], 0.013694713935881 - 0.060293291844495im, rtol = 1E-12)
@test isapprox(coef.t₁₂[1,1], -0.115673040865732 - 0.229890263554481im, rtol = 1E-12)
@test isapprox(coef.r₂₁[1,1], 0.143881724333044 - 0.098603264296048im, rtol = 1E-12)
@test isapprox(coef.t₂₁[1,1], -0.347405450713043 - 0.690438584712133im, rtol = 1E-12)

# Test imaginary wave input frenes coefficient (s waves)
sx = range(1.5/√2, 1.5/√2, length = 2)
mls = MultilayerStructure(n, h, ReferenceFrame(0.,0,0,0.,0.))
field = FieldAngularSpectrumScalar_gaussian(sx, sx, 50E-6, 1550E-9, n[1], 1, ReferenceFrame(0.,0,0,0,0))

coef = Jolab.coefficient_general(mls, field)
@test isapprox(coef.r₁₂[1,1], -0.150961489254006 - 0.142653073710589im, rtol = 1E-12)
@test isapprox(coef.t₁₂[1,1], -0.181631767295952 - 0.079386334674067im, rtol = 1E-12)
@test isapprox(coef.r₂₁[1,1], 0.120595195598718 - 0.142697052795419im, rtol = 1E-12)
@test isapprox(coef.t₂₁[1,1], -0.797645047292507 - 0.348629084097925im, rtol = 1E-12)

# Test imaginary wave input frenes coefficient (s waves)
h = [100E-9, 200E-9, 300E-9]
n = [2+im, 3, 4+im, 5, 6]

sx = range(1.5/√2, 1.5/√2, length = 2)
mls = MultilayerStructure(n, h, ReferenceFrame(0.,0,0,0.,0.))
field = FieldAngularSpectrumScalar_gaussian(sx, sx, 50E-6, 1550E-9, n[1], 1, ReferenceFrame(0.,0,0,0,0))

coef = Jolab.coefficient_general(mls, field)
@test isapprox(coef.r₁₂[1,1], 0.064506848490957 + 0.203429135269147im, rtol = 1E-12)
@test isapprox(coef.t₁₂[1,1], -0.179599985581118 - 0.173227690399885im, rtol = 1E-12)
@test isapprox(coef.r₂₁[1,1], 0.142113245185404 - 0.209783024882484im, rtol = 1E-12)
@test isapprox(coef.t₂₁[1,1], -0.716578548332634 - 0.053567166641940im, rtol = 1E-12)

# Test imaginary wave input if wave is propagating backwards
field = FieldAngularSpectrumScalar_gaussian(sx, sx, 50E-6, 1550E-9, n[end], -1, ReferenceFrame(0.,0,sum(h),0,0))
coef = Jolab.coefficient_general(mls, field)
@test isapprox(coef.r₁₂[1,1], 0.064506848490957 + 0.203429135269147im, rtol = 1E-12)
@test isapprox(coef.t₁₂[1,1], -0.179599985581118 - 0.173227690399885im, rtol = 1E-12)
@test isapprox(coef.r₂₁[1,1], 0.142113245185404 - 0.209783024882484im, rtol = 1E-12)
@test isapprox(coef.t₂₁[1,1], -0.716578548332634 - 0.053567166641940im, rtol = 1E-12)
return true

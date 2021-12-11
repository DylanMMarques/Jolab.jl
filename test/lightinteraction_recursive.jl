using Jolab, Test

mls = MultilayerStructure([1., 2, 3+im, 1], [100E-9, 150E-9], ReferenceFrame(0.,0,0,0.,0.))
sx = range(-0.1, 0.1, length = 11)

for dir in [1, -1]
    local field
    local fieldt
    local fieldr
    field = FieldAngularSpectrumScalar_gaussian(sx, sx, 50E-6, 1500E-9, 1., dir, ReferenceFrame(0.,0,0,0,0))

    mls_r = [MultilayerStructure([1., 2], zeros(0), ReferenceFrame(0.,0,0,0.,0.)),
        MultilayerStructure([2., 3+im], zeros(0), ReferenceFrame(0.,0,100E-9,0.,0.)),
        MultilayerStructure([3+im, 1.], zeros(0), ReferenceFrame(0.,0,250E-9,0.,0.))]

    (fieldr, fieldt) = lightinteraction(mls, field)
    (fieldr2, fieldt2) = lightinteraction(mls_r, field)

    changereferenceframe!(fieldt2, fieldt.ref)
    changereferenceframe!(fieldr2, fieldr.ref)

    @test all(isapprox.(fieldr.e_SXY, fieldr2.e_SXY, rtol = 1E-12))
    @test all(isapprox.(fieldt.e_SXY, fieldt2.e_SXY, rtol = 1E-12))

    (fieldr2, fieldt2) = lightinteraction_recursivegridded(mls_r, field, rtol = 1E-9)

    changereferenceframe!(fieldt2, fieldt.ref)
    changereferenceframe!(fieldr2, fieldr.ref)

    @test all(isapprox.(fieldr.e_SXY, fieldr2.e_SXY, rtol = 1E-8))
    @test all(isapprox.(fieldt.e_SXY, fieldt2.e_SXY, rtol = 1E-8))
end

return true

using Jolab, Test

# test_opt fails due to kernel abstractions
nsr = range(0, 0.25, length = 1000)
r = range(0, 2E-3, length = 1000)
axicon = Axicon(5*π/180, Medium(1.44), Medium(1), ReferenceFrame((0,0,0), (0,0,0)); solver = (nsr = nsr, r = r))

field = MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, 1E-3, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
(rfield, tfield) = light_interaction(axicon, field)
@test tfield isa Jolab.MeshedBeam{<:Any, <:Any, Jolab.NSR_NSθ_λ}
@test isapprox(intensity(tfield), intensity(field), rtol = 1E-4)
@test iszero(intensity(rfield))
@test isapprox(nsr[findmax(abs, tfield.e)[2]], (axicon.axicon_medium.n - axicon.medium.n) * axicon.α, rtol = 1E-2)

(field_r_2, ori_field) = light_interaction(axicon, tfield)
@test iszero(intensity(field_r_2))
@test isapprox(intensity(ori_field), intensity(tfield), rtol = 1E-4)

# Check that the direction with peak amplitude is the expected one based on the axicon angle and refractive indices
@test isapprox(nsr[findmax(abs, tfield.e)[2]], (axicon.axicon_medium.n - axicon.medium.n) * axicon.α, rtol = 1E-2)

function test_autodiff(λ)
    nsr = range(0, 0.5, length = 100)
    axicon = Axicon(5*π/180, Medium(1.44), Medium(1), ReferenceFrame((0,0,0), (0,0,0)); solver = (nsr = nsr,))
    
    r = range(0, 2E-3, length = 100)
    field = MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, 1E-3, λ, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
    (rfield, tfield) = light_interaction(axicon, field)
    intensity(tfield)
end

autodiff(Enzyme.Forward, test_autodiff, Duplicated, Duplicated(1500E-9, 1.0))

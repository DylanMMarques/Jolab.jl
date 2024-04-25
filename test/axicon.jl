using Jolab

nsr = range(0, 0.5, length = 1000)
axicon = Axicon(5*π/180, Medium(1.44), Medium(1), ReferenceFrame((0,0,0), (0,0,0)); solver = (nsr = nsr,))

r = range(0, 2E-3, length = 10000)
field = MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, 1E-3, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
(rfield, tfield) = light_interaction(axicon, field)
@test tfield isa Jolab.MeshedBeam{<:Any, <:Any, <:Jolab.NSR_NSθ_t}
intensity(tfield)
intensity(field)

function test_autodiff(λ)
    nsr = range(0, 0.5, length = 100)
    axicon = Axicon(5*π/180, Medium(1.44), Medium(1), ReferenceFrame((0,0,0), (0,0,0)); solver = (nsr = nsr,))
    
    r = range(0, 2E-3, length = 100)
    field = MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, 1E-3, λ, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
    (rfield, tfield) = light_interaction(axicon, field)
    intensity.((rfield, tfield))
end
test_autodiff(1500E-9)
autodiff(Enzyme.Forward, test_autodiff, Duplicated, Duplicated(1500E-9, 1.0))
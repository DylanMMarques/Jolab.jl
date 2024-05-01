using Jolab

nsr = range(0, 0.5, length = 1000)
axicon = Axicon(5*π/180, Medium(1.44), Medium(1), ReferenceFrame((0,0,0), (0,0,0)); solver = (nsr = nsr,))

r = range(0, 2E-3, length = 10000)
field = MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, 1E-3, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
(rfield, tfield) = light_interaction(axicon, field)
@test tfield isa Jolab.MeshedBeam{<:Any, <:Any, <:Jolab.NSR_NSθ_t}
intensity(tfield) 
intensity(field)

mirror = Mirror(Medium.((1, 1.44)), ReferenceFrame((0,0,2f), (0,0,0)); reflectivity = 1)

r = range(0, 2E-3, length = 1000)
nsr = range(0, 0.2, length = 1000)
field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 10E-6, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
f = 10E-3
lens = Lens(f, 1, (Medium(1), Medium(1)), ReferenceFrame((0,0,f), (0,0,0)))
axicon = Axicon(5*π/180, Medium(1.44), Medium(1), ReferenceFrame((0,0,2f), (0,0,0)); solver = (nsr = nsr, r = r))

function coherent_sum(beam)
    f(e, ind) = e * (Jolab.volume(beam.mesh, ind) * beam.medium.n)
    abs2(mapreduce(f, +, vec(beam.e), eachindex(beam.e)))
end
coherent_sum(field)

fieldt = light_interaction(lens, field)[2]
fieldt2 = light_interaction(axicon, fieldt)[2]
rfield = light_interaction(mirror, fieldt2)[1]
rfield_tmp = deepcopy(rfield)
rfield_tmp.e .= rfield_tmp.e.^2
coherent_sum(rfield_tmp)
rfield2 = light_interaction(axicon, rfield)[1]
intensity(rfield2)
rfield3 = light_interaction(lens, rfield2)[1]
intensity(rfield3)
coherent_sum(rfield_tmp)


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
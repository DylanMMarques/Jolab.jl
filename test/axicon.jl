using Jolab, Test

nsr = range(0, 0.25, length = 1000)
axicon = Axicon(5*π/180, Medium(1.44), Medium(1), ReferenceFrame((0,0,0), (0,0,0)); solver = (nsr = nsr, r = r))

r = range(0, 2E-3, length = 1000)
field = MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, 1E-3, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
(rfield, tfield) = light_interaction(axicon, field)
@test tfield isa Jolab.MeshedBeam{<:Any, <:Any, Jolab.NSR_NSθ_λ}
@test isapprox(intensity(tfield), intensity(field), rtol = 1E-4)
@test iszero(intensity(rfield))
@test isapprox(nsr[findmax(abs, tfield.e)[2]], (axicon.axicon_medium.n - axicon.medium.n) * axicon.α, rtol = 1E-2)

(field_r_2, ori_field) = light_interaction(axicon, tfield)
@test iszero(intensity(field_r_2))
@test isapprox(intensity(ori_field), intensity(tfield), rtol = 1E-4)

# The amplitude of the field should be approximately the same as only a phase screen 
# was applied
@test abs(sum(abs.(field.e) .- abs.(ori_field.e))) <= sum(abs.(field.e)) * 1E-3

@test isapprox(nsr[findmax(abs, tfield.e)[2]], (axicon.axicon_medium.n - axicon.medium.n) * axicon.α, rtol = 1E-2)

focal = 10E-3
mirror = Mirror(Medium.((1, 1.44)), ReferenceFrame((0,0,2focal), (0,0,0)); reflectivity = 1)

r = range(0, 2E-3, length = 500)
nsr = range(0, 0.2, length = 500)
field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 10E-6, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
lens = Lens(focal, 1, (Medium(1), Medium(1)), ReferenceFrame((0.0,0.0,focal), (0.0,0.0,0.0)))
axicon = Axicon(5*π/180, Medium(1.44), Medium(1), ReferenceFrame((0,0,2focal), (0,0,0)); solver = (nsr = nsr, r = r))

function coherent_sum(beam)
    f(e, ind) = e * (Jolab.volume(beam.mesh, ind) * beam.medium.n)
    abs2(mapreduce(f, +, vec(beam.e), eachindex(beam.e)))
end
coherent_sum(field)
focal = 10E-3
h = 500E-6

function tmp_r(λ, z)
    fp = (Mirror(Medium.((1, 1.44)), ReferenceFrame((0,0,2focal + z), (0,0,0)); reflectivity = 0.95), 
        Mirror(Medium.((1.44, 1)), ReferenceFrame((0,0, 2focal + h + z), (0,0,0)); reflectivity = 0.95))
    mirror = Mirror(Medium.((1, 1.44)), ReferenceFrame((0,0,2focal + z), (0,0,0)); reflectivity = 1)
    field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 12.5E-6, λ, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
    fieldt = light_interaction(lens, field)[2]
    fieldt2 = light_interaction(axicon, fieldt)[2]

    ## With detection propagation
    ifield_tmp = deepcopy(fieldt2)
    ifield_tmp.e .= ifield_tmp.e.^2

    ref_field_tmp = translate_referenceframe(light_interaction(mirror, translate_referenceframe(ifield_tmp, (0,0,2focal + z))[2])[1], (0,0,2focal))[1]
    rfield_tmp = translate_referenceframe(light_interaction(fp, translate_referenceframe(ifield_tmp, (0,0,2focal + z))[2])[1], (0,0,2f))[1]
    
    rfield = light_interaction(fp, translate_referenceframe(fieldt2, (0,0,2focal + z))[2])[1]
    rfield2 = light_interaction(axicon, translate_referenceframe(rfield, (0,0,2focal))[1])[1]
    rfield3 = light_interaction(lens, rfield2)[1]
    rfield3.e .*= field.e
    
    ref_field = translate_referenceframe(light_interaction(mirror, translate_referenceframe(fieldt2, (0,0,2focal + z))[2])[1], (0,0,2f))[1]
    ref_field = light_interaction(axicon, ref_field)[1]
    ref_field = light_interaction(lens, ref_field)[1]
    ref_field.e .*= field.e
    coherent_sum(rfield3) / coherent_sum(ref_field), coherent_sum(rfield_tmp) / coherent_sum(ref_field_tmp) 
end

function reflectance_axi_mirror(fieldt2, z)
    mirror = Mirror(Medium.((1, 1.44)), ReferenceFrame((0,0,2focal + z), (0,0,0)); reflectivity = 1)
    fieldt = translate_referenceframe(fieldt2, (0,0,2focal + z))[2]
    ref_field = light_interaction(mirror, fieldt)[1]

    ref_field2 = light_interaction(axicon, translate_referenceframe(ref_field, (0,0,2focal))[1])[1]
    ref_field3 = light_interaction(lens, ref_field2)[1]
    ref_field3.e .*= field.e
    coherent_sum(ref_field3)
end
function field_focal_plane(λ)
    field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 12.5E-6, λ, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
    fieldt = light_interaction(lens, field)[2]
    fieldt2 = light_interaction(axicon, fieldt)[2]
end
z_vec = range(0, 5E-2, length = 2)
fieldi = field_focal_plane(1550E-9)
axial = map(z -> reflectance_axi_mirror(fieldi, z), z_vec)
z_max = z_vec[findmax(axial)[2]]
# λ_b = 1546.9568023494082nm
# tmp_r(ustrip(m, λ_b), z_max)
# λ = ustrip.(m, λ_b .+ range(-1nm, 1nm, length= 100))
# itf_tmp = tmap(λ -> tmp_r(λ, z_max), λ)

function test_autodiff(λ)
    nsr = range(0, 0.5, length = 100)
    axicon = Axicon(5*π/180, Medium(1.44), Medium(1), ReferenceFrame((0,0,0), (0,0,0)); solver = (nsr = nsr,))
    
    r = range(0, 2E-3, length = 100)
    field = MonochromaticSpatialBeamRadialSymmetric_gaussian(Forward, r, 1E-3, λ, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
    (rfield, tfield) = light_interaction(axicon, field)
    intensity(tfield)
end
@time test_autodiff(1500E-9)

autodiff(Enzyme.Forward, test_autodiff, Duplicated, Duplicated(1500E-9, 1.0))

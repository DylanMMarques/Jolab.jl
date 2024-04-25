mirror = Mirror((Medium(1), Medium(1)), ReferenceFrame((0,0,0), (0,0,0)); reflectivity = 0.99)
mirror2 = Mirror((Medium(1), Medium(1)), ReferenceFrame((0,0,0E-9), (0,0,0)); reflectivity = 0.99)

comps = (mirror, mirror2)

nsr = range(0, 0.5, length = 1000)
field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 10E-6, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))

light_interaction(mirror, field)
ScatteringMatrix(mirror, field)

scat = ScatteringMatrix(comps, field)

mls = DielectricStack(Medium.([1,2,3]), [100E-9], ReferenceFrame((0,0,0), (0,0,0)))
a = ScatteringMatrix(mls, field)

int_1 = DielectricStack(Medium.([1,2]), Float64[], ReferenceFrame((0,0,0), (0,0,0)))
int_2 = DielectricStack(Medium.([2,3]), Float64[], ReferenceFrame((0,0,100E-9), (0,0,0)))
comps = (int_1, int_2)
b = ScatteringMatrix(comps, field)

@test a.mat_itob ≈ b.mat_itob
@test a.mat_itof ≈ b.mat_itof


mls = DielectricStack(Medium.([1,2,3]), [100E-9], ReferenceFrame((0,0,0), (0,0,0)))

function tmp(λ)
    # mirror_1 = Mirror((Medium(1), Medium(1)), ReferenceFrame((0,0,0), (0,0,0)); reflectivity = 0.99)
    # mirror_2 = Mirror((Medium(1), Medium(1)), ReferenceFrame((0,0,100E-6), (0,0,0)); reflectivity = 0.99)
    # fp = (mirror_1, mirror_2)

    fp = DielectricStack(Medium.([1,2,3]), [100E-9], ReferenceFrame((0,0,0), (0,0,0)))
    nsr = LinRange(0, 0.5, 1000)
    field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 10E-6, λ, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
    scat = ScatteringMatrix(fp, field)
    res = light_interaction(scat, field)
    intensity(res[1])
end

tmp(1550E-9);
sens(λ) = autodiff_deferred(Enzyme.Reverse, tmp, DuplicatedNoNeed, Duplicated(λ, 1.0))[1]
sens(1550E-9)
dsens_dλ = autodiff_deferred(Enzyme.Forward, sens, Duplicated, Duplicated(1550E-9, 1.0))

using Jolab, Test

nsr = range(0, 0.5, length = 1000)
field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 10E-6, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))

mls = DielectricStack(Medium.([1,2,3]), [100E-6], ReferenceFrame((0,0,0), (0,0,0)))
mls_scat = ScatteringMatrix(mls, field)

int_1 = DielectricStack(Medium.([1,2]), Float64[], ReferenceFrame((0,0,0), (0,0,0)))
int_2 = DielectricStack(Medium.([2,3]), Float64[], ReferenceFrame((0,0,100E-6), (0,0,0)))
comps = (int_1, int_2)
inv_scat = ScatteringMatrix(comps, field)

@test mls_scat.mat_itob ≈ inv_scat.mat_itob
@test mls_scat.mat_itof ≈ inv_scat.mat_itof

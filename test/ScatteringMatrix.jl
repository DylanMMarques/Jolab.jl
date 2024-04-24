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
a[1] ≈ b[1]
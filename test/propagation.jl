using Jolab, Test
nsr = range(0, 0.5, length = 1000)
field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))

prop = Propagation((ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0.0,0.0,100E-9), (0,0,0))), Medium(2))

@test all((0, intensity(field)) .≈ intensity.(light_interaction(prop, field)))

field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
mls = DielectricStack([Medium(2), Medium(2), Medium(2)], [100E-9], ReferenceFrame((0,0,0), (0,0,0)))
field_mls = light_interaction(mls, field)
prop = Propagation((ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0.0,0.0,100E-9), (0,0,0))), Medium(2))
fields_prop = light_interaction(prop, field)
@test_opt light_interaction(mls, field)
@test all((field_mls) .≈ (fields_prop))

field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Backward, nsr, 10E-6, 1550E-9, Medium(2), ReferenceFrame((0,0,100E-9), (0,0,0E-9)))
mls = DielectricStack([Medium(2), Medium(2), Medium(2)], [100E-9], ReferenceFrame((0,0,0), (0,0,0)))
field_mls = light_interaction(mls, field)
prop = Propagation((ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0.0,0.0,100E-9), (0,0,0))), Medium(2))
fields_prop = light_interaction(prop, field)
@test all((field_mls) .≈ (fields_prop))

scat_prop = ScatteringMatrix(prop, field)
scat_mls = ScatteringMatrix(mls, field)
@test scat_mls.mat_itof ≈ scat_prop.mat_itof
@test scat_mls.mat_itob ≈ scat_prop.mat_itob

field = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, 10E-6, 1550E-9, Medium(2+im), ReferenceFrame((0,0,0), (0,0,0)))
mls = DielectricStack([Medium(2 + im), Medium(2 + im), Medium(2 + im)], [100E-9], ReferenceFrame((0,0,0), (0,0,0)))
field_mls = light_interaction(mls, field)
prop = Propagation((ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0.0,0.0,100E-9), (0,0,0))), Medium(2 + im))
fields_prop = light_interaction(prop, field)
@test all((field_mls) .≈ (fields_prop))

scat_prop = ScatteringMatrix(prop, field)
@test_opt ScatteringMatrix(prop, field)
scat_mls = ScatteringMatrix(mls, field)
@test scat_mls.mat_itof ≈ scat_prop.mat_itof
@test scat_mls.mat_itob ≈ scat_prop.mat_itob

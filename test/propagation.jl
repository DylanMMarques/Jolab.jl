using Jolab, Test, JET
nsr = range(0, 0.5, length = 100)
nst = range(0, 2π, length = 100)[1:(end-1)]
λ = 1550E-9
grid_ang_radial = Jolab.CylindricalGrid(Jolab.NSR_NSθ_λ, nsr, nst, λ)
field = Jolab.AngularSpectrum(
    Jolab.Forward,
    grid_ang_radial,
    Jolab.Gaussian(10E-6),
    Medium(2),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)

prop = Propagation(
    (ReferenceFrame((0, 0, 0), (0, 0, 0)), ReferenceFrame((0.0, 0.0, 100E-9), (0, 0, 0))),
    Medium(2),
)

@test all((0, intensity(field)) .≈ intensity.(light_interaction(prop, field)))

nsr = range(0, 0.5, length = 100)
nst = range(0, 2π, length = 100)[1:(end-1)]
grid_ang_radial = Jolab.CylindricalGrid(Jolab.NSR_NSθ_λ, nsr, nst, λ)
field = Jolab.AngularSpectrum(
    Jolab.Forward,
    grid_ang_radial,
    Jolab.Gaussian(10E-6),
    Medium(2),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)
mls = DielectricStack(
    [Medium(2), Medium(2), Medium(2)],
    [100E-9],
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)
field_mls = light_interaction(mls, field)
prop = Propagation(
    (ReferenceFrame((0, 0, 0), (0, 0, 0)), ReferenceFrame((0.0, 0.0, 100E-9), (0, 0, 0))),
    Medium(2),
)
fields_prop = light_interaction(prop, field)
@test_opt light_interaction(mls, field)
@test all((field_mls) .≈ (fields_prop))

nsr = range(0, 0.5, length = 100)
nst = range(0, 2π, length = 100)[1:(end-1)]
grid_ang_radial = Jolab.CylindricalGrid(Jolab.NSR_NSθ_λ, nsr, nst, λ)
field = Jolab.AngularSpectrum(
    Jolab.Backward,
    grid_ang_radial,
    Jolab.Gaussian(10E-6),
    Medium(2),
    ReferenceFrame((0, 0, 100E-9), (0, 0, 0E-9)),
)
mls = DielectricStack(
    [Medium(2), Medium(2), Medium(2)],
    [100E-9],
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)
field_mls = light_interaction(mls, field)
prop = Propagation(
    (ReferenceFrame((0, 0, 0), (0, 0, 0)), ReferenceFrame((0.0, 0.0, 100E-9), (0, 0, 0))),
    Medium(2),
)
fields_prop = light_interaction(prop, field)
@test all((field_mls) .≈ (fields_prop))

scat_prop = ScatteringMatrix(prop, field)
scat_mls = ScatteringMatrix(mls, field)
@test scat_mls.mat_itof ≈ scat_prop.mat_itof
@test scat_mls.mat_itob ≈ scat_prop.mat_itob

nsr = range(0, 0.5, length = 100)
nst = range(0, 2π, length = 100)[1:(end-1)]
grid_ang_radial = Jolab.CylindricalGrid(Jolab.NSR_NSθ_λ, nsr, nst, λ)
field = Jolab.AngularSpectrum(
    Jolab.Forward,
    grid_ang_radial,
    Jolab.Gaussian(10E-6),
    Medium(2+im),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)
mls = DielectricStack(
    [Medium(2 + im), Medium(2 + im), Medium(2 + im)],
    [100E-9],
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)
field_mls = light_interaction(mls, field)
prop = Propagation(
    (ReferenceFrame((0, 0, 0), (0, 0, 0)), ReferenceFrame((0.0, 0.0, 100E-9), (0, 0, 0))),
    Medium(2 + im),
)
fields_prop = light_interaction(prop, field)
@test all((field_mls) .≈ (fields_prop))

scat_prop = ScatteringMatrix(prop, field)
@test_opt ScatteringMatrix(prop, field)
scat_mls = ScatteringMatrix(mls, field)
@test scat_mls.mat_itof ≈ scat_prop.mat_itof
@test scat_mls.mat_itob ≈ scat_prop.mat_itob

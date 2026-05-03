using Jolab, Test

nsr = range(0, 0.5, length = 50)
nθ = 0.0
λ = 1550E-9
w0 = 10E-6
grid = Jolab.CylindricalGrid(Jolab.NSR_NSθ_λ, nsr, nθ, λ)
field = Jolab.AngularSpectrum(
    Jolab.Forward,
    grid,
    Jolab.Gaussian(w0),
    Medium(1),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)

mls = DielectricStack(Medium.([1, 2, 3]), [100E-6], ReferenceFrame((0, 0, 0), (0, 0, 0)))
@test_opt ScatteringMatrix(mls, field)
mls_scat = ScatteringMatrix(mls, field)

int_1 = DielectricStack(Medium.([1, 2]), Float64[], ReferenceFrame((0, 0, 0), (0, 0, 0)))
int_2 =
    DielectricStack(Medium.([2, 3]), Float64[], ReferenceFrame((0, 0, 100E-6), (0, 0, 0)))
comps = (int_1, int_2)
inv_scat = ScatteringMatrix(comps, field)

@time ScatteringMatrix(comps, field);

@test mls_scat.mat_itob ≈ inv_scat.mat_itob
@test mls_scat.mat_itof ≈ inv_scat.mat_itof

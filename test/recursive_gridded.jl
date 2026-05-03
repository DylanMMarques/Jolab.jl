using Jolab, Test

nsx = range(-0.1, 0.1, length = 64)
nsy = range(-0.1, 0.1, length = 64)
λ = 1550E-9
grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsy, λ)
beam = Jolab.AngularSpectrum(
    Float64,
    Forward,
    grid_ang,
    Jolab.Gaussian(50E-6),
    Medium(1.0),
    ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
);

mls = DielectricStack(
    Medium.([1, 1.5, 2]),
    [100E-9],
    ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
)
(r, t) = light_interaction(mls, beam)

prop = Propagation(mls.frames, Medium(1.5))
mls_layer = (
    DielectricStack(Medium.([1, 1.5]), zeros(0), ReferenceFrame((0, 0, 0.0), (0, 0, 0.0))),
    prop,
    DielectricStack(
        Medium.([1.5, 2]),
        zeros(0),
        ReferenceFrame((0, 0, 100E-9), (0, 0, 0.0)),
    ),
)

@test_opt (broken = VERSION < v"1.11") Jolab.lightinteraction_recursivegridded(
    mls_layer,
    beam,
    rtol = 1E-9,
)
(r2, t2) = Jolab.lightinteraction_recursivegridded(mls_layer, beam, rtol = 1E-9)

@test isapprox(r2, r, rtol = 1E-7)
@test isapprox(t2, t, rtol = 1E-7)


###
mls = DielectricStack(Medium.([1, 1, 1]), ones(1), ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)))
rec_mls = (
    Propagation(
        (ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)), ReferenceFrame((0, 0, 0.1), (0, 0, 0))),
        Medium(1.0),
    ),
    Propagation(
        (ReferenceFrame((0, 0, 0.1), (0, 0, 0.0)), ReferenceFrame((0, 0, 0.5), (0, 0, 0))),
        Medium(1.0),
    ),
    Propagation(
        (ReferenceFrame((0, 0, 0.5), (0, 0, 0.0)), ReferenceFrame((0, 0, 1.0), (0, 0, 0))),
        Medium(1.0),
    ),
)
(r3, t3) =
    Jolab.lightinteraction_recursivegridded(rec_mls, beam, rtol = 1E-20, printBool = false)
(r4, t4) = light_interaction(mls, beam)
@test isapprox(r3, r4, rtol = 1E-9)
@test isapprox(t3, t4, rtol = 1E-9)


mls = DielectricStack(
    Medium.([1, 1.5, 2, 2.5, 1]),
    [100E-9, 500E-9, 200E-9],
    ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
)
(r, t) = light_interaction(mls, beam)

topography_zero(x, y) = 0.0

mls_layer = (
    Jolab.RoughInterface(
        Medium.((1, 1.5)),
        topography_zero,
        ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
    ),
    Propagation(
        (
            ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
            ReferenceFrame((0, 0, 100E-9), (0, 0, 0.0)),
        ),
        Medium(1.5),
    ),
    Jolab.RoughInterface(
        Medium.((1.5, 2)),
        topography_zero,
        ReferenceFrame((0, 0, 100E-9), (0, 0, 0.0)),
    ),
    Propagation(
        (
            ReferenceFrame((0, 0, 100E-9), (0, 0, 0.0)),
            ReferenceFrame((0, 0, 600E-9), (0, 0, 0.0)),
        ),
        Medium(2.0),
    ),
    Jolab.RoughInterface(
        Medium.((2, 2.5)),
        topography_zero,
        ReferenceFrame((0, 0, 600E-9), (0, 0, 0.0)),
    ),
    Propagation(
        (
            ReferenceFrame((0, 0, 600E-9), (0, 0, 0.0)),
            ReferenceFrame((0, 0, 800E-9), (0, 0, 0.0)),
        ),
        Medium(2.5),
    ),
    Jolab.RoughInterface(
        Medium.((2.5, 1)),
        topography_zero,
        ReferenceFrame((0, 0, 800E-9), (0, 0, 0.0)),
    ),
)

(r2, t2) = Jolab.lightinteraction_recursivegridded(
    mls_layer,
    beam,
    rtol = 1E-20,
    printBool = false,
)
@test isapprox(r2, r, rtol = 1E-6)
@test isapprox(t2, t, rtol = 1E-6)

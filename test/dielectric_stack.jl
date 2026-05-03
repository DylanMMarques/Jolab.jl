using Jolab, StaticArrays, Test

## Test errors 

stack_test = DielectricStack(
    Medium.((@SVector [1, 1.5, 1])),
    (@SVector [100E-9]),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)

grid_pw = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, 0.0, 0.1, 1550E-9)
pw = Jolab.AngularSpectrum(
    Jolab.Forward,
    grid_pw,
    reshape([1.0 + 0.0im], 1, 1),
    Medium(2.0),
    ReferenceFrame((1, 0, 0), (0, 0, 0)),
)
@test_throws ArgumentError light_interaction(stack_test, pw)

grid_pw = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, 0.0, 0.1, 1550E-9)
pw = Jolab.AngularSpectrum(
    Jolab.Forward,
    grid_pw,
    reshape([1.0 + 0.0im], 1, 1),
    Medium(1.0),
    ReferenceFrame((2, 0, 0), (0, 0, 0)),
)
@test_throws ArgumentError light_interaction(stack_test, pw)

grid_pw = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, 0.0, 0.1, 1550E-9)
pw = Jolab.AngularSpectrum(
    Jolab.Forward,
    grid_pw,
    reshape([1.0 + 0.0im], 1, 1),
    Medium(1.0),
    ReferenceFrame((2, 0, 0), (0, 0.5, 0.1)),
)
@test_throws ArgumentError light_interaction(stack_test, pw)

grid_pw = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, 0.0, 0.1, 1550E-9)
pw = Jolab.AngularSpectrum(
    Jolab.Forward,
    grid_pw,
    reshape([1.0 + 0.0im], 1, 1),
    Medium(1.0),
    ReferenceFrame((1, 0, 0), (0, 0.5, 0.1)),
)
@test_throws ArgumentError light_interaction(stack_test, pw)
## Test values
function test_reflection_coeffiecient(nsx, nsy, λ, mls)
    grid_pw = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsy, λ)
    pw = Jolab.AngularSpectrum(
        Jolab.Forward,
        grid_pw,
        reshape([1.0 + 0.0im], 1, 1),
        first(mls.mat),
        ReferenceFrame((0, 0, 0), (0, 0, 0)),
    )
    (rpw, tpw) = light_interaction(mls, pw)
    (rpw.e[1], tpw.e[1])
end
ref = ReferenceFrame((0, 0, 0), (0, 0, 0))
stack_test = DielectricStack(Medium.((@SVector [1, 2, 1])), (@SVector [100E-9]), ref)
@test all(
    test_reflection_coeffiecient(0.0, 0.0, 1550E-9, stack_test) .≈ (
        -0.3801572124640627 + 0.28909310142234657im,
        0.5318174730142268 + 0.6993395798310894im,
    ),
)
@test all(
    test_reflection_coeffiecient(0.0, 0.1, 1550E-9, stack_test) .≈
    (-0.3817505426194731 + 0.2902356449628659im, 0.531090735792542 + 0.6985502300894781im),
)
@test all(
    test_reflection_coeffiecient(0.1, 0.0, 1550E-9, stack_test) .≈
    (-0.3817505426194731 + 0.2902356449628659im, 0.531090735792542 + 0.6985502300894781im),
)
@test all(
    test_reflection_coeffiecient(0.1, 0.1, 1550E-9, stack_test) .≈ (
        -0.38336673222871664 + 0.29138460541823613im,
        0.5303414130920262 + 0.6977556491387847im,
    ),
)

stack_test = DielectricStack(Medium.((@SVector [1, 2 + im, 1])), (@SVector [200E-9]), ref)
@test all(
    test_reflection_coeffiecient(0.0, 0.0, 1550E-9, stack_test) .≈ (
        -0.4726748478534631 - 0.2266907860346169im,
        0.06294761742913318 + 0.3843200058489115im,
    ),
)
@test all(
    test_reflection_coeffiecient(0.1, 0.0, 1550E-9, stack_test) .≈ (
        -0.4748859437971371 - 0.22622476562202212im,
        0.06387585643028817 + 0.3830637619377823im,
    ),
)
@test all(
    test_reflection_coeffiecient(0.0, 0.1, 1550E-9, stack_test) .≈ (
        -0.4748859437971371 - 0.22622476562202212im,
        0.06387585643028817 + 0.3830637619377823im,
    ),
)
@test all(
    test_reflection_coeffiecient(0.1, 0.1, 1550E-9, stack_test) .≈ (
        -0.4771130976569241 - 0.22574918641516817im,
        0.06479985900217344 + 0.38179405532679594im,
    ),
)

stack_test = DielectricStack(Medium.((@SVector [2, 1, 2])), (@SVector [200E-9]), ref)
@test all(
    test_reflection_coeffiecient(0.0, 0.0, 1550E-9, stack_test) .≈
    (0.3801572124640627 - 0.28909310142234657im, 0.5318174730142268 + 0.6993395798310894im),
)
@test all(
    test_reflection_coeffiecient(0.1, 0.0, 1550E-9, stack_test) .≈
    (0.38004014247437345 - 0.2907048932717478im, 0.5335016668393845 + 0.6974497305293892im),
)
@test all(
    test_reflection_coeffiecient(0.0, 0.1, 1550E-9, stack_test) .≈
    (0.38004014247437345 - 0.2907048932717478im, 0.5335016668393845 + 0.6974497305293892im),
)
@test all(
    test_reflection_coeffiecient(0.1, 0.1, 1550E-9, stack_test) .≈
    (0.3799169568043824 - 0.2923210356280086im, 0.5351813183829595 + 0.6955519207907793im),
)

## Water to air for total internal reflection testing
stack_test =
    DielectricStack(Medium.((@SVector [1.33, 1, 1])), (@SVector Float64[10E-6]), ref)

@test all(
    .≈(
        test_reflection_coeffiecient(0.0, 0.8 * 1.333, 1550E-9, stack_test),
        (
            0.6431032383925089 - 0.765779488344437im,
            4.949051151775065e-7 - 2.3065390964139887e-7im,
        ),
        atol = 1E-13,
    ),
)
@test all(
    .≈(
        test_reflection_coeffiecient(0.0, 0.8 * 1.333, 1550E-9, stack_test),
        (
            0.6431032383925089 - 0.765779488344437im,
            4.949051151775065e-7 - 2.3065390964139887e-7im,
        ),
        atol = 1E-13,
    ),
)
@test all(
    .≈(
        test_reflection_coeffiecient(0.8 * 1.333, 0.0, 1550E-9, stack_test),
        (
            0.6431032383925089 - 0.765779488344437im,
            4.949051151775065e-7 - 2.3065390964139887e-7im,
        ),
        atol = 1E-13,
    ),
)
@test all(
    .≈(
        test_reflection_coeffiecient(0.7 * 1.333, 0.0, 1550E-9, stack_test),
        (0.4498594097636552 + 0.0im, -0.6182242523375912 + 1.3114461795673662im),
        atol = 1E-12,
        rtol = 1E-6,
    ),
)
@test all(
    .≈(
        test_reflection_coeffiecient(0.0, 0.7 * 1.333, 1550E-9, stack_test),
        (0.4498594097636552 + 0.0im, -0.6182242523375912 + 1.3114461795673662im),
        atol = 1E-12,
        rtol = 1E-6,
    ),
)

## Test radial with cartesian beam
stack_test = DielectricStack(Medium.((@SVector [1, 2, 3])), (@SVector Float64[10E-6]), ref)


## Backward mode
stack_f =
    DielectricStack(Medium.((@SVector [1, 4, 2+im, 1])), (@SVector [500E-9, 200E-9]), ref)
stack_b =
    DielectricStack(Medium.((@SVector [1, 2+im, 4, 1])), (@SVector [200E-9, 500E-9]), ref)
function test_reflection_coeffiecient(stack_forward, stack_backward, nsx, nsy, λ)
    grid_pw_f = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsy, λ)
    pw_f = Jolab.AngularSpectrum(
        Jolab.Forward,
        grid_pw_f,
        reshape([1.0 + 0.0im], 1, 1),
        first(stack_forward.mat),
        ref,
    )
    grid_pw_b = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsy, λ)
    pw_b = Jolab.AngularSpectrum(
        Jolab.Backward,
        grid_pw_b,
        reshape([1.0 + 0.0im], 1, 1),
        first(stack_backward.mat),
        last(stack_backward.frames),
    )
    (rpw_b, tpw_b) = light_interaction(stack_backward, pw_b)
    (rpw_f, tpw_f) = light_interaction(stack_forward, pw_f)
    rpw_f.e[1] ≈ tpw_b.e[1] && tpw_f.e[1] ≈ rpw_b.e[1]
end
@test test_reflection_coeffiecient(stack_b, stack_f, 0.1, 0.2, 1500E-9)

function f_beam(λ)
    mls = DielectricStack(
        Medium.((@SVector [1, 1.5, 1])),
        (@SVector [100E-9]),
        ReferenceFrame((0, 0, 0), (0, 0, 1)),
    )
    nsx = range(-0.95, 0.95, length = 50)
    grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, λ)
    beam = Jolab.AngularSpectrum(
        Jolab.Forward,
        grid_ang,
        reshape(nsx .* nsx', (length(nsx), length(nsx))),
        Medium(1.0),
        ReferenceFrame((0, 0, 0), (0, 0, 1)),
    )
    (rbeam, tbeam) = light_interaction(mls, beam)
end
@test_opt f_beam(1550E-9)

function test_scatmat_f(mls)
    nsx = range(0, 0.95, length = 50)
    grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, 1500E-9)
    beam = Jolab.AngularSpectrum(
        Float64,
        Jolab.Forward,
        grid_ang,
        reshape(nsx .* nsx', (length(nsx), length(nsx))),
        Medium(1.0),
        ReferenceFrame((0, 0, 0), (0, 0, 0)),
    )
    mat = ScatteringMatrix(mls, beam)
    (mat_r, mat_t) = light_interaction(mat, beam)
    (aux_r, aux_t) = light_interaction(mls, beam)
    mat_r ≈ aux_r && mat_t ≈ aux_t
end
mls = DielectricStack(
    Medium.((@SVector [1, 1.5, 1])),
    (@SVector [100E-9]),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)
mirror = Mirror(
    (Medium(1.0), Medium(1.0)),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
    reflectivity = 0.9,
)
@test test_scatmat_f(mls)
@test_opt test_scatmat_f(mls)
@test test_scatmat_f(mirror)

function test_scatmat_b(mls)
    nsx = range(0, 0.95, length = 50)
    grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, 1550E-9)
    beam = Jolab.AngularSpectrum(
        Float64,
        Jolab.Backward,
        grid_ang,
        reshape(nsx .* nsx', (length(nsx), length(nsx))),
        Medium(1.0),
        last(mls.frames),
    )
    mat = ScatteringMatrix(mls, beam)
    (mat_r, mat_t) = mat * beam
    (aux_r, aux_t) = light_interaction(mls, beam)
    mat_r ≈ aux_r && mat_t ≈ aux_t
end
@test test_scatmat_b(mls)
@test test_scatmat_b(mirror)

nsr = range(0, 0.5, length = 50)
nθ = 0.0
λ = 1550E-9
grid_field = Jolab.CylindricalGrid(Jolab.NSR_NSθ_λ, nsr, nθ, λ)
field = Jolab.AngularSpectrum(
    Jolab.Forward,
    grid_field,
    Jolab.Gaussian(10E-6),
    Medium(1),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)
mls = DielectricStack(
    Medium.((@SVector [1, 1.5, 1])),
    (@SVector [100E-9]),
    ReferenceFrame((0, 0, 100E-9), (0, 0, 0)),
)
@test_throws ArgumentError light_interaction(mls, field)
@test_throws ArgumentError ScatteringMatrix(mls, field)

using Jolab, Test
using Enzyme
using FiniteDiff
import FiniteDiff: finite_difference_derivative

# Lens tests 
#
function field_after_lens(focal_len, nsx, nsy, ω, λ, medium)
    lens = Lens(
        focal_len,
        0.5,
        (Medium(1), Medium(1.0)),
        ReferenceFrame((0, 0, focal_len), (0, 0, 0)),
    )
    grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsy, λ)
    field = Jolab.AngularSpectrum(
        Jolab.Forward,
        grid_ang,
        Jolab.Gaussian(ω),
        medium,
        ReferenceFrame((0, 0, 0), (0, 0, 0)),
    )
    (rfield, tfield) = light_interaction(lens, field)
    return maximum(abs, tfield.e)
end

at(ω) = field_after_lens(focal_len, nsx, nsx, ω, 1500E-9, Medium(1))
@test autodiff(Enzyme.Forward, at, Duplicated, Duplicated(10E-6, 1.0))[1] ≈
      finite_difference_derivative(at, [10E-6], Val{:central}, Float64, relstep = 1E-9)[1]

at(f) = field_after_lens(f, nsx, nsx, 10E-6, 1500E-9, Medium(1))
@test autodiff(Enzyme.Forward, at, Duplicated, Duplicated(10E-6, 1.0))[1] ≈
      finite_difference_derivative(
    x -> at(x),
    [10E-6],
    Val{:central},
    Float64,
    relstep = 1E-9,
)[1]

function field_after_lens(focal_len, x, y, ω, λ, medium)
    lens = Lens(
        focal_len,
        0.5,
        (Medium(1), Medium(1.0)),
        ReferenceFrame((0, 0, focal_len), (0, 0, 0)),
    )
    grid_spatial = Jolab.CartesianGrid(Jolab.X_Y_λ, x, y, λ)
    field = Jolab.SpatialBeam(
        Jolab.Forward,
        grid_spatial,
        Jolab.Gaussian(ω),
        medium,
        ReferenceFrame((0, 0, 0), (0, 0, 0)),
    )
    (rfield, tfield) = light_interaction(lens, field)
    return maximum(abs, tfield.e)
end

at(ω) = field_after_lens(focal_len, x, x, ω, 1500E-9, Medium(1))
@test isapprox(
    autodiff(Enzyme.Forward, at, Duplicated, Duplicated(1E-6, 1.0))[1],
    finite_difference_derivative(
        x -> at(x),
        [1E-6],
        Val{:central},
        Float64,
        relstep = 1E-9,
    )[1],
    rtol = 1E-3,
)

at(focal_len) = field_after_lens(focal_len, x, x, 10E-6, 1500E-9, Medium(1))
@test autodiff(Enzyme.Forward, at, Duplicated, Duplicated(10E-6, 1.0))[1] ≈
      finite_difference_derivative(
    x -> at(x),
    [10E-6],
    Val{:central},
    Float64,
    relstep = 1E-9,
)[1]

# Test with Axicon
function test_autodiff(λ)
    nsr = range(0, 0.5, length = 100)
    axicon = Axicon(
        5*π/180,
        Medium(1.44),
        Medium(1),
        ReferenceFrame((0, 0, 0), (0, 0, 0));
        solver = (nsr = nsr,),
    )

    r = range(0, 2E-3, length = 100)
    nθ = 0.0
    grid_field = Jolab.CylindricalGrid(Jolab.R_θ_λ, r, nθ, λ)
    field = Jolab.SpatialBeam(
        Jolab.Forward,
        grid_field,
        Jolab.Gaussian(1E-3),
        Medium(1),
        ReferenceFrame((0, 0, 0), (0, 0, 0)),
    )
    (rfield, tfield) = light_interaction(axicon, field)
    intensity(tfield)
end

autodiff(Enzyme.Forward, test_autodiff, Duplicated, Duplicated(1500E-9, 1.0))

# ReferenceFrame tests

# using Enzyme, Jolab, StaticArrays, FiniteDiff
# import FiniteDiff: finite_difference_derivative
# Enzyme.API.runtimeActivity!(true)
# 
# function f(in_x)
#     nsx, nsy, e, λ, n, k, x, y, z = in_x
#     pw = PlaneWaveScalar(nsx, nsy, e + im, λ * 1E-9, Medium(n + k*im), ReferenceFrame((x,y,z), (0,0,0)))
#     pw2 = translate_referenceframe(pw, (2x, 2y, 2z))
#     pw2.e
# end
# 
# in_x = [0.1, 0.1, 1.0, 1550, 1, 0, 10E-9, 10E-9, 10E-9]
# enz = enz = enz = enz = Enzyme.jacobian(Forward, f, in_x)
# fin = FiniteDiff.finite_difference_jacobian(f, in_x, Val{:central}, ComplexF64, relstep = 1E-5)
# @test all(isapprox.(enz, fin, atol = 1E-10))
# 
# function f(in_x)
#     nsx, nsy, e, λ, n, x, y, z, θ, ϕ = in_x
#     pw = PlaneWaveScalar(nsx, nsy, e + 0im, λ * 1E-9, Medium(n), ReferenceFrame((x * 1E-9, y * 1E-9, z * 1E-9), (0,0,0)))
#     pw2 = translate_referenceframe(pw, (2x, 2y, 2z))
#     pw3 = rotate_referenceframe(pw2, (θ, 0, ϕ))
#     complex(pw3.nsx), complex(pw3.nsy), pw3.e
# end
# in_x = [0.1, 0.1, 1.0, 1550, 1, 10, 10, 10, π/4, π/8]
# f(in_x)
# enz = Enzyme.jacobian(Forward, f, in_x)
# fin_nsx = FiniteDiff.finite_difference_jacobian(i -> f(i)[1], in_x, Val{:central}, ComplexF64, relstep = 1E-5)
# fin_nsy = FiniteDiff.finite_difference_jacobian(i -> f(i)[2], in_x, Val{:central}, ComplexF64, relstep = 1E-5)
# fin_e = FiniteDiff.finite_difference_jacobian(i -> f(i)[3], in_x, Val{:central}, ComplexF64, relstep = 1E-12)
# @test all(isapprox.(first.(enz), fin_nsx, atol = 1E-10))
# @test all(isapprox.(getindex.(enz,2), fin_nsy, atol = 1E-10))
# @test all(isapprox.(last.(enz), fin_e, rtol = 1E-3))
#
# DielectricStack tests
## Jacobian tests
function jac(x)
    n, k, h, nsx, nsy, λ = x
    stack_test = DielectricStack(
        Float64,
        Medium.((@SVector [1 + k*im, n, 1])),
        (@SVector [h]),
        ReferenceFrame((0, 0, 0), (0, 0, 1)),
    )
    Jolab.rtss(stack_test, Forward, sqrt(nsx^2 + nsy^2), λ)
end
enz_jac = Enzyme.jacobian(Enzyme.Forward, jac, [2, 1, 100E-9, 0.1, 0.1, 1550E-9])

fin_jac_r = FiniteDiff.finite_difference_jacobian(
    x -> jac(x)[1],
    [2, 1, 100E-9, 0.1, 0.1, 1550E-9],
    Val{:central},
    ComplexF64,
    relstep = 1E-9,
)
fin_jac_t = FiniteDiff.finite_difference_jacobian(
    x -> jac(x)[2],
    [2, 1, 100E-9, 0.1, 0.1, 1550E-9],
    Val{:central},
    ComplexF64,
    relstep = 1E-9,
)

@test all(isapprox.(first.(enz_jac[1]), vec(fin_jac_r), rtol = 1E-4))
@test all(isapprox.(last.(enz_jac[1]), vec(fin_jac_t), rtol = 1E-4))

## Plane Wave tests
# Wavelength dependency

stack_test = DielectricStack(
    Medium.((@SVector [1+0im, 2.0, 1])),
    (@SVector [100E-9]),
    ReferenceFrame((0, 0, 0), (0, 0, 1)),
)

ad_diff = autodiff(
    Enzyme.Forward,
    Jolab.rtss,
    Duplicated,
    Const(stack_test),
    Const(Forward),
    Const(0.1),
    Duplicated(1550E-9, 1.0),
)[1]
num_diff = (
    finite_difference_derivative(
        i -> Jolab.rtss(stack_test, Forward, 0.1, i)[1],
        1550E-9;
        absstep = 1E-18,
    ),
    finite_difference_derivative(
        i -> Jolab.rtss(stack_test, Forward, 0.1, i)[2],
        1550E-9;
        absstep = 1E-18,
    ),
)
@test all((num_diff) .≈ ad_diff)

# direction dependency
ad_diff = autodiff(
    Enzyme.Forward,
    Jolab.rtss,
    Duplicated,
    Const(stack_test),
    Const(Forward),
    Duplicated(0.1, 1.0),
    Const(1550E-9),
)[1]
num_diff = (
    finite_difference_derivative(
        i -> Jolab.rtss(stack_test, Forward, i, 1550E-9)[1],
        0.1;
        absstep = 1E-18,
    ),
    finite_difference_derivative(
        i -> Jolab.rtss(stack_test, Forward, i, 1550E-9)[2],
        0.1;
        absstep = 1E-10,
    ),
)
@test all((num_diff) .≈ ad_diff)

# Refractive index dependency
function rtss_f(n)
    mls = DielectricStack(
        Medium.((@SVector [1+0im, n, 1])),
        (@SVector [100E-9]),
        ReferenceFrame((0, 0, 0), (0, 0, 1)),
    )
    Jolab.rtss(mls, Forward, 0.1, 1550E-9)
end
ad_diff = autodiff(Enzyme.Forward, rtss_f, Duplicated, Duplicated(2.0, 1.0))[1]
num_diff = (
    finite_difference_derivative(i -> rtss_f(i)[1], 2.0; absstep = 1E-18),
    finite_difference_derivative(i -> rtss_f(i)[2], 2.0; absstep = 1E-18),
)
@test all(ad_diff .≈ num_diff)

# Thickness dependency
function thickness_dependency(h)
    mls = DielectricStack(
        Medium.((@SVector [1+0im, 2, 1])),
        (@SVector [h]),
        ReferenceFrame((0, 0, 0), (0, 0, 1)),
    )
    Jolab.rtss(mls, Forward, 0.1, 1550E-9)
end
ad_diff =
    autodiff(Enzyme.Forward, thickness_dependency, Duplicated, Duplicated(10E-9, 1.0))[1]
num_diff = (
    finite_difference_derivative(i -> thickness_dependency(i)[1], 10E-9; absstep = 1E-20),
    finite_difference_derivative(i -> thickness_dependency(i)[2], 10E-9; absstep = 1E-20),
)
@test all(ad_diff .≈ num_diff)


function f_rbeam(mls, nsx, nsy, waist, λ, medium, frame)
    grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsy, λ)
    beam = Jolab.AngularSpectrum(
        Float64,
        Forward,
        grid_ang,
        Jolab.Gaussian(waist),
        medium,
        frame,
    )
    (rbeam, tbeam) = light_interaction(mls, beam)
    intensity(rbeam)
end

function f_tbeam(mls, nsx, nsy, waist, λ, medium, frame)
    grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsy, λ)
    beam = Jolab.AngularSpectrum(
        Float64,
        Forward,
        grid_ang,
        Jolab.Gaussian(waist),
        medium,
        frame,
    )
    (rbeam, tbeam) = light_interaction(mls, beam)
    intensity(tbeam)
end

medium = Medium(1.0)
frame = ReferenceFrame((0, 0, 0), (0, 0, 1))
nsx = range(-0.5, 0.5, length = 100)
nsy = nsx
waist = 10E-6
T = Float64
mls = DielectricStack(
    Medium.((@SVector [1, 1.5, 1])),
    (@SVector [100E-9]),
    ReferenceFrame((0, 0, 0), (0, 0, 1)),
)

ar(λ) = f_rbeam(mls, nsx, nsy, waist, λ * 1E-9, medium, frame)
at(λ) = f_tbeam(mls, nsx, nsy, waist, λ * 1E-9, medium, frame)

num_diff = finite_difference_derivative(ar, 1550.0; absstep = 1E-50)
val = ar(1550.0)
ad_diff = Tuple(
    autodiff(
        set_runtime_activity(Enzyme.ForwardWithPrimal),
        ar,
        Duplicated,
        Duplicated(1550.0, 1.0),
    ),
)
@test all((ad_diff) .≈ (num_diff, val))

num_diff = finite_difference_derivative(at, 1550.0; absstep = 1E-20)
val = at(1550.0)
ad_diff = Tuple(
    autodiff(
        set_runtime_activity(Enzyme.ForwardWithPrimal),
        at,
        Duplicated,
        Duplicated(1550.0, 1.0),
    ),
)
@test all((ad_diff) .≈ (num_diff, val))

# Broken tests
# mirror = Mirror((Medium(1.0), Medium(1.0)), ReferenceFrame((0,0,0), (0,0,0)); reflectivity = 0.99)
# ar(λ) = f_rbeam(mirror, nsx, nsy, e, medium, frame, λ)
# at(λ) = f_tbeam(mirror, nsx, nsy, e, medium, frame, λ)
#
# num_diff = finite_difference_derivative(i -> f_rbeam(mirror, nsx, nsy, e, medium, frame, i), 1550E-9; absstep = 1E-20)
# val = f_rbeam(mirror, nsx, nsy, e, medium, frame, 1550E-9)
# ad_diff = Tuple(autodiff(Enzyme.Forward, ar, Duplicated, Duplicated(1550E-9, 1.0)))
# @test all((ad_diff) .≈ (val, num_diff))
#
# num_diff = finite_difference_derivative(i -> f_tbeam(mirror, nsx, nsy, e, medium, frame, i), 1550E-9; absstep = 1E-20)
# val = f_tbeam(mirror, nsx, nsy, e, medium, frame, 1550E-9)
# ad_diff = Tuple(autodiff(Enzyme.Forward, at, Duplicated, Duplicated(1550E-9, 1.0)))
# @test all((ad_diff) .≈ (val, num_diff))
#
#
# ## Auto diff on Radial symmetric beams
# function f_beam(mls, nsr, ω, medium, frame, λ)
#     beam2 = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Float64, Forward, nsr, ω, λ, medium, frame)
#     (rbeam, tbeam) = light_interaction(mls, beam2)
#     (intensity(tbeam), intensity(rbeam))
# end
# aux_f(λ) = f_beam(mls, nsx, 10E-6, medium, frame, λ)
# ad_diff = autodiff(Enzyme.Forward, aux_f, Duplicated, Duplicated(1550E-9, 1.0))
# num_diff_t = finite_difference_derivative(i -> aux_f(i)[1], 1550E-9; absstep = 1E-20)
# num_diff_r = finite_difference_derivative(i -> aux_f(i)[2], 1550E-9; absstep = 1E-20)
# @test all(ad_diff[1] .≈ (num_diff_t, num_diff_r))

using Jolab, StaticArrays, Test
import Jolab: Forward, Backward

## Test errors 

stack_test = DielectricStack(Medium.((@SVector [1, 1.5, 1])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,0)))

pw = PlaneWaveScalar(Forward, 0, 0.1, 1, 1550E-9, Medium(2.0), ReferenceFrame((1,0,0), (0,0,0)))
@test_throws ArgumentError light_interaction(stack_test, pw)

pw = PlaneWaveScalar(Forward, 0, 0.1, 1, 1550E-9, Medium(1.0), ReferenceFrame((2,0,0), (0,0,0)))
@test_throws ArgumentError light_interaction(stack_test, pw)

pw = PlaneWaveScalar(Forward, 0, 0.1, 1, 1550E-9, Medium(1.0), ReferenceFrame((2,0,0), (0,.5,.1)))
@test_throws ArgumentError light_interaction(stack_test, pw)

pw = PlaneWaveScalar(Forward, 0, 0.1, 1, 1550E-9, Medium(1.0), ReferenceFrame((1,0,0), (0,.5,.1)))
@test_throws ArgumentError light_interaction(stack_test, pw)

## Test types 

function test_type(T)
    stack_test = DielectricStack(T, Medium.((@SVector [1, 1.5, 1])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,0)))
    pw = PlaneWaveScalar(T, Forward, 0, 0.1, 1, 1550E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))
    (rpw, tpw) = light_interaction(stack_test, pw)
    typeof(rpw.e) == Complex{T} && typeof(rpw.nsx) == T && typeof(tpw.e) == Complex{T} && typeof(tpw.nsx) == T
end
@test test_type(BigFloat)
@test test_type(Float64)
@test test_type(Float32)
@test test_type(Float16)

function test_type_complex(T)
    stack_test = DielectricStack(T, Medium.((@SVector [1 + im, 1.5 + im, 1 + im])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,0)))
    pw = PlaneWaveScalar(T, Forward, 0, 0.1, 1, 1550E-9, Medium(1 + im), ReferenceFrame((0,0,0), (0,0,0)))
    (rpw, tpw) = light_interaction(stack_test, pw)
    typeof(rpw.e) == Complex{T} && typeof(rpw.nsx) == T && typeof(tpw.e) == Complex{T} && typeof(tpw.nsx) == T
end
@test test_type_complex(BigFloat)
@test test_type_complex(Float64)
@test test_type_complex(Float32)
@test test_type_complex(Float16)

## Test values
function test_reflection_coeffiecient(nsx, nsy, λ, mls)
    pw = PlaneWaveScalar(Forward, nsx, nsy, 1, λ, first(mls.mat), ReferenceFrame((0,0,0), (0,0,0)))
    (rpw, tpw) = light_interaction(mls, pw)
    (rpw.e, tpw.e)
end
ref = ReferenceFrame((0,0,0), (0,0,0))
stack_test = DielectricStack(Medium.((@SVector [1, 2, 1])), (@SVector [100E-9]), ref)
@test all(test_reflection_coeffiecient(0, 0, 1550E-9, stack_test) .≈ (-0.3801572124640627 + 0.28909310142234657im, 0.5318174730142268 + 0.6993395798310894im))
@test all(test_reflection_coeffiecient(0, 0.1, 1550E-9, stack_test) .≈ (-0.3817505426194731 + 0.2902356449628659im, 0.531090735792542 + 0.6985502300894781im))
@test all(test_reflection_coeffiecient(0.1, 0, 1550E-9, stack_test) .≈ (-0.3817505426194731 + 0.2902356449628659im, 0.531090735792542 + 0.6985502300894781im))
@test all(test_reflection_coeffiecient(0.1, 0.1, 1550E-9, stack_test) .≈ (-0.38336673222871664 + 0.29138460541823613im, 0.5303414130920262 + 0.6977556491387847im))

stack_test = DielectricStack(Medium.((@SVector [1, 2 + im, 1])), (@SVector [200E-9]), ref)
@test all(test_reflection_coeffiecient(0, 0, 1550E-9, stack_test) .≈ (-0.4726748478534631 - 0.2266907860346169im, 0.06294761742913318 + 0.3843200058489115im))
@test all(test_reflection_coeffiecient(0.1, 0, 1550E-9, stack_test) .≈ (-0.4748859437971371 - 0.22622476562202212im, 0.06387585643028817 + 0.3830637619377823im))
@test all(test_reflection_coeffiecient(0, 0.1, 1550E-9, stack_test) .≈ (-0.4748859437971371 - 0.22622476562202212im, 0.06387585643028817 + 0.3830637619377823im))
@test all(test_reflection_coeffiecient(0.1, 0.1, 1550E-9, stack_test) .≈ (-0.4771130976569241 - 0.22574918641516817im, 0.06479985900217344 + 0.38179405532679594im))

stack_test = DielectricStack(Medium.((@SVector [2, 1, 2])), (@SVector [200E-9]), ref)
@test all(test_reflection_coeffiecient(0, 0, 1550E-9, stack_test) .≈ (0.3801572124640627 - 0.28909310142234657im, 0.5318174730142268 + 0.6993395798310894im)
)
@test all(test_reflection_coeffiecient(0.1, 0, 1550E-9, stack_test) .≈ (0.38004014247437345 - 0.2907048932717478im, 0.5335016668393845 + 0.6974497305293892im)
)
@test all(test_reflection_coeffiecient(0, 0.1, 1550E-9, stack_test) .≈ (0.38004014247437345 - 0.2907048932717478im, 0.5335016668393845 + 0.6974497305293892im)
)
@test all(test_reflection_coeffiecient(0.1, 0.1, 1550E-9, stack_test) .≈ (0.3799169568043824 - 0.2923210356280086im, 0.5351813183829595 + 0.6955519207907793im)
)

## Water to air for total internal reflection testing
stack_test = DielectricStack(Medium.((@SVector [1.33, 1, 1])), (@SVector Float64[10E-6]), ref)
@test all(test_reflection_coeffiecient(0, 0.8 * 1.333, 1550E-9, stack_test) .≈ (0.6431032383925089 - 0.765779488344437im, 4.949051151775065e-7 - 2.3065390964139887e-7im))
@test all(test_reflection_coeffiecient(0.8 * 1.333, 0, 1550E-9, stack_test) .≈ (0.6431032383925089 - 0.765779488344437im, 4.949051151775065e-7 - 2.3065390964139887e-7im))
@test all(test_reflection_coeffiecient(0.7 * 1.333, 0, 1550E-9, stack_test) .≈ (0.4498594097636552 + 0.0im, -0.6182242523375912 + 1.3114461795673662im))
@test all(test_reflection_coeffiecient(0, 0.7 * 1.333, 1550E-9, stack_test) .≈ (0.4498594097636552 + 0.0im, -0.6182242523375912 + 1.3114461795673662im))

## Backward mode
stack_f = DielectricStack(Medium.((@SVector [1, 4, 2+im, 1])), (@SVector [500E-9, 200E-9]), ref)
stack_b = DielectricStack(Medium.((@SVector [1, 2+im, 4, 1])), (@SVector [200E-9, 500E-9]), ref)
function test_reflection_coeffiecient(stack_forward, stack_backward, nsx, nsy, λ)
    pw_f = PlaneWaveScalar(Forward, nsx, nsy, 1, λ, first(stack_forward.mat), ref)
    pw_b = PlaneWaveScalar(Backward, nsx, nsy, 1, λ, first(stack_backward.mat), last(stack_backward.frames))
    (rpw_b, tpw_b) = light_interaction(stack_backward, pw_b)
    (rpw_f, tpw_f) = light_interaction(stack_forward, pw_f)
    rpw_f.e ≈ tpw_b.e && tpw_f.e ≈ rpw_b.e
end
@test test_reflection_coeffiecient(stack_b, stack_f, 0.1, 0.2, 1500E-9)

function f_beam(λ)
    mls = DielectricStack(Medium.((@SVector [1, 1.5, 1])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,1)))
    nsx = range(0, 0.95, length = 100) .+ zeros(100)'
    beam = Jolab.monochromatic_angularspectrum(Float64, Forward, nsx, nsx', (nsx .* nsx)', λ, Medium(1.0), ReferenceFrame((0,0,0), (0,0,1)));
    (rbeam, tbeam) = light_interaction(mls, beam)
end
f_beam(1550E-9);

function test_scatmat_f()
    mls = DielectricStack(Medium.((@SVector [1, 1.5, 1])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,0)))
    nsx = range(0, 0.95, length = 100) .+ zeros(100)'
    beam = Jolab.monochromatic_angularspectrum(Float64, Forward, nsx, nsx', (nsx .* nsx)', 1500E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)));
    mat = Jolab.ScatteringMatrix(mls, beam)
    (mat_r, mat_t) = light_interaction(mat, beam)
    (aux_r, aux_t) = light_interaction(mls, beam)
    mat_r ≈ aux_r && mat_t ≈ aux_t
end
@test test_scatmat_f()

function test_scatmat_b()
    mls = DielectricStack(Medium.((@SVector [1, 1.5, 1])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,0)))
    nsx = range(0, 0.95, length = 100) .+ zeros(100)'
    beam = Jolab.monochromatic_angularspectrum(Float64, Backward, nsx, nsx', (nsx .* nsx)', 1550E-9, Medium(1.0), last(mls.frames));
    mat = Jolab.ScatteringMatrix(mls, beam)
    (mat_r, mat_t) = light_interaction(mat, beam)
    (aux_r, aux_t) = light_interaction(mls, beam)
    mat_r ≈ aux_r && mat_t ≈ aux_t
end
@test test_scatmat_b()

## Enzyme tests on Dielectric Stack
using Enzyme, Jolab, StaticArrays, FiniteDiff
import FiniteDiff: finite_difference_derivative
Enzyme.API.runtimeActivity!(true)

## Jacobian tests
function jac(x)
    n, k, h, nsx, nsy, λ = x
    stack_test = DielectricStack(Float64, Medium.((@SVector [1 + k*im, n, 1])), (@SVector [h]), ReferenceFrame((0,0,0), (0,0,1)))
    Jolab.rtss(stack_test, Forward, nsx, nsy, λ)
end
enz_jac = Enzyme.jacobian(Enzyme.Forward, jac, [2, 1, 100E-9, .1, .1, 1550E-9])

fin_jac_r = FiniteDiff.finite_difference_jacobian(x -> jac(x)[1], [2, 1, 100E-9, .1, .1, 1550E-9], Val{:central}, ComplexF64, relstep = 1E-9)
fin_jac_t = FiniteDiff.finite_difference_jacobian(x -> jac(x)[2], [2, 1, 100E-9, .1, .1, 1550E-9], Val{:central}, ComplexF64, relstep = 1E-9)

@test all(isapprox.(first.(enz_jac), fin_jac_r, rtol = 1E-4))
@test all(isapprox.(last.(enz_jac), fin_jac_t, rtol = 1E-4))

## Plane Wave tests
# Wavelength dependency

stack_test = DielectricStack(Medium.((@SVector [1+0im, 2.0, 1])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,1)))

ad_diff = autodiff(Enzyme.Forward, Jolab.rtss, Duplicated, Const(stack_test), Const(Forward), Const(0.1), Const(0.1), Duplicated(1550E-9, 1.0))[2]
num_diff = (finite_difference_derivative(i -> Jolab.rtss(stack_test, Forward, 0.1, 0.1, i)[1], 1550E-9; absstep = 1E-18), finite_difference_derivative(i -> Jolab.rtss(stack_test, Forward, 0.1, 0.1, i)[2], 1550E-9; absstep = 1E-18))
@test all((num_diff) .≈ ad_diff)

# direction dependency
ad_diff = autodiff(Enzyme.Forward, Jolab.rtss, Duplicated, Const(stack_test), Const(Forward), Duplicated(0.1, 1.0), Const(0.1), Const(1550E-9))[2]
num_diff = (finite_difference_derivative(i -> Jolab.rtss(stack_test, Forward, i, 0.1, 1550E-9)[1], 0.1; absstep = 1E-18),
    finite_difference_derivative(i -> Jolab.rtss(stack_test, Forward, i, 0.1, 1550E-9)[2], 0.1; absstep = 1E-10))
@test all((num_diff) .≈ ad_diff)

# Refractive index dependency
function f(n)
    mls = DielectricStack(Medium.((@SVector [1+0im, n, 1])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,1)))
    Jolab.rtss(mls, Forward, 0.1, 0.1, 1550E-9)
end
ad_diff = autodiff(Enzyme.Forward, f, Duplicated, Duplicated(2.0, 1.0))[2]
num_diff = (finite_difference_derivative(i -> f(i)[1], 2.0; absstep = 1E-18),
    finite_difference_derivative(i -> f(i)[2], 2.0; absstep = 1E-18))
@test all(ad_diff .≈ num_diff)

# Thickness dependency
function f(h)
    mls = DielectricStack(Medium.((@SVector [1+0im, 2, 1])), (@SVector [h]), ReferenceFrame((0,0,0), (0,0,1)))
    Jolab.rtss(mls, Forward, 0.1, 0.1, 1550E-9)
end
ad_diff = autodiff(Enzyme.Forward, f, Duplicated, Duplicated(10E-9, 1.0))[2]
num_diff = (finite_difference_derivative(i -> f(i)[1], 10E-9; absstep = 1E-20),
    finite_difference_derivative(i -> f(i)[2], 10E-9; absstep = 1E-20))
@test all(ad_diff .≈ num_diff)


function f_rbeam(mls, nsx, nsy, e, medium, frame, λ)
    beam2 = Jolab.monochromatic_angularspectrum(Float64, Forward, nsx, nsy, e, λ, medium, frame)
    (rbeam, tbeam) = light_interaction(mls, beam2)
    intensity(rbeam)
end

function f_tbeam(mls, nsx, nsy, e, medium, frame, λ)
    beam2 = Jolab.monochromatic_angularspectrum(Float64, Forward, nsx, nsy, e, λ, medium, frame)
    (rbeam, tbeam) = light_interaction(mls, beam2)
    intensity(tbeam)
end

medium = Medium(1.0)
frame = ReferenceFrame((0,0,0), (0,0,1))
nsx = range(0, 0.5, length = 100) .+ zeros(100)'
nsy = collect(nsx')
e = deepcopy(nsx)
T = Float64
mls = DielectricStack(Medium.((@SVector [1, 1.5, 1])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,1)))

ar(λ) = f_rbeam(mls, nsx, nsy, e, medium, frame, λ)
at(λ) = f_tbeam(mls, nsx, nsy, e, medium, frame, λ)

num_diff = finite_difference_derivative(i -> f_rbeam(mls, nsx, nsy, e, medium, frame, i), 1550E-9; absstep = 1E-20)
val = f_rbeam(mls, nsx, nsy, e, medium, frame, 1550E-9)
ad_diff = Tuple(autodiff(Enzyme.Forward, ar, Duplicated, Duplicated(1550E-9, 1.0)))
@test all((ad_diff) .≈ (val, num_diff))

num_diff = finite_difference_derivative(i -> f_tbeam(mls, nsx, nsy, e, medium, frame, i), 1550E-9; absstep = 1E-20)
val = f_tbeam(mls, nsx, nsy, e, medium, frame, 1550E-9)
ad_diff = Tuple(autodiff(Enzyme.Forward, at, Duplicated, Duplicated(1550E-9, 1.0)))
@test all((ad_diff) .≈ (val, num_diff))
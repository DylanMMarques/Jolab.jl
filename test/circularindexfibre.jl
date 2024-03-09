import Pkg
Pkg.activate("/home/dylan/.julia/dev/Jolab/test")
using Jolab, Test, Enzyme, Optim
import Jolab: Forward, Backward
Enzyme.API.runtimeActivity!(true)

profile = CircularStepIndexProfile(10E-6, 0.2, Medium(1.55))

fibre = Fibre(profile, 1, Medium.((1,1)), (ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0,0,1), (0,0,0))))

Jolab.findmodes!(fibre, 1500E-9)

x = range(-100E-6, 100E-6, length = 20)

field = MonochromaticSpatialBeam(Float64, Forward, x, x, x .* x', 1500E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))

(back, forw) = light_interaction(fibre, field)

na = 0.05
λ = 1.5e-6

V(r, λ, na) = 2π * r / λ * na
r(na, λ, V) = V / 2π * λ / na

nclad = Jolab.nclad(1.5, na)
k = 2π / λ
b_cte(β) = (β - k * nclad) / k / (profile.ncore.n - nclad)

r_s = 5E-6:5E-6:200E-6
V_s = V.(r_s, λ, na)
n_modes(r_i) = begin 
    profile = CircularStepIndexProfile(r_i, na, Medium(1.5))
    fibre = Fibre(profile, 1, Medium.((1,1)), (ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0,0,1), (0,0,0))))
    Jolab.findmodes!(fibre, λ)
    modes = fibre.modes[λ]
    (modes, length(modes))
end
val = n_modes.(r_s)

using Enzyme, Jolab, StaticArrays, FiniteDiff
import FiniteDiff: finite_difference_derivative

function coupling_field_derivative(x, fibre, e, λ, medium,ref)
    e_m = zeros(ComplexF64, length(x), length(x))
    e_m[1] = e
    field = MonochromaticSpatialBeam(Forward, x, x, e_m, λ, medium, ref)
    (back, forw) = light_interaction(fibre, field)
    forw.modes.e
end
x = range(-150E-6, 150E-6, length = 10)
const x2 = x
e = ones(length(x), length(x)) .+ eps()
b = ones(size(e))
const fibre2 = fibre 

frame = ReferenceFrame((0,0,0), (0,0,0))
e = 1.0
tmp_(e, λ) = coupling_field_derivative(x2, fibre, e, λ, Medium(1.0), frame)

enz = autodiff(Enzyme.Forward, tmp_, Duplicated, Duplicated(e, deepcopy(e)), Const(1500E-9))
fd_val = finite_difference_derivative(i -> tmp_(i, 1500E-9), 1.0)
@test all(isapprox.(fd_val, enz[2]))

enz = autodiff(Enzyme.Forward, tmp_, Duplicated, Const(e), Duplicated(1500E-9, 1.0))

function coupling_field_derivative(x, fibre, e, λ, medium,ref)
    e_2 = reshape(e, length(x), length(x))
    field = MonochromaticSpatialBeam(Forward, x, x, e_2, λ, medium, ref)
    (back, forw) = light_interaction(fibre, field)
    forw.modes.e
end
tmp_2(e) = coupling_field_derivative(x2, fibre, e, 1500E-9, Medium(1.0), frame)
e = ones(length(x), length(x)) .+ eps()
tmp_2(e)
enz_2 = autodiff(Enzyme.Forward, tmp_2, DuplicatedNoNeed, Duplicated(vec(e), ones(length(e))))
enz = jacobian(Enzyme.Reverse, tmp_2, vec(e), Val(17))
all(isapprox.(enz * vec(e), enz_2.var"1"; atol = 1E-5))

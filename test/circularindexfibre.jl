import Pkg
Pkg.activate("/home/dylan/JolabAaTest/")
using Jolab, Test, Enzyme, Optim
import Jolab: Forward, Backward

profile = CircularStepIndexProfile(100E-6, 0.2, Medium(1.55))

fibre = Fibre(profile, 1, Medium.((1,1)), (ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0,0,1), (0,0,0))))

Jolab.findmodes!(fibre, 1500E-9)

x = -100E-6:1E-6:100E-6

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


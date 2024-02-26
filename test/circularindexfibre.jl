import Pkg
Pkg.activate("/home/dylan/JolabADTest/")
using Jolab, Test, Enzyme, Optim
import Jolab: Forward, Backward

profile = CircularStepIndexProfile(100E-6, 0.2, Medium(1.55))

fibre = Fibre(profile, 1, Medium.((1,1)), (ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0,0,1), (0,0,0))))

Jolab.findmodes!(fibre, 1500E-9)

_x = -100E-6:1E-6:100E-6
x = _x .+ 0(_x')
y = x'

field = Jolab.monochromatic_fieldspace(Float64, Forward, x, y, ones(size(x)), 1500E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))

(back, forw) = light_interaction(fibre, field)


na = 0.1
λ = 1.5e-6

V(r, λ, na) = 2π * r / λ * na
r(na, λ, V) = V / 2π * λ / na

nclad = Jolab.nclad(1.5, 1)
k = 2π / λ
b_cte(n3) = (k*n3 - k * nclad) / k / (profile.ncore.n - nclad)

r_s = 1E-6:1E-6:35E-6
V_s = V.(r_s, λ, na)
n_modes(r_i) = begin 
    profile = CircularStepIndexProfile(r_i, na, Medium(1.5))
    fibre = Fibre(profile, 1, Medium.((1,1)), (ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0,0,1), (0,0,0))))
    Jolab.findmodes!(fibre, λ)
    modes = fibre.modes[λ]
    length(modes)
end

for r_i in r_s
end

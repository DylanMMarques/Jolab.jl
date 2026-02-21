using Jolab, Test

z(x,y) = 0.0
int = Jolab.RoughInterface(Medium.((1, 2)), z, ReferenceFrame((0,0,0.0), (0,0,0.0)))
p_int = Jolab.DielectricStack(Medium.([1, 2]), zeros(0), ReferenceFrame((0,0,0.0), (0,0,0.0)))

nsx = range(-0.1, 0.1, length=64)
nsy = range(-0.1, 0.1, length=64)
λ = 1550E-9
beam = MonochromaticAngularSpectrum(Float64, Forward, nsx, nsx, ones(ComplexF64, 64, 64),λ, Medium(1.0), int.frame);

(r, t) = light_interaction(int, beam)
(r_p, t_p) = light_interaction(p_int, beam)
@test r ≈ r_p
@test t ≈ t_p

z_dist = 100E-9
z(x,y) = z_dist
int = Jolab.RoughInterface(Medium.((1, 2)), z, ReferenceFrame((0,0,0.0), (0,0,0.0)))
p_int = Jolab.DielectricStack(Medium.([1, 2]), zeros(0), ReferenceFrame((0,0,0.0), (0,0,0.0)))
nsz2 = sqrt.(complex.(2^2 .- nsx.^2 .- nsy'.^2));
nsz1 = sqrt.(complex.(1^2 .- nsx.^2 .- nsy'.^2));
(r, t) = light_interaction(int, beam)
intensity.((r, t))
(r_p, t_p) = light_interaction(p_int, beam)


## Compairison agaisnt arctan perturbation on a filter
data = [0.00408221220653889,
    0.0011402066408047376,
    0.000684671687514954,
    0.0026399636502819055,
    0.006679540081315143,
    0.01230492617854357]

ref1 = ReferenceFrame((0,0,0.0),(0,0., 0))
b = 1500E-9 / 2
ref2 = ReferenceFrame((0,0,b/1.5),(0,0.,0))
ref3 = ReferenceFrame((0,0,ref2.origin[3] + b),(0,0.,0))
ref4 = ReferenceFrame((0,0,ref3.origin[3] + b/1.5),(0,0.,0))
sx = range(-.7, .7, length = 64+1)
sx = sx[2:end]
λ = range(1450E-9, 1550E-9, length = 6)
zp(x,y) = 5E-9 * (atan(-1E7 * x) / π + .5)
zm(x,y) = -5E-9 * (atan(-1E7 * x) / π + .5)

RoughInterface(n1,n2,z, ref) = RoughMultilayerStructure([n1, n2], zeros(0), [z], ref)

rmls = [Jolab.RoughInterface(Medium.((1, 1.5)), zp, ref1), Jolab.RoughInterface(Medium.((1.5, 1)), zm, ref2), Jolab.RoughInterface(Medium.((1,1.5)),zp,ref3), Jolab.RoughInterface(Medium.((1.5,1)),zm,ref4)]
int = zeros(length(λ))

for i in 1:length(λ)
    local field
    local fieldr
    local fieldt
    field = MonochromaticAngularSpectrum_gaussian(Forward, nsx, nsy, 10E-6, 1550E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
    (fieldr, fieldt) = lightinteraction_recursivegridded(rmls, field, rtol = 1E-9)
    int[i] = intensity(fieldr)
    (fieldr, fieldt) = lightinteraction(rmls_r, field)
    int2[i] = intensity(fieldr)
end

@test all(isapprox.(int, data, rtol = 1E-8))
@test all(isapprox.(int2, data, rtol = 5E-2))
return true

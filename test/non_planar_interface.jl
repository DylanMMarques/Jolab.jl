using Jolab, Test, Unitful

zero_offset(x, y) = 0.0
int = Jolab.RoughInterface(
    Medium.((1, 2)),
    zero_offset,
    ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
)
p_int = Jolab.DielectricStack(
    Medium.([1, 2]),
    zeros(0),
    ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
)

@test dimension.(Jolab.rr_1(1.0, 1.0, 1.0, 1.0, 1550 * u"nm")) ==
      (NoDims, NoDims, dimension(u"nm^-1"), NoDims, dimension(u"nm^-1"), NoDims)

nsx = range(-0.5, 0.5, length = 64)
nsy = range(-0.5, 0.5, length = 64)
λ = 1550E-9
grid_beam = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, λ)
beam = Jolab.AngularSpectrum(
    Float64,
    Jolab.Forward,
    grid_beam,
    ones(ComplexF64, 64, 64),
    Medium(1.0),
    int.frame,
)

(r, t) = light_interaction(int, deepcopy(beam))
@test (@allocated light_interaction(int, deepcopy(beam))) < 1E7
(r_p, t_p) = light_interaction(p_int, deepcopy(beam))
@test isapprox(r, r_p, rtol = 1E-10)
@test isapprox(t, t_p, rtol = 1E-10)

grid_beam_b = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, λ)
beam = Jolab.AngularSpectrum(
    Float64,
    Jolab.Backward,
    grid_beam_b,
    ones(ComplexF64, 64, 64),
    Medium(2.0),
    int.frame,
)
(r, t) = light_interaction(int, deepcopy(beam))
(r_p, t_p) = light_interaction(p_int, deepcopy(beam))
@test isapprox(r, r_p, rtol = 1E-10)
@test isapprox(t, t_p, rtol = 1E-10)

grid_beam_f = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, λ)
beam = Jolab.AngularSpectrum(
    Float64,
    Jolab.Forward,
    grid_beam_f,
    ones(ComplexF64, 64, 64),
    Medium(1.0),
    int.frame,
)
z_dist = 1E-9
topography_dist(x, y) = z_dist
int = Jolab.RoughInterface(
    Medium.((1, 2)),
    topography_dist,
    ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
)
p_int = Jolab.DielectricStack(
    Medium.([1, 1, 2]),
    [z_dist],
    ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
)
(r, t) = light_interaction(int, beam)
(r_p, t_p) = light_interaction(p_int, beam)
(r_p, _) = translate_referenceframe(r_p, (0, 0, 0.0))
(_, t_p) = translate_referenceframe(t_p, (0, 0, 0.0))
@test isapprox(r, r_p, rtol = 1E-4)
@test isapprox(t, t_p, rtol = 1E-4)

## Compairison agaisnt arctan perturbation on a filter
data = [
    0.00408221220653889,
    0.0011402066408047376,
    0.000684671687514954,
    0.0026399636502819055,
    0.006679540081315143,
    0.01230492617854357,
]

ref1 = ReferenceFrame((0, 0, 0.0), (0, 0.0, 0))
b = 1500E-9 / 2
ref2 = ReferenceFrame((0, 0, b/1.5), (0, 0.0, 0))
ref3 = ReferenceFrame((0, 0, ref2.origin[3] + b), (0, 0.0, 0))
ref4 = ReferenceFrame((0, 0, ref3.origin[3] + b/1.5), (0, 0.0, 0))
sx = range(-0.7, 0.7, length = 64+1)[2:65]
λ = range(1450E-9, 1550E-9, length = 6)
topography_tan(x, y) = 5E-9 * (atan(-1E7 * x) / π + 0.5)
topography_tan_minus(x, y) = -topography_tan(x, y)


rmls = [
    Jolab.RoughInterface(Medium.((1, 1.5)), topography_tan, ref1),
    Propagation((ref1, ref2), Medium(1.5)),
    Jolab.RoughInterface(Medium.((1.5, 1)), topography_tan_minus, ref2),
    Propagation((ref2, ref3), Medium(1.0)),
    Jolab.RoughInterface(Medium.((1, 1.5)), topography_tan, ref3),
    Propagation((ref3, ref4), Medium(1.5)),
    Jolab.RoughInterface(Medium.((1.5, 1)), topography_tan_minus, ref4),
]
int = zeros(length(λ))

λ_base = 1500E-9
grid_ang_base = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, sx, sx, λ_base)
field = Jolab.AngularSpectrum(
    Jolab.Forward,
    grid_ang_base,
    Jolab.Gaussian(10E-6),
    Medium(1.0),
    ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
)
@time (fieldr, fieldt) =
    Jolab.lightinteraction_recursivegridded(rmls, field, printBool = false);

for i = 1:length(λ)
    grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, sx, sx, λ[i])
    field = Jolab.AngularSpectrum(
        Jolab.Forward,
        grid_ang,
        Jolab.Gaussian(10E-6),
        Medium(1.0),
        ReferenceFrame((0, 0, 0.0), (0, 0, 0.0)),
    )
    (fieldr, fieldt) =
        Jolab.lightinteraction_recursivegridded(rmls, field, rtol = 1E-9; printBool = false)
    int[i] = intensity(fieldr) / intensity(field)
end

@test all(isapprox.(int, data, rtol = 1E-6))

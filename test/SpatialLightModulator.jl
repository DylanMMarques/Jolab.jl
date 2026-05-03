using Jolab, Test
import Jolab: Forward

aperture = Jolab.SpatialLightModulator_aperture(
    Float64,
    1E-1,
    Medium.((1, 1)),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)
disk = Jolab.SpatialLightModulator_diskaperture(
    Float64,
    0.0,
    1E-1,
    Medium.((1, 1)),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)

nsx = range(-0.1, 0.1, length = 64)
nsy = range(-0.1, 0.1, length = 64)
λ = 1550E-9
grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsy, λ)
angspe = Jolab.AngularSpectrum(
    Float64,
    Forward,
    grid_ang,
    Jolab.Gaussian(10E-2),
    Medium(1),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)
@test_throws MethodError light_interaction(disk, angspe)

x = range(-0.1, 0.1, length = 64)
y = range(-0.1, 0.1, length = 64)
grid_spatial = Jolab.CartesianGrid(Jolab.X_Y_λ, x, y, λ)
beam = Jolab.SpatialBeam(
    Float64,
    Forward,
    grid_spatial,
    Jolab.Gaussian(10E-2),
    Medium(1),
    ReferenceFrame((0, 0, 0), (0, 0, 0)),
)

(_, t1) = light_interaction(aperture, beam)
(_, t2) = light_interaction(disk, beam)
@test t1 ≈ t2

@test_opt light_interaction(disk, beam)
@test_opt light_interaction(aperture, beam)

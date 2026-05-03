using Jolab, Test

mfd = 10E-6
frames = (ReferenceFrame((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), ReferenceFrame((0,0,1E-2), (0,0,0.0)))
media = Medium.((1,1))
fibre = SingleModeFibre(mfd, 1E-2, media, frames)
Jolab.findmodes!(fibre, 1500E-9)

x = LinRange(-50E-6, 50E-6, 100)
λ = 1500E-9
grid_spatial = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ)
spa = Jolab.SpatialBeam(Jolab.Forward, grid_spatial, Jolab.Gaussian(mfd), Medium(1.0), frames[1])
(fieldr, fieldt) = light_interaction(fibre, spa)
@test intensity(fieldt) ≈ 1 rtol = 1E-3

r = range(0, 50E-6, 1000)
t_sym = range(0, step = 2π, length = 1)
grid_radial = Jolab.CylindricalGrid(Jolab.R_θ_λ, r, t_sym, λ)
spa = Jolab.SpatialBeam(Jolab.Forward, grid_radial, Jolab.Gaussian(mfd), Medium(1.0), frames[1])
(fieldr, fieldt) = light_interaction(fibre, spa)
@test iszero(intensity(fieldr)) 
@test intensity(fieldt) ≈ 1 rtol = 1E-3

nsx = LinRange(-0.5, 0.5, 100)
grid_ang = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, λ)
field = Jolab.AngularSpectrum(Jolab.Forward, grid_ang, Jolab.Gaussian(mfd), Medium(1.0), frames[1])
(fieldr, fieldt) = light_interaction(fibre, field)
@test iszero(intensity(fieldr))
@test intensity(fieldt) ≈ 1 rtol = 1E-3

nsr = LinRange(0, 0.5, 1000)
grid_ang_rad = Jolab.CylindricalGrid(Jolab.NSR_NSθ_λ, nsr, t_sym, λ)
field = Jolab.AngularSpectrum(Jolab.Forward, grid_ang_rad, Jolab.Gaussian(mfd), Medium(1.0), frames[1])
@test_opt broken = (VERSION < v"1.11") light_interaction(fibre, field)
(fieldr, fieldt) = light_interaction(fibre, field)
@test intensity(fieldt) ≈ 1 rtol = 1E-3

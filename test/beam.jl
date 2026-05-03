using Jolab, Test
import Test: @inferred

# Test AngularSpectrum with NSX_NSY_λ coordinates
nsx = range(-.5, 0.5, length = 100)
λ = 1550E-9
grid_angular = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, λ)

beam = Jolab.AngularSpectrum(Jolab.Forward, grid_angular, (nsx .* nsx'), Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))

## Check same definition
@test isapprox(beam, beam)
@test (beam ≈ deepcopy(beam))

# Test with different nsx range
grid_angular_2 = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, 2nsx, 2nsx, λ)
beam2 = Jolab.AngularSpectrum(Jolab.Forward, grid_angular_2, (nsx .* nsx'), Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))
@test !isapprox(beam, beam2)

# Test with different wavelength
λ_2 = 1500E-9
grid_angular_3 = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, 2nsx, 2nsx, λ_2)
beam2 = Jolab.AngularSpectrum(Jolab.Forward, grid_angular_3, (nsx .* nsx'), Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))
@test !isapprox(beam, beam2)

# Test with different medium
beam2 = Jolab.AngularSpectrum(Jolab.Forward, grid_angular_3, (nsx .* nsx'), Medium(2.0), ReferenceFrame((0,0,0), (0,0,0)))
@test !isapprox(beam, beam2)

# Test with different reference frame
beam2 = Jolab.AngularSpectrum(Jolab.Forward, grid_angular_3, (nsx .* nsx'), Medium(2.0), ReferenceFrame((0,0,0), (0,0,1)))
@test !isapprox(beam, beam2)

# Test type-stability with random data
ns = range(-.1, 0.1, length = 100)
grid_ns = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, ns, ns, 1550E-9)
angspe = Jolab.AngularSpectrum(Jolab.Forward, grid_ns, rand(length(ns), length(ns)), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test_opt Jolab.AngularSpectrum(Jolab.Forward, grid_ns, rand(length(ns), length(ns)), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test_opt intensity(angspe)

# Test SpatialBeam type-stability
x = range(-100E-6, 100E-6, length = 100)
λ_spatial = 1550E-9
grid_spatial = Jolab.CartesianGrid(Jolab.X_Y_λ, x, x, λ_spatial)
@test_opt Jolab.SpatialBeam(Jolab.Forward, grid_spatial, rand(length(x), length(x)), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
beam_spatial = Jolab.SpatialBeam(Jolab.Forward, grid_spatial, rand(length(x), length(x)), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test_opt intensity(beam_spatial)

# Test SpatialBeam with Gaussian
x_gauss = range(-100E-6, 100E-6, length = 100)
grid_gauss = Jolab.CartesianGrid(Jolab.X_Y_λ, x_gauss, x_gauss, 1550E-9)
gaussian_spatial = Jolab.Gaussian(10E-6)
beam_spatial_gauss = Jolab.SpatialBeam(Jolab.Forward, grid_gauss, gaussian_spatial, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam_spatial_gauss) ≈ 1

# Test AngularSpectrum with Gaussian
x_ang = range(-0.5, 0.5, length = 100)
grid_angular_gauss = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, x_ang, x_ang, 1550E-9)
gaussian_angular = Jolab.Gaussian(10E-6)
beam_angular_gauss = Jolab.AngularSpectrum(Jolab.Forward, grid_angular_gauss, gaussian_angular, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam_angular_gauss) ≈ 1

# Test that isapprox returns false for beams with different coordinate types
x_spatial = range(-100E-6, 100E-6, length = 50)
y_spatial = range(-100E-6, 100E-6, length = 50)
λ_test = 1550E-9
grid_xy = Jolab.CartesianGrid(Jolab.X_Y_λ, x_spatial, y_spatial, λ_test)
beam_xy = Jolab.SpatialBeam(Jolab.Forward, grid_xy, ones(ComplexF64, 50, 50), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))

r_spatial = range(0, 100E-6, length = 50)
θ_spatial = range(0, 2π, length = 50)[1:end-1]
grid_rθ = Jolab.CylindricalGrid(Jolab.R_θ_λ, r_spatial, θ_spatial, λ_test)
beam_rθ = Jolab.SpatialBeam(Jolab.Forward, grid_rθ, ones(ComplexF64, 50, 49), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))

@test !(beam_xy ≈ beam_rθ)

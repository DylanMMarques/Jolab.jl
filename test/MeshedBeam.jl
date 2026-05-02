using Jolab, Test

# Test SpatialBeam construction with grid and explicit array
x = range(-100E-6, 100E-6, length=50)
y = range(-100E-6, 100E-6, length=50)
λ = 1550E-9
grid_spatial = Jolab.CartesianGrid(Jolab.X_Y_λ, x, y, λ)

e = ones(ComplexF64, 50, 50)
beam = Jolab.SpatialBeam(Forward, grid_spatial, e, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test size(beam.e) == (50, 50, 1, 1)

# Test with different medium
beam2 = Jolab.SpatialBeam(Forward, grid_spatial, e, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test size(beam2.e) == (50, 50, 1, 1)

# Test AngularSpectrum construction with grid and explicit array
nsr = range(0, 0.5, length=50)
nst = range(0, 2π, length=50)
λ_ang = 1550E-9
grid_angular = Jolab.CylindricalGrid(Jolab.NSR_NSθ_λ, nsr, nst, λ_ang)

e_ang = ones(ComplexF64, 50, 50)
beam_ang = Jolab.AngularSpectrum(Forward, grid_angular, e_ang, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test size(beam_ang.e) == (50, 50, 1, 1)

# Test view functionality
view_beam = @view beam[:,:,1:1,1:1]
@test view_beam.e ≈ beam.e
@test view_beam.mesh ≈ beam.mesh

# Test AngularSpectrum with Gaussian
gaussian = Jolab.Gaussian(10E-6)
beam_ang_gauss = Jolab.AngularSpectrum(Forward, grid_angular, gaussian, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test size(beam_ang_gauss.e) == (50, 50, 1, 1)

# Test SpatialBeam with R_θ_λ coordinates
r = range(0, 100E-6, length=50)
θ = range(0, 2π, length=50)[1:end-1]
λ_cyl = 1550E-9
grid_spatial_cyl = Jolab.CylindricalGrid(Jolab.R_θ_λ, r, θ, λ_cyl)

e_cyl = ones(ComplexF64, 50, 49)
beam_cyl = Jolab.SpatialBeam(Forward, grid_spatial_cyl, e_cyl, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test size(beam_cyl.e) == (50, 49, 1, 1)

# Test SpatialBeam R_θ_λ Gaussian intensity normalization
r_fine = range(0, 200E-6, length=200)
θ_fine = range(0, 2π, length=200)[1:end-1]
λ_cyl_fine = 1550E-9
grid_spatial_cyl_fine = Jolab.CylindricalGrid(Jolab.R_θ_λ, r_fine, θ_fine, λ_cyl_fine)
gaussian_cyl = Jolab.Gaussian(50E-6)
beam_cyl_gauss = Jolab.SpatialBeam(Forward, grid_spatial_cyl_fine, gaussian_cyl, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam_cyl_gauss) ≈ 1 rtol=0.1

# Test SpatialBeam Gaussian intensity normalization
x_fine = range(-50E-6, 50E-6, length=200)
y_fine = range(-50E-6, 50E-6, length=200)
λ_spatial = 1550E-9
grid_spatial_fine = Jolab.CartesianGrid(Jolab.X_Y_λ, x_fine, y_fine, λ_spatial)
gaussian_spatial = Jolab.Gaussian(10E-6)
beam_spatial_gauss = Jolab.SpatialBeam(Forward, grid_spatial_fine, gaussian_spatial, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam_spatial_gauss) ≈ 1

# Test AngularSpectrum Gaussian intensity normalization
nsr_fine = range(0, 0.5, length=500)
nst_fine = range(0, 2π, length=500)[1:end-1]
λ_ang_fine = 1550E-9
grid_angular_fine = Jolab.CylindricalGrid(Jolab.NSR_NSθ_λ, nsr_fine, nst_fine, λ_ang_fine)
gaussian_ang = Jolab.Gaussian(50E-6)
beam_ang_gauss_fine = Jolab.AngularSpectrum(Forward, grid_angular_fine, gaussian_ang, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam_ang_gauss_fine) ≈ 1 rtol=0.01

# Test AngularSpectrum with NSX_NSY_λ coordinates
nsx = range(-0.5, 0.5, length=50)
nsy = range(-0.5, 0.5, length=50)
λ_cart_ang = 1550E-9
grid_angular_cart = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsy, λ_cart_ang)

e_cart_ang = ones(ComplexF64, 50, 50)
beam_ang_cart = Jolab.AngularSpectrum(Forward, grid_angular_cart, e_cart_ang, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test size(beam_ang_cart.e) == (50, 50, 1, 1)

# Test AngularSpectrum NSX_NSY_λ Gaussian intensity normalization
nsx_fine = range(-0.3, 0.3, length=200)
nsy_fine = range(-0.3, 0.3, length=200)
λ_cart_ang_fine = 1550E-9
grid_angular_cart_fine = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx_fine, nsy_fine, λ_cart_ang_fine)
gaussian_cart_ang = Jolab.Gaussian(50E-6)
beam_ang_cart_gauss = Jolab.AngularSpectrum(Forward, grid_angular_cart_fine, gaussian_cart_ang, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test intensity(beam_ang_cart_gauss) ≈ 1 rtol=0.1

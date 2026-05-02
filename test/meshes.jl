using Jolab, Test
using Meshes, Unitful

x = LinRange(-.1, 0.15, 100)
y = LinRange(-.2, 0.3, 100)
z = LinRange(-4, 5, 100)
xyz = (x, y, z)
mesh = Jolab.CartesianGrid(length.(xyz), first.(xyz), step.(xyz))

@test_throws BoundsError Jolab.centroid(mesh, CartesianIndex(101,1,1))
@test_throws BoundsError Jolab.centroid(mesh, 100 + 100*100 + 100 * 100 *100)
@test_throws BoundsError Jolab.centroid(mesh, 0)
@test Jolab.centroid(mesh, 1) == first.(Jolab.get_ranges(mesh))
@test_opt Jolab.centroid(mesh, CartesianIndex(2,2,2))
@test_opt Jolab.volume(mesh, CartesianIndex(2,2,2))
@test_opt Jolab.centroid(mesh, 1)
@test_opt Jolab.volume(mesh, 1)

mesh_xyz = Jolab.CartesianGrid(Jolab.X_Y_Z, x, y, z)
@test mesh_xyz.lengths == length.(xyz)
@test mesh_xyz.spacing == step.(xyz)
@test mesh_xyz.origin == Jolab.X_Y_Z(first.(xyz)...)
@test all(Jolab.centroid(mesh_xyz, CartesianIndex(2, 3, 4)) .≈ (x[2], y[3], z[4]))
@test Jolab.get_ranges(mesh_xyz) == xyz

wav = 1550E-9
field = Jolab.MonochromaticSpatialBeam_gaussian(Forward, x, y, 10E-6, wav, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test (x,y) == Jolab.get_ranges(field.mesh)[1:2]
@test Jolab.centroid(field.mesh, 1)[3] == wav


grid = Jolab.CylindricalGrid((10, 1, 10), Jolab.R_θ_λ(0.0,0.0,0.0), (0.1, float(2π), 0.1))
@test_throws BoundsError Jolab.centroid(grid, CartesianIndex(11,1,1))
@test_throws BoundsError Jolab.centroid(grid, 10 + 10*10 + 1)
@test_throws BoundsError Jolab.centroid(grid, 0)
@test_throws ArgumentError Jolab.CylindricalGrid((10, 10, 10), Jolab.R_θ_λ(0.0, -0.1, 0.0), (0.1, 2π / 10, 0.1))

r = 0.1
@test Jolab.centroid(grid, 2) == 
    Jolab.centroid(grid, CartesianIndex(2,1,1)) == 
    [0.1, 0.0, 0.0]

@test Jolab.centroid(grid, CartesianIndex(10,1,10)) == 
    [0.9, 0.0, 0.9]
    
@test Jolab.volume(grid, 2) == 
    Jolab.volume(grid, CartesianIndex(2,1,1)) ==
    ustrip(u"m^3", measure(Cylinder(0.15)) - measure(Cylinder(0.05))) * 0.1
@test Jolab.volume(grid, 1) == 
    Jolab.volume(grid, CartesianIndex(1,1,1)) ==
    0.0

grid = Jolab.CylindricalGrid((10, 10, 10), Jolab.R_θ_λ(0.0,0,0), (0.1, 2π / 10, 0.1))
@test Jolab.centroid(grid, CartesianIndex(2,2,2)) == [0.1, 2π / 10, 0.1]

@test all(Jolab.centroid(grid, CartesianIndex(10,10,10)) .≈ (1 - 0.1, 2π - 2π/10, 1 - 0.1))

@test sum(i -> Jolab.volume(grid, CartesianIndex(2, i, 1)), 1:10)u"m^3" ≈ (measure(Cylinder(0.15)) - measure(Cylinder(0.05))) * 0.1
@test sum(i -> Jolab.volume(grid, CartesianIndex(10, i, 1)), 1:10)u"m^3" ≈ (measure(Cylinder(0.95)) - measure(Cylinder(0.85))) * 0.1

grid = Jolab.CylindricalGrid((10, 10, 10), Jolab.R_θ_λ(0.05,0,0), (0.1, 2π / 10, 0.1))
@test sum(i -> Jolab.volume(grid, i), eachindex(grid))u"m^3" ≈ measure(Cylinder(1))
@test_opt Jolab.centroid(grid, CartesianIndex(2,2,2))
@test_opt Jolab.volume(grid, CartesianIndex(2,2,2))
@test_opt Jolab.centroid(grid, 1)
@test_opt Jolab.volume(grid, 1)

# Test coordinate validation for CartesianGrid
@test_throws ArgumentError Jolab.CartesianGrid(Jolab.R_θ_λ, range(0, 1), range(0, 2π), 1550E-9)
@test_throws ArgumentError Jolab.CartesianGrid(Jolab.NSR_NSθ_λ, range(0, 0.5), range(0, 2π), 1550E-9)

# Test coordinate validation for CylindricalGrid
@test_throws ArgumentError Jolab.CylindricalGrid(Jolab.X_Y_λ, range(0, 1), range(0, 2π), 1550E-9)
@test_throws ArgumentError Jolab.CylindricalGrid(Jolab.NSX_NSY_λ, range(-0.5, 0.5), range(-0.5, 0.5), 1550E-9)

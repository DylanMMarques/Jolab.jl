using Jolab, Test
using Meshes, Unitful

grid = Jolab.CylindricalGrid((10, 1, 10), Jolab.R_θ_λ(0.0, 0.0, 0.0), (0.1, float(2π), 0.1))

r = 0.1
Jolab.centroid(grid, 2)
@test Jolab.volume(grid, 2) ==
      Jolab.volume(grid, CartesianIndex(2, 1, 1)) ==
      ustrip(u"m^3", measure(Cylinder(0.15)) - measure(Cylinder(0.05))) * 0.1
@test Jolab.volume(grid, 1) == Jolab.volume(grid, CartesianIndex(1, 1, 1)) == 0.0

grid = Jolab.CylindricalGrid((10, 10, 10), Jolab.R_θ_λ(0.0, 0, 0), (0.1, 2π / 10, 0.1))
@test Jolab.centroid(grid, CartesianIndex(2, 2, 2)) == [0.1, 2π / 10, 0.1]

@test all(
    Jolab.centroid(grid, CartesianIndex(10, 10, 10)) .≈ (1 - 0.1, 2π - 2π/10, 1 - 0.1),
)

@test sum(i -> Jolab.volume(grid, CartesianIndex(2, i, 1)), 1:10)u"m^3" ≈
      (measure(Cylinder(0.15)) - measure(Cylinder(0.05))) * 0.1
@test sum(i -> Jolab.volume(grid, CartesianIndex(10, i, 1)), 1:10)u"m^3" ≈
      (measure(Cylinder(0.95)) - measure(Cylinder(0.85))) * 0.1

grid = Jolab.CylindricalGrid((10, 10, 10), Jolab.R_θ_λ(0.05, 0, 0), (0.1, 2π / 10, 0.1))
@test sum(i -> Jolab.volume(grid, i), eachindex(grid))u"m^3" ≈ measure(Cylinder(1))

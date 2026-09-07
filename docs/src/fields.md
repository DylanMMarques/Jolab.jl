# Representation of optical fields

To represent an optical field in a simulation, Jolab uses various
A optical field can be represented in different manners in Jolab. Probably the most common representation is based on the spatial profile of the field in the spatial domain, where the electric field values are defined at specific physical locations in space. Such a field can be created in Jolab using the `SpatialBeam` function and requires a mesh grid defined in the spatial domain, the `X_Y_λ` coordinates.

```@example fields


using Jolab
# Create a spatial beam with a Gaussian profile
x = range(-50E-6, 50E-6, length=100) # 100 points from -50 microns to 50 microns
y = range(-50E-6, 50E-6, length=100)
wavelength = 500E-9 # 500 nm
grid = CartesianGrid(X_Y_λ, x, y, wavelength) # Create a Cartesian grid in the spatial domain
field = SpatialBeam(Forward, grid, Gaussian(10E-6), Medium(1.0), ReferenceFrame((0,0,0), (0,0,0))) # Create a spatial beam with a Gaussian profile
```

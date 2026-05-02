# First simulation

Jolab can be used to simulate complex optical systems. The basis of Jolab is to decompose the system into a sequence of optical elements, and then propagate the light through all the elements. Internally, Jolab uses dedicated solvers for each optical element which allow to achieve high accuracy and speed. The example shows a simple example of a collimated beam focused through a microscope slide using a lens. The aim of the simulation is to measure the reflectance and transmittance by the system. First, the lens is defined:

```julia first_simulation
using Jolab

# Define the optical system

# Refractive index before and after the lens
lens_media = (Medium(1.0), Medium(1.0))
focal_length = 10E-3 # 10 mm
numerical_aperture = 0.5

# Define the reference frame specifying the position and orientation of the lens
lens_frame = ReferenceFrame((0,0,0), (0,0,0)) 

lens = Lens(focal_length, numerical_aperture, lens_media, lens_frame)
```

With the lens defined, we can now define the microscope slide placed on the focal plane of the lens. The microscope slide is defined as a dielectric stack with a thickness of 1 mm and a refractive index of 1.5 (glass).

```julia first_simulation
# Define the microscope slide
thickness = [1E-3] # 1 mm

# The first and lens_media[2] are the same medium as the lens (air in this case), the second is the glass of the microscope slide
stack_refractive_index = [lens_media[2], Medium(1.5), lens_media[2]]  # glass

slide_frame = lens_frame + ReferenceFrame((0,0,focal_length), (0,0,0)) # Place the slide at the focal plane of the lens
slide = DielectricStack(stack_refractive_index, thickness, slide_frame)
```

The code currently defines the optical elements of the system. The next step is to define the light source illuminating the optical system. In this example, we will use a monochromatic scalar collimated beam. 

## Defining the coordinate grid

In Jolab, optical fields are represented on a computational grid of coordinates. For spatial beams, this is a 3D Cartesian grid defined by three coordinates: **x** (horizontal), **y** (vertical), and **λ** (wavelength). The `CartesianGrid` class manages this computational domain:

- **x and y ranges**: Define the spatial extent and resolution of the field in the transverse plane
- **wavelength (λ)**: Typically a single wavelength for monochromatic simulations, but can span a range for multi-wavelength analysis

Each point on the grid represents the field value at that spatial location and wavelength. The grid stores both the field values and the coordinates, allowing efficient numerical propagation and interaction calculations.

```julia first_simulation
x = range(-5E-3, 5E-3, length=100) # 10 mm beam diameter
y = range(-5E-3, 5E-3, length=100)

# Define the coordinate system for the simulation (x, y, and wavelength)
grid= Jolab.CartesianGrid(Jolab.X_Y_λ, x, y, range(500E-9, step=1E-15, length = 1)) 

field_frame = ReferenceFrame((0,0,-focal_length), (0,0,0)) # Place the field in the back focal plane of the lens

gaussian_e_field(coord) = exp(-((coord[1]^2 + coord[2]^2)/(2*(2.5E-3)^2))) # Gaussian beam with a waist of 2.5 mm
e = gaussian_e_field.(centroid.(grid))

field = SpatialBeam(grid, e, lens_media[1], field_frame)
```

The field is now defined, and we can propagate it through the system. The propagation is performed in order of the optical elements defined in the system. In this example, the field is first propagated through the lens to calculate the field at the focal plane. For this, we use the light_interaction function which calculates the reflected and transmitted field by an optical component for a given incident field.

```julia first_simulation
(rfield, tfield) = light_interaction(lens, field)
```

The light_interaction function returns the reflected and transmitted fields by the lens. In this case, no field is reflected by the lens because the model of lens does not include reflections (simulates an ideal lens). The transmitted field is the field at the focal plane of the lens, which illuminates the microscope slide. The next step is to propagate the transmitted field through the microscope slide.

```julia first_simulation
(rfield_slide, tfield_slide) = light_interaction(slide, tfield)
```





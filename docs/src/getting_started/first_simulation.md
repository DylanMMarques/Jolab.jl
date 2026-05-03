# First simulation

Jolab can be used to simulate complex optical systems. The basis of Jolab is to decompose the system into a sequence of optical elements, and then propagate the light through all the elements. Internally, Jolab uses dedicated solvers for each optical element which allow to achieve high accuracy and speed. The example shows a simple example of a collimated beam focused through a microscope slide using a lens. The aim of the simulation is to measure the reflectance and transmittance by the system. First, the lens is defined:

```@example first_simulation
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

```@example first_simulation
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

```@example first_simulation
x = range(-5E-3, 5E-3, length=100) # 10 mm beam diameter
y = range(-5E-3, 5E-3, length=100)

# Define the coordinate system for the simulation (x, y, and wavelength)
λ = 500E-9  # 500 nm
grid = Jolab.CartesianGrid(Jolab.X_Y_λ, x, y, λ)
```

## Probing grid positions

Once the grid is defined, the coordinates at each grid point can be retrieved using the `centroid` function. For example, the first grid point has (x,y,wavelength) coordinates of:
```@example first_simulation
Jolab.centroid(grid, 1) # returns the coordinates of the first grid point
```

The `centroid` can be used to evaluate to calculate the electric field of a Gaussian beam at each grid point, as:

```@example first_simulation
gaussian_e_field(coord) = exp(-((coord[1]^2 + coord[2]^2)/(2*(2.5E-3)^2))) # Gaussian beam with a waist of 2.5 mm
e = [gaussian_e_field(Jolab.centroid(grid, i)) for i in eachindex(grid)]
```

Jolab is a framework to simulate complex optical systems and therefore, the position of the field relative to each optical element is important. In Jolab, the position and orientation of optical elements and fields are defined using reference frames. A reference frame is defined by a position vector and an orientation (rotation) vector. The position vector specifies the location of the reference frame in space, while the orientation vector specifies the rotation of the reference frame relative to a global coordinate system.

```@example first_simulation
# Define the reference frame for the field
field_frame = ReferenceFrame((0,0,-focal_length), (0,0,0)) # Place the field in the back focal plane of the lens
```

Finally, we can combine the grid, electric field, and reference frame to define the spatial beam representing the light source in the simulation. The `SpatialBeam` class represents a spatial beam of light defined on a Cartesian grid, with a specified electric field and reference frame. We specify `Forward` to indicate the direction of propagation (the light travels forward through the optical system).
```@example first_simulation
field = SpatialBeam(Forward, grid, e, lens_media[1], field_frame)
```

The field is now defined, and we can propagate it through the system. The propagation is performed in order of the optical elements defined in the system. In this example, the field is first propagated through the lens to calculate the field at the focal plane. For this, we use the light_interaction function which calculates the reflected and transmitted field by an optical component for a given incident field.

```@example first_simulation
(rfield, tfield) = light_interaction(lens, field)
```

The light_interaction function returns the reflected and transmitted fields by the lens. In this case, the field reflected by the lens is null because the model of lens does not include reflections (simulates an ideal lens).

```@example first_simulation
intensity(rfield) 
```

Let's now focus on the transmitted field `tfield`. The intensity of the transmitted field is the same as the incident field, which is expected for an ideal lens that does not introduce losses.

```@example first_simulation
intensity(tfield) / intensity(field) 
```

Looking at the reference frame of the transmitted field, we can see that the reference frame changed to the focal plane of the lens. This is due to the solver of the lens, which propagates the field from the back focal plane to the focal plane of the lens. 

```@example first_simulation
tfield.frame
```

The field incident upon the lens was defined in the spatial domain, meaning that the electric field values were defined at specific physical locations in space. The transmitted field is now represented in the angular spectrum domain, which is the Fourier-space representation of the field. The representation of the field can be checked by looking at the coordinates of the mesh grid of the field (`X_Y_λ` for spatial domain, `NSX_NSY_λ` for angular spectrum domain - more info in field representations).
```@example first_simulation
(typeof(field.mesh), typeof(tfield.mesh))
```

Similarly to the reference frames, the solver of each element determine the type of field representation of the incident, transmitted, and reflected fields. The information about the field representations needed for each solver is available in the documentation of each solver.

## Converting the field back to spatial domain with the Fourier operator

The transmitted field `tfield` is currently represented in the angular spectrum domain (`NSX_NSY_λ`). To analyze or visualize the field at specific spatial locations, we can convert it back to the spatial domain (`X_Y_λ`). This is accomplished using a Fourier transform operator, which performs an inverse FFT to compute the spatial field distribution from the angular spectrum.

The `FourierFFT` operator performs this transformation:

```@example first_simulation
# Apply Fourier transform to convert from angular spectrum to spatial domain
(_, field_spatial) = light_interaction(Jolab.FourierFFT(), tfield)
```

The returned `field_spatial` is now a field in spatial domain coordinates, showing the electric field at specific physical locations. We can verify this by checking the coordinate types:

```@example first_simulation
# The input was in angular spectrum coordinates (NSX_NSY_λ)
# The output is in spatial coordinates (X_Y_λ)  
(typeof(tfield.mesh), typeof(field_spatial.mesh))
```

For example, we can compute the intensity distribution in spatial domain to see how the Gaussian beam is focused by the lens:

```@example first_simulation
using CairoMakie

heatmap(abs2.(field_spatial.e[:,:,1,1]))
```

## Propagating through the microscope slide
The transmitted field is now ready to interact with the next optical element in the system. The next step is to propagate the transmitted field through the microscope slide.

```@example first_simulation
(rfield_slide, tfield_slide) = light_interaction(slide, tfield)
```

We could for example analyse the phase profile introduced by the microscope slide by looking at the phase of the transmitted field `tfield_slide`:

```@example first_simulation
heatmap(angle.(tfield_slide.e[:,:,1,1]))
```

In some scenarios, it may be interesting to propagate the reflected field `rfield_slide` back through the lens to analyze the reflected light at the back focal plane. This can be done by applying the `light_interaction` function again, using the reflected field as input and the lens as the optical element:

```@example first_simulation
(rfield_slide_backward), rfield_slide_forward) = light_interaction(lens, rfield_slide)
```

Here, the `rfield_slide_backward` is the field propagating backward thought the optical system, while `rfield_slide_forward` is the field propagating forward. Consequently, the `rfield_slide_backward` is the field transmitted by the lens while the `rfield_slide_forward` is the field reflected by the lens. This is the convention used in Jolab, the `light_interaction` function always returns two outputs, the first one which is the field propagating backward and the second one which is the field propagating forward.


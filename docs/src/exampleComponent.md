# Sequential modelling
This is an example of using Jolab to simulate the light propagation on optical system. Jolab has models implemented to simulate various optical components but to simulate optical system the user must apply propagate the light by the multiple optical components, i.e. sequential modelling. We start by modelling a simple optical system: two lenses placed sequentially on the 4f configuration.

The first step is to load the Jolab into Julia
```julia
using Jolab
using Plots; plotly()
```

## Defining the optical setup
We will now create the first lens which focus the collimated beam. A lens is defined based on its focal length (f), its numerical aperture (na) and its referenceframe, i.e. the position in space where the lens is placed.
```julia
f1 = 30E-3 # focal length in meters of the lens
na1 = 1 # Numerical aperture of the lens
ref1 = ReferenceFrame(0, 0, f1) # Position in meters of the lens in meters (x,y,z)
lens1 = Lens(f1, na1, ref1); # Creation of the lens
```

The next step is to create the second lens which will re-collimate the beam. As before, we define its focal length, na and referenceframe. We place the lens such that the focal plane of both lenses is the same.
```julia
f2 = 60E-3 # focal length in meters of the lens
na2 = 1 # Numerical aperture of the lens
ref2 = ReferenceFrame(0, 0, 2f1 + f2) # Position in meters of the lens
lens2 = Lens(f2, na2, ref2); # Creation of the lens
```
We also need to specify the collimated beam incident upon the lens. In this example, we will consider a Gaussian beam in the back focal plane of the lens. A Gaussian beam is defined based on the beam waist (σ), wavelength (λ) and referenceframe:
```julia
λ = 1550E-9 # wavelength in meters
σ = 1000E-6 # Gaussian beam waist in meters. Defined based on where the intensity decreases to 1/e²
x = σ .* range(-1.5, 1.5, length = 128); # Position x in meters where the field is evaluated
y = σ .* range(-1.5, 1.5, length = 128); # Position y in meters where the field is evaluated
n = Jolab.refractiveindex_air(printBool= false); # Refractive index of the medium
dir = 1; # Direction of propagation of the field (forward 1 or backward -1)
ref = ReferenceFrame(0 , 0, 0); # Position in mirrors of the field (x,y,z)
fieldBackFocalPlane = FieldSpace_gaussian(x, y, σ, λ, n(λ), dir, ref); # Creation of the field
plot(fieldBackFocalPlane) # Display the field intensity
```

## Propagating the field by the optical setup

At this point, all the optical components and initial fields of the simulations are created. The next step is to propagate the field by the optical components. This is done with the function ```(backward_field, forward_field) = lightinteraction(Optical_component, incident_field)``` where ```backward_field``` is the field propagating backward after the interaction with the optical component and ```forward_field``` is the field propagating forward. For the case the an optical component does not reflect light, for example as an ideal lens, one of the field will be nothing.

Following the script example, the simulation of the first lens is made as:
```julia
(tmp, angspeFocalPlane) = lightinteraction(lens1, fieldBackFocalPlane); # Calculates the field in the focal plane of the lens
```

We have calculated the field in the focal plane of both lenses. We can now propagate the field by the second lens. We need to repeat the same approach as before:
```julia
(tmp, fieldBackFocalPlane_2) = lightinteraction(lens2, angspeFocalPlane);
plot(fieldBackFocalPlane_2)
```
We have finally calculated how light propagates by two lenses. We can see that the beam is magnified as it propagates by the two lens due to the different focal length.
This example shows the basic concept of the toolbox - the modularity. This toolbox is built to allow its users to calculate how light interacts with single optical components and, by putting code together, how light propagates on optical setups. This simple example only simulates lens but any lens could be replace by any other optical component such as, fibre, mirror, interferometers, etc...

# Multimode Fibre
This example shows the usage of the toolbox to calculate light propagation in a cilindre. The normal profile of the refractive index is defined by the refractive index of the cilindre and the refractive index of the surronding media. The toolbox allows to simulate light propagation inside the cilindre assuming the beam incident upon its face. This example will show how to use the toolbox to simulate light propagation by a Multi Mode Fibre. As a aim, we will consider an field incident upon one tip of the fibre an propagate it to the other side of the fibre.
We start by loading the toolbox and a plotting package.
```julia
using Jolab
using Plots; plotly()
```
## Defining the optical setup
We need now to define the fibre itself based on the refractive index of the core as the fibre Numerical Aperture.
```julia
r_fib = 50E-6; # radius of the fibre core in meters
na_fib = .1; # numerical aperture of the fibre
ref_fib = ReferenceFrame(0,0,0,0,0); # position and orientaion of the fibre tip
length_fib= 1; # length of the fibre. The referenceframe second fibre tip is on the other side of the cilindre. I.e., in this case the referenceframe of the second tip is ReferenceFrame(0,0,length_fib,0,0);
nambient = Jolab.refractiveindex_air(printBool = false); # refractive index before the fibre
ncore_fib = Jolab.refractiveindex_fusedsilica(printBool = false); # refractive index of the fibre core
fibre = CircularStepIndexFibre(r_fib, na_fib, ncore_fib, nambient, ref_fib, nambient, length_fib); # Initialization of the fibre
```
In the previous lines we create the multi mode fibre.

## Modes calculation
Light inside the fibre propagate as modes, i.e. the beam can be decomposed as a sum of modes which replicate in space along the fibre. The modes can be calculated, for a given wavelength, using the function findmodes! which will solve numerically the wave equation for the given fibre and wavelength.

```julia
λ = 1550E-9; # Defines the wavelength for the simulation
Jolab.findmodes!(fibre, λ); # Calculates the possible modes inside the fibre
modes = fibre.modes[1]; # We choose the modes for the wavelength considered. If we calculate the modes for more than 1 wavelength they would be saved on fibre.modes[i].
plot(modes[50], type = :real, fill = true, c = :bluesreds) # plot the mode

```
We can plot the shape of the modes using the function plot as above. The number of modes found for a given fibre/wavelength can be calculated using the function numberofmodes(modes).
The toolbox saves the modes inside the variable fibre and therefore they need to be calculated only once for a given fibre and wavelength. The modes are automatically calculated (if not precalculated) when using the function lightinteraction(fibre, field) and therefore it is not needed to call the function findmodes!.

## Coupling light inside the fibre
To calculate the light coupling inside the fibre we first need to define the field insident upon the fibre. In this case, we will consider a Gaussian beam at the tip of the fibre but not centered
```julia
x = range(-50E-6, 50E-6, length = 100); # x spatial coordinates where the field and modes are evaluated in meters. Make sure that the sampling of the field and modes is good enough within the region of interest, i.e. where the field is incident
y = range(-50E-6, 50E-6, length = 100); # y spatial coordinates where the field and modes are evaluated in meters. Make sure that the sampling of the field and modes is good enough within the region of interest, i.e. where the field is incident
σ = 50E-6; # Beam waist of the Gaussian beam defined by 1/e^2 of the Gaussian intensity profile
ref_field = ReferenceFrame(50E-6,0,0,0,0); # ReferenceFrame of the field with the position x,y,z in meters. The beams is offset by 50 μm to the center of the fibre
field = FieldSpaceScalar_gaussian(x, y, σ, λ, nambient(λ), 1, ref_field); # Initialization of the field
plot(field, fill=true, c = :amp, type = :abs) # plot the incident field upon the fibre
```
We have defined all the necessary parameters to perform the simulation. The next line of code calculates the light coupling inside the fibre.
```julia
(tmp, fieldmodes) = lightinteraction(fibre, field);
```
fieldmodes is a struct FieldModes which describes a field based into the modes of the multi mode fibre. This allows to propagate the light inside the fibre until the second tip. In order to calculate the field that goes out of the fibre (in the order side) we need to use the function lightinteraction again with the spatial coordinates where the field will be evaluated. The field propagation is performed automatically.

```julia
plot(fieldmodes, c = :amp, fill = true, type = :abs) # Plot hte field outside the fibre
```
We can be see that field_tip2 is defined on the other side of the tip (look at the referenceframe field_tip2.ref).
```julia
intensity(fieldmodes) / intensity(field) # Difference between the intensity of the light incident and the light transmitted.
```
We can also calculate the part of the light that is coupled inside the fibre by making the ratio between the transmitted and incident beam intensity.

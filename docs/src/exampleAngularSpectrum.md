# Angular Spectrum
This is an example from Jolab to familiarise an user with the concept of angular spectrum.

The first step is to load the Jolab and a plotting package
```julia
using Jolab
using Plots; plotly()
```
## Defining a light field in space
We start by defining a focused Gaussian beam in space.
```julia
λ = 1550E-9 # wavelength in meters
σ = 5E-6 # Gaussian beam waist in meters. Defined based on where the intensity decreases to 1/e²
x = σ .* collect(range(-1.5, 1.5, length = 128)); # Position x in meters where the field is evaluated
y = σ .* collect(range(-1.5, 1.5, length = 128)); # Position y in meters where the field is evaluated
n = Jolab.refractiveindex_air(printBool= false); # Refractive index of the medium
dir = 1; # Direction of propagation of the field (forward 1 or backward -1)
ref = ReferenceFrame(0 , 0, 0); # Position in mirrors of the field (x,y,z)
fieldspace = FieldSpace_gaussian(x, y, σ, λ, n(λ), dir, ref); # Creation of the field
plot(fieldspace, fill = true, c = :amp) # Display the field intensity
```
The field calculated on the section above is defined based on the field distribution in space, i.e. we know the amplitude and phase of the field in the positions (x,y). If we would put a camera in front of the field, the data acquired would be the graph plotted above.

## Transforming to an Angular Spectrum
An optical field can be interpreted in different ways. For example, instead of assuming our field as a distribution in space we can assume the field as a distribution of plane waves with different directions of propagation, i.e. the field is composed as a distribution of plane waves with complex amplitude E(n sx, n sy) propagating with direction (nsx,nsy). nsx and nsy are the spatial frequencies and related with the angle of propagation (θ,ϕ) by (n sin(θ) cos(ϕ), n sin(θ)sin(ϕ)) where n is the medium refractive index. This way of interpreting the field is commonly called angular spectrum.

The angular spectrum of a field is calculated by 2D Fourrier Transform the field distribution in space. Jolab has already a function to calculate the angular spectrum based on a field distribution in space.
```julia
nsx = range(-.5, .5, length = 128); #x spatial frequencies where the field is evaluated
nsy = range(-.5, .5, length = 128); #y spatial frequencies where the field is evaluated
angspe = FieldAngularSpectrum_fromspace(fieldspace, nsx, nsy); # calculate the angular spectrum from the field distribution in space
plot(angspe, fill = true, c = :amp)
```
The color scale of the graph above represent the amplitude of the plane wave and the x,y axis the nsx and nsy spatial frequencies, i.e. the direction of the respective plane wave. This can be illustrator by changing the distribution of the plane wave amplitude so that the field just has one angular component.

```julia
angspeTest = deepcopy(angspe) # create a copy of the angspe
angspeTest.e_SXY .= 0.; # Set all plane waves amplitude to 0
angspeTest.e_SXY[1,64,64] = 1; # Choose a plane a set amplitude 1. For the index [64,64] the plane wave is propagating parallel to the z axis
z = copy(x) # Create the z coordinates where the field is evaluated
fieldspaceTest = angspeto3Dspace(angspeTest, x, [0.], z); # evaluate the field in space. It returns a 4D matrix where the first dimension is the vectorial dimension and then x,y,z
plot(x, z, real.(fieldspaceTest[1,:,1,:]), fill = true, c = :bluesreds) # Plot the real part of the data
```
The code above shows that if we only have a single component defined on the angular spectrum the respective field in space corresponds to a plane wave. We recommend the user to play with the code and changing the index which will change the direction of the plane wave.

## Transforming back to the field in space
From an angular spectrum we can calculate the field in space using the respective inverse transformation.
```julia
fieldspace2 = FieldSpace_fromangspe(angspe, x, y);
plot(fieldspace2, fill = true, c = :amp)
```
The code above transformed the field back to its distribution in space and ,as we can see, fieldspace2 is identical to fieldspace which was initially calculated.
It is important to understand that defining a field in space or as a angular spectrum is irrelevant as the information about the field is the same but just in another format. The user can always switch between the spatial field and the angular spectrum using the respective function (angspetofield and fieldtoangspe). The only errors induced are due to evaluating the 2D fourrier transform numerically.

## The application of the Angular Spectrum

The big advantages of the angular spectrum is that there are theories that allow to calculate how a plane wave interacts with different optical components. The most basic example is the plane wave propagation in free space. It is challenging to calculate how a focused beam propagates in free space based on its distribution in space. However, it is very well known how a plane wave propagates in space and therefore, by considering the propagation of all plane waves, how an angular spectrum propagates.
```julia
z = 2 .* x; #  position in the z axis where the field is evaluated
fieldpropagation = angspeto3Dspace(angspe, x, [0.], z); # evaluate the field in space. It returns a 4D matrix where the first dimension is the vectorial dimension and then x,y,z
plot(x, z, real.(fieldpropagation[1,:,1,:]), fill = true, c = :bluesreds) # Plot the real part of the data
```
The code above shows the field propagation using the angular spectrum in which we can see the Gaussian field first converging and then diverging as it propagates. This approach is valid for any arbitrary beam.

## Using both formalism
Jolab was developed to use both formalism as it increases the versability. For some optical components, such multilayer structures, spheres, how light interacts with them is know for an incident plane wave and therefore the angular spectrum is used. On other scenarios, such as cameras and detectors, the toolbox requires the field in space to calculate the signal measured.

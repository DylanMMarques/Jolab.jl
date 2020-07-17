using Jolab
using Plots
plotly()

x_slm = range(-10E-3, 10E-3, length = 25); # x coordinates of the SLM phase mask pixels in meters
y_slm = range(-10E-3, 10E-3, length = 25); # y coordinates of the SLM phase mask pixels in meters
mask_slm = ones(ComplexF64, length(x_slm), length(y_slm)); # initial mask of the SLM. Each position of the matrix phase[iX, iY] is the modulation value for the SLM pixel at position x[iX], y[iY].
ref_slm = ReferenceFrame(0.,0,0); # x,y,z position of the SLM center in meters
slm = SpatialLightModulator(x_slm, y_slm, mask_slm, ref_slm); # Initialization of the SLM

f_1 = 50E-3
na_1 = 1.
ref_1 = ReferenceFrame(0,0,f_1)
lens_1 = Lens(f_1, na_1, ref_1)

r_fib = 50E-6; # radius of the fibre core in meters
na_fib = .1; # numerical aperture of the fibre
ref_fib = ReferenceFrame(0,0,2f_1); # position and orientaion of the fibre tip
length_fib= .1; # length of the fibre. The referenceframe second fibre tip is on the other side of the cilindre. I.e., in this case the referenceframe of the second tip is ReferenceFrame(0,0,length_fib,0,0);
nambient = Jolab.refractiveindex_air(printBool = false); # refractive index before the fibre
ncore_fib = Jolab.refractiveindex_fusedsilica(printBool = false); # refractive index of the fibre core
fibre = CircularStepIndexFibre(r_fib, na_fib, ncore_fib, nambient, ref_fib, nambient, length_fib); # Initialization of the fibre

ref_3 = ReferenceFrame(0,0,.11)
phatomLayer = MultilayerStructure([1,1], zeros(0), ref_3)
setup = [lens_1, Jolab.FourierTransform(), fibre, Jolab.FourierTransform(), phatomLayer]

λ = 1550E-9; # wavelength in meters of the field
nambient = Jolab.refractiveindex_air(printBool = false); # refrative index of the ambient medium
dir = 1; # direction of the field
x = range(-10E-3, 10E-3, length = 100); # x coordinates in meters where the field is evaluated
y = range(-10E-3, 10E-3, length = 100); # y coordinates in meters where the field is evaluated
ref_f = ReferenceFrame(0,0,0.)
field = FieldSpace_uniform(x, y, λ, nambient(λ), dir, ref_f); # Initialization of field with uniform intensity.

@time coef = Jolab.coefficient_general(setup, field)

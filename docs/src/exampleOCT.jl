using Jolab, Plots

### Building the OCT system:

sr = range(0, 0.2, length = 100)

n_air = Jolab.refractiveindex_air()

# SingleModeFibre(Mode field diamenter, refractive index outside the fibre, direction at which the fibre is pointing (1 or -1))
smf = SingleModeFibre(9.2E-6, n_air, 1, ReferenceFrame(0,0,0))

## Reference OCT branch
f_col_ref = 18.67E-3
z_scan = 0E-6 # Parameter to scan the mirror
col_ref = Lens(f_col_ref, 0.14, ReferenceFrame(0, 0, f_col_ref))
mirror_ref = Mirror(1, n_air, n_air, ReferenceFrame(0, 0, 2f_col_ref + z_scan))

r_R = sr .* f_col_ref
sr_colimated = range(0, 0.001, length = 100)
fourier_ref = FourierTransform(r_R, r_R, sr_colimated, sr_colimated)

## Sample arm
f_col_sam = 18.67E-3
col_sam = Lens(f_col_sam, 0.14, ReferenceFrame(0, 0, f_col_sam))
f_obj_sam = 50E-3
obj_sam = Lens(f_obj_sam, 1, ReferenceFrame(0,0,2f_col_sam + f_obj_sam))

# The sample will be a mirror with a pinhole in front. This will give you the point spread function of the OCT
mirror_sam = Mirror(1, n_air, n_air, ReferenceFrame(0,0,2f_col_sam + 2f_obj_sam))

sr_focus = sr .* f_col_sam / f_obj_sam
r_focused = range(0, 100E-6, length = 100)
fourier_sam = FourierTransform(r_focused, r_focused, sr_focus, sr_focus)

# Calculating the optical field
λ = 1310E-9

# Reference arm 
field = Jolab.FieldAngularSpectrumScalarRadialSymmetric_fromfibre(smf, sr, λ)
(tmp, tfield) = lightinteraction(col_ref, field)
(tmp, tfield) = lightinteraction(fourier_ref, tfield)
(rfield, tmp) = lightinteraction(mirror_ref, tfield)
(rfield, tmp) = lightinteraction(fourier_ref, rfield)
(rfield, tmp) = lightinteraction(col_ref, rfield)
mode_coupling_ref = Jolab.signal_complex(smf, rfield)


function apply_pihnole(field::Jolab.AbstractFieldSpace)
    field.e_SXY[3:end] .*= 0
    field
end

# Sample arm
field = Jolab.FieldAngularSpectrumScalarRadialSymmetric_fromfibre(smf, sr, λ)
(tmp, tfield) = lightinteraction(col_sam, field)
(tmp, tfield) = lightinteraction(obj_sam, tfield)
(tmp, tfield) = lightinteraction(fourier_sam, tfield)
# tfield = apply_pihnole(tfield) # You could add a pinhole here to study the Point spread function
(rfield, tmp) = lightinteraction(mirror_sam, tfield)
(rfield, tmp) = lightinteraction(fourier_sam, rfield)
(rfield, tmp) = lightinteraction(obj_sam, rfield)
(rfield, tmp) = lightinteraction(col_sam, rfield)
mode_coupling_sam = Jolab.signal_complex(smf, rfield)

measurement = abs2(mode_coupling_ref + mode_coupling_sam)

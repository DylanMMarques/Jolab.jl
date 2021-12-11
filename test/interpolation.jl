using Jolab, Interpolations, Test

nsx_X = range(-.2, .2, length = 10)
nsy_Y = range(-.2, .2, length = 10)
nsx_itp = range(-.22, .22, length = 20)
nsy_itp = range(-.22, .22, length = 20)

field = FieldAngularSpectrumScalar_gaussian(nsx_X, nsy_Y, 50E-6, 1550E-9, 1, 1, ReferenceFrame(0,0,0.))
changereferenceframe!(field, ReferenceFrame(10E6, 0, 0,))
field_itp = Jolab.interpolate(field, nsx_itp, nsy_itp)
itpref = LinearInterpolation((nsx_X, nsy_Y), reshape(field.e_SXY, 10, 10), extrapolation_bc = 0.)
@test all(field_itp.e_SXY .== vec(itpref.(nsx_itp, nsy_itp')))

field = FieldAngularSpectrumScalar_gaussian(nsx_X, nsy_Y, 50E-6, 1550E-9, 1, -1, ReferenceFrame(0,0,0.))
changereferenceframe!(field, ReferenceFrame(10E6, 0, 0,))
field_itp = Jolab.interpolate(field, nsx_itp, nsy_itp)
itpref = LinearInterpolation((nsx_X, nsy_Y), reshape(field.e_SXY, 10, 10), extrapolation_bc = 0.)
@test all(field_itp.e_SXY .== vec(itpref.(nsx_itp, nsy_itp')))

x_X = range(-50E-6, 50E-6, length = 10)
y_Y = range(-50E-6, 50E-6, length = 10)
x_itp = range(-55E-6, 55E-6, length = 20)
y_itp = range(-55E-6, 55E-6, length = 20)

field = FieldSpaceScalar_gaussian(x_X, y_Y, 50E-6, 1550E-9, 1, 1, ReferenceFrame(0,0,0.))
field_itp = Jolab.interpolate(field, x_itp, y_itp)
itpref = LinearInterpolation((x_X, y_Y), reshape(field.e_SXY, 10, 10), extrapolation_bc = 0.)
@test all(field_itp.e_SXY .== vec(itpref.(x_itp, y_itp')))

field = FieldSpaceScalar_gaussian(x_X, y_Y, 50E-6, 1550E-9, 1, -1, ReferenceFrame(0,0,0.))
field_itp = Jolab.interpolate(field, x_itp, y_itp)
itpref = LinearInterpolation((x_X, y_Y), reshape(field.e_SXY, 10, 10), extrapolation_bc = 0.)
@test all(field_itp.e_SXY .== vec(itpref.(x_itp, y_itp')))
return true
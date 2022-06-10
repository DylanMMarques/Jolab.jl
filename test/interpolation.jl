stop
using Jolab, Interpolations, Test

nsx_X = range(-.2, .2, length = 10)
nsy_Y = range(-.2, .2, length = 11)
nsx_itp = range(-.22, .22, length = 15)
nsy_itp = collect(range(-.22, .22, length = 14))

x_X = collect(range(0E-6, 100, length = 10))
x_itp = collect(range(0E-6, 100, length = 10))

for (field_type, x, xitp,) in ((:FieldSpaceScalarRadialSymmetric_gaussian, :x_X, :x_itp), (:FieldAngularSpectrumScalarRadialSymmetric_gaussian, :nsx_X, :nsx_itp))
    eval(quote
    @show $xitp
        itp_coef = Jolab.Interpolation($xitp, nsy_itp)
        field = Jolab.$field_type($x, 25E-6, 1550E-9, 1, 1, ReferenceFrame(0,0,0))
        field_itp = interpolate(field, $xitp, nsy_Y)
        itpref = LinearInterpolation(($x,), reshape(field.e_SXY, 10), extrapolation_bc = 0.)
        @test all(field_itp.e_SXY .== vec(itpref.($xitp)))
        coef = Jolab.coefficient_general(itp_coef, field)
        (fieldl, fieldr) = lightinteraction(coef, field)
        @show fieldr.e_SXY
        @show field_itp.e_SXY
        stop
        # @test all(field_itp.e_SXY .== fieldr.e_SXY)
    end)
end

using Interpolations
x1 = range(0, 20, length = 11)
x2 = range(0, 20, length = 11)
y = x1 .+ x2'
itp = interpolate(y, BSpline(Cubic(Line(OnCell()))))
# itp = scale(itp, x1, x2)


x = (5.4, 5.4,)
wis = Interpolations.weightedindexes((Interpolations.value_weights,), itp, x)
abs(y[wis...] - itp(x...)) < 1E-8

A_x = 1.0:1.1:40.0
A_x = A_x.^5
A = rand(length(A_x), length(A_x))
nodes = (A_x, A_x)
itp = interpolate(nodes, A, BSpline(Linear()))
wis = Interpolations.weightedindexes((Interpolations.gradient_weights,), itp, x)
itp(x...)
A[wis...]

Interpolations.gradient_weights
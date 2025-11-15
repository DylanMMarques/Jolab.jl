using Jolab, Test

focal_len = 1e-3
lens = Lens(focal_len, .5, (Medium(1), Medium(1.0)), ReferenceFrame((0,0, focal_len), (0,0,0)))

@test_throws ArgumentError Lens(1E-3, .5, (Medium(1 + im), Medium(1.0)), ReferenceFrame((0,0,0), (0,0,0)))

nsx = range(-1, 1, length=10)
field = MonochromaticAngularSpectrum(Forward, nsx, nsx, rand(ComplexF64, 10, 10), 1500E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

field = MonochromaticAngularSpectrum(Forward, nsx, nsx, rand(ComplexF64, 10, 10), 1500E-9, Medium(1), ReferenceFrame((0,0,1), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

field = MonochromaticAngularSpectrum(Forward, nsx, nsx, rand(ComplexF64, 10, 10), 1500E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
(rfield, tfield) = light_interaction(lens, field)
@test tfield isa Jolab.MeshedSpatialBeam
@test rfield isa Jolab.MeshedAngularSpectrum


x = range(-1E-3, 1E-3, length = 10)
field = MonochromaticSpatialBeam(Forward, x, x, rand(ComplexF64, 10, 10), 1500E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)
field = MonochromaticSpatialBeam(Forward, x, x, rand(ComplexF64, 10, 10), 1500E-9, Medium(1), ReferenceFrame((0,0,1), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

field = MonochromaticSpatialBeam(Backward, x, x, rand(ComplexF64, 10, 10), 1500E-9, Medium(2), ReferenceFrame((0,0,0), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)
field = MonochromaticSpatialBeam(Backward, x, x, rand(ComplexF64, 10, 10), 1500E-9, Medium(1), ReferenceFrame((0,0,1), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

field = MonochromaticSpatialBeam(Backward, x, x, rand(ComplexF64, 10, 10), 1500E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

field = MonochromaticSpatialBeam(Forward, x, x, rand(ComplexF64, 10, 10), 1500E-9, Medium(1), ReferenceFrame((0,0,2focal_len), (0,0,0)))
@test_throws ArgumentError light_interaction(lens, field)

x = LinRange(-50E-6, 50E-6, 100)
field = MonochromaticSpatialBeam_gaussian(Forward, x, x, 10E-6, 1500E-9, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
(bfield, ffield) = light_interaction(lens, field)
@test bfield isa Jolab.MeshedSpatialBeam
@test ffield isa Jolab.MeshedAngularSpectrum
@test iszero(intensity(bfield))
@test isapprox(intensity(ffield), intensity(field); rtol = 1E-5)

x = LinRange(-100E-6, 100E-6, 100)
field = MonochromaticSpatialBeam_gaussian(Backward, x, x, 10E-6, 1500E-9, Medium(1), ReferenceFrame((0,0,2focal_len), (0,0,0)))
(bfield, ffield) = light_interaction(lens, field)
intensity(bfield)
@test bfield isa Jolab.MeshedAngularSpectrum
@test ffield isa Jolab.MeshedSpatialBeam
@test iszero(intensity(ffield))
@test isapprox(intensity(bfield), intensity(field); rtol = 1E-5)

sca = ScatteringMatrix(lens, field)
@test all(light_interaction(sca, field) .≈ light_interaction(lens, field))

x = LinRange(-.5, .5, 1000)
field = MonochromaticAngularSpectrum_gaussian(Backward, x, x, 10E-6, 1500E-9, Medium(1), ReferenceFrame((0,0,2focal_len), (0,0,0)))
(bfield, ffield) = light_interaction(lens, field)
@test bfield isa Jolab.MeshedSpatialBeam
@test ffield isa Jolab.MeshedAngularSpectrum
@test iszero(intensity(ffield))
@test isapprox(intensity(bfield), intensity(field); rtol = 1E-2) # Not sure if correct

sca = ScatteringMatrix(lens, field)
@test all(light_interaction(sca, field) .≈ light_interaction(lens, field))

using Enzyme
using FiniteDiff
import FiniteDiff: finite_difference_derivative

function field_after_lens(focal_len, nsx, nsy, ω, λ, medium)
    lens = Lens(focal_len, 0.5, (Medium(1), Medium(1.0)), ReferenceFrame((0,0, focal_len), (0,0,0)))
    field = MonochromaticAngularSpectrum_gaussian(Forward, nsx, nsy, ω, λ, medium, ReferenceFrame((0,0,0), (0,0,0)))
    (rfield, tfield) = light_interaction(lens, field)
    return maximum(abs, tfield.e)
end
 
at(ω) = field_after_lens(focal_len, nsx, nsx, ω, 1500E-9, Medium(1))
@test autodiff(Enzyme.Forward, at, Duplicated, Duplicated(10E-6, 1.0))[2] ≈
    finite_difference_derivative(x -> at(x), [10E-6], Val{:central}, Float64, relstep = 1E-9)[1]

at(f) = field_after_lens(f, nsx, nsx, 10E-6, 1500E-9, Medium(1))
@test autodiff(Enzyme.Forward, at, Duplicated, Duplicated(10E-6, 1.0))[2] ≈
    finite_difference_derivative(x -> at(x), [10E-6], Val{:central}, Float64, relstep = 1E-9)[1]

function field_after_lens(focal_len, x, y, ω, λ, medium)
    lens = Lens(focal_len, 0.5, (Medium(1), Medium(1.0)), ReferenceFrame((0,0, focal_len), (0,0,0)))
    field = MonochromaticSpatialBeam_gaussian(Forward, x, y, ω, λ, medium, ReferenceFrame((0,0,0), (0,0,0)))
    (rfield, tfield) = light_interaction(lens, field)
    return maximum(abs, tfield.e)
end
 
at(ω) = field_after_lens(focal_len, x, x, ω, 1500E-9, Medium(1))
@test autodiff(Enzyme.Forward, at, Duplicated, Duplicated(10E-6, 1.0))[2] ≈
    finite_difference_derivative(x -> at(x), [10E-6], Val{:central}, Float64, relstep = 1E-9)[1]

at(f) = field_after_lens(f, x, x, 10E-6, 1500E-9, Medium(1))
@test autodiff(Enzyme.Forward, at, Duplicated, Duplicated(10E-6, 1.0))[2] ≈
    finite_difference_derivative(x -> at(x), [10E-6], Val{:central}, Float64, relstep = 1E-9)[1]

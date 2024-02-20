using Test, Jolab


θ, ϕ = 0.015, π/4
pw = PlaneWaveScalar(Forward, 0, 0, 1, 1550E-9, Medium(1.0), ReferenceFrame((0,0,0), (0, 0, 0)))
pw2 = rotate_referenceframe(pw, (θ, 0, ϕ))
@test pw2.nsx ≈ sin(θ)*cos(ϕ)
@test pw2.nsy ≈ sin(θ)*sin(ϕ)
@test pw.e ≈ pw2.e

θ, ϕ = 0,0
p1 = (1, .5, .1)
p2 = (1.2, .6, .0)
pw = PlaneWaveScalar(Forward, .1, .35, 5.6 + 2.5im, 1550E-9, Medium(1.0), ReferenceFrame(p1, (θ, 0, ϕ)))
pw2 = Jolab.translate_referenceframe(pw, p2)
@test pw2.nsx ≈ .1 atol = 1E-15
@test pw2.nsy ≈ .35 atol = 1E-15
@test pw2.e ≈ pw.e * exp(im * 2π/pw.wavelength * sum((p2 .- p1) .* (pw.nsx, pw.nsy, √(pw.medium.n^2 - pw.nsx^2 - pw.nsy^2))))

nsx = range(-1, 1, length=2000) .+ zeros(2000)'
λ = 1550E-9
beam = Jolab.monochromatic_angularspectrum(Float64, nsx, nsx', (nsx .* nsx)' .+ 0im, λ, Medium(1.0 + im), ReferenceFrame((0,0,0), (0,0,0)))
ref2 = ReferenceFrame((100E-9,50E-9,10E-9), (0,0,0))

Jolab.translate_referenceframe(beam, ref2.origin)
Jolab.translate_referenceframe(beam, (1,1,1))


## Enzyme test

using Enzyme, Jolab, StaticArrays, FiniteDiff
import FiniteDiff: finite_difference_derivative
Enzyme.API.runtimeActivity!(true)


function f(in_x)
    nsx, nsy, e, λ, n, k, x, y, z = in_x
    pw = PlaneWaveScalar(nsx, nsy, e + im, λ * 1E-9, Medium(n + k*im), ReferenceFrame((x,y,z), (0,0,0)))
    pw2 = translate_referenceframe(pw, (2x, 2y, 2z))
    pw2.e
end

in_x = [0.1, 0.1, 1.0, 1550, 1, 0, 10E-9, 10E-9, 10E-9]
enz = enz = enz = enz = Enzyme.jacobian(Forward, f, in_x)
fin = FiniteDiff.finite_difference_jacobian(f, in_x, Val{:central}, ComplexF64, relstep = 1E-5)
@test all(isapprox.(enz, fin, atol = 1E-10))

function f(in_x)
    nsx, nsy, e, λ, n, x, y, z, θ, ϕ = in_x
    pw = PlaneWaveScalar(nsx, nsy, e + 0im, λ * 1E-9, Medium(n), ReferenceFrame((x * 1E-9, y * 1E-9, z * 1E-9), (0,0,0)))
    pw2 = translate_referenceframe(pw, (2x, 2y, 2z))
    pw3 = rotate_referenceframe(pw2, (θ, 0, ϕ))
    complex(pw3.nsx), complex(pw3.nsy), pw3.e
end
in_x = [0.1, 0.1, 1.0, 1550, 1, 10, 10, 10, π/4, π/8]
f(in_x)
enz = Enzyme.jacobian(Forward, f, in_x)
fin_nsx = FiniteDiff.finite_difference_jacobian(i -> f(i)[1], in_x, Val{:central}, ComplexF64, relstep = 1E-5)
fin_nsy = FiniteDiff.finite_difference_jacobian(i -> f(i)[2], in_x, Val{:central}, ComplexF64, relstep = 1E-5)
fin_e = FiniteDiff.finite_difference_jacobian(i -> f(i)[3], in_x, Val{:central}, ComplexF64, relstep = 1E-12)
@test all(isapprox.(first.(enz), fin_nsx, atol = 1E-10))
@test all(isapprox.(getindex.(enz,2), fin_nsy, atol = 1E-10))
@test all(isapprox.(last.(enz), fin_e, rtol = 1E-3))
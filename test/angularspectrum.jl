using Jolab, Test

nsx = range(-.05, .05, length = 101)
angspe = FieldAngularSpectrumScalar_gaussian(nsx, nsx, 50E-6, 1550E-9, 1., 1, ReferenceFrame(0,0,0,0,0))
int_i = intensity(angspe)
cart = Jolab.CartesianIndices(angspe)

θ, ϕ = 0.015, π/4
changereferenceframe!(angspe, ReferenceFrame(0,0,0,θ,ϕ))
(~,arg) = findmax(abs.(angspe.e_SXY))
@test isapprox(angspe.nsx_X[cart[arg][2]], -sin(θ)*cos(ϕ), atol = 1E-4)
@test isapprox(angspe.nsy_Y[cart[arg][3]], -sin(θ)*sin(ϕ), atol = 1E-4)
@test isapprox(intensity(angspe), int_i, atol = 1E-3)

θ, ϕ = 0, 0
changereferenceframe!(angspe, ReferenceFrame(0,0,0,θ,ϕ))
(~,arg) = findmax(abs.(angspe.e_SXY))
@test isapprox(angspe.nsx_X[cart[arg][2]], -sin(θ)*cos(ϕ), atol = 1E-4)
@test isapprox(angspe.nsy_Y[cart[arg][3]], -sin(θ)*sin(ϕ), atol = 1E-4)
@test isapprox(intensity(angspe), int_i, atol = 1E-3)

changereferenceframe!(angspe, ReferenceFrame(0,0,0,0,0))
(~,arg) = findmax(abs.(angspe.e_SXY))
nsx_max = angspe.nsx_X[cart[arg][2]]
nsy_max = angspe.nsy_Y[cart[arg][3]]

θ, ϕ = π/3, π/6
(x,y,z) = Jolab.rotatecoordinatesfromto(0,0,1., 0,0, θ, ϕ)
@test (x,y,z) == (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))
(x,y,z) = Jolab.rotatecoordinatesfromto(x,y,z, θ,ϕ, -2θ, -2ϕ)
@test (x,y,z) == (sin(-2θ) * cos(-2ϕ), sin(-2θ) * sin(-2ϕ), cos(-2θ))
(x,y,z) = Jolab.rotatecoordinatesfromto(x,y,z, -2θ,-2ϕ, 0,0)
@test all(isapprox.((x,y,z),(0,0,1), atol = 1E-16))

nsx = range(-.05, .05, length = 11)
angspe = FieldAngularSpectrumScalar_uniform(nsx, nsx, 1550E-9, 2., 1, ReferenceFrame(0,0,0,0,0))
Δx, Δy, Δz = 100E-9, 200E-9, 300E-9
angspe2 = changereferenceframe(angspe, ReferenceFrame(Δx, Δy, Δz, 0, 0))
k = 2π / angspe.λ
@test all(isapprox.(reshape(angspe.e_SXY ./ angspe2.e_SXY, 11, 11), exp.(-im .* k .* (angspe.nsx_X .* Δx .+ angspe.nsy_Y' .* Δy .+ .√(2^2 .- angspe.nsx_X.^2 .- angspe.nsy_Y'.^2) .* Δz)), rtol = 1E-15))

angspe = FieldAngularSpectrumScalar_uniform(nsx, nsx, 1550E-9, 2. + im, -1, ReferenceFrame(0,0,0,0,0))
Δx, Δy, Δz = 100E-9, 200E-9, 300E-9
angspe2 = changereferenceframe(angspe, ReferenceFrame(Δx, Δy, Δz, 0, 0))
k = 2π / angspe.λ
@test all(isapprox.(reshape(angspe.e_SXY ./ angspe2.e_SXY, 11, 11), exp.(-im .* k .* (angspe.nsx_X .* Δx .+ angspe.nsy_Y' .* Δy .- .√((2+im)^2 .- angspe.nsx_X.^2 .- angspe.nsy_Y'.^2) .* Δz)), rtol = 1E-15))

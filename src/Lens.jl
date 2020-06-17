mutable struct Lens{T} <: AbstractOpticalComponent{T}
	f::JolabFunction1D{T,T}
	na::T
	ref::ReferenceFrame{T}
	function Lens{A}(f, na, ref) where A
		return new{A}(f, na, ref)
	end
end
Lens(f, na, ref) = Lens{Float64}(f, na, ref)


function lightinteraction(lens::Lens, angspe::AbstractFieldAngularSpectrum{T}) where T
	lensref = (angspe.dir > 0) ? ref1(lens, angspe.λ) : ref2(lens, angspe.λ)
	# bestref = planelineintersection(lensref, angspe.ref)
	# checkposition(lensref, angspe.ref) || error("need to do this")
	# bestref = ReferenceFrame(lensref.x, lensref.y, lensref.z, angspe.ref.θ, angspe.ref.ϕ)
	angsperef = changereferential(angspe, lensref);

	if size(angsperef.e_SXY, 1) == 3
		(x_X, y_Y, e_SXY) = lightinteractionvectorial(lens, angsperef);
	else
		(x_X, y_Y, e_SXY) = lightinteractionscalar(lens, angsperef);
	end
	# (Δx, Δy, Δz) = rotatecoordinatesto(0, 0, 1, angsperef.ref.θ, angsperef.ref.ϕ)
	# (Δx, Δy, Δz) = rotatecoordinatesfrom(Δx, Δy, Δz, -lensref.θ, -lensref.ϕ)
	# x = lens.f(angspe.λ) * Δx
	# y = lens.f(angspe.λ) * Δy

	lensref2 = (angspe.dir > 0) ? ref2(lens, angsperef.λ) : ref1(lens, angsperef.λ)
	# (x, y, z) = rotatecoordinatesto(x, y, 0, lensref2.θ, lensref2.ϕ)

	# ref = ReferenceFrame(x + lensref2.x, y + lensref2.y, z + lensref2.z, lensref2.θ, lensref2.ϕ)

	angspe isa FieldAngularSpectrumSymmetric ? (return FieldSpaceSymmetric{T}(x_X, e_SXY, angspe.λ, angspe.n, angspe.dir, lensref2)) : (return FieldSpace{T}(x_X, y_Y, e_SXY, angspe.λ, angspe.n, angspe.dir, lensref2))
end

function lightinteractionvectorial(lens::Lens, angspe::FieldAngularSpectrum)
	@show("This must be remade. Do not trust result results")
	e_SXY = Array{Complex{Float64}, 3}(undef, size(angspe.e_SXY));
	(abs(imag(angspe.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
	sx_X = angspe.nsx_X;
	sy_Y = reshape(angspe.nsy_Y, 1,:)
	θ = acos.(.√(complex.(1 .- sx_X.^2 .- sy_Y.^2)))
	ϕ = atan.(sx_X, sy_Y)

	cosθ = cos.(θ); sinθ = sin.(θ); cosϕ = cos.(ϕ); sinϕ = sin.(ϕ);

	e_SXY[1,:,:] .= cosϕ .* (angspe.e_SXY[1,:,:] .* cosϕ .* cosθ .- angspe.e_SXY[3,:,:] .* sinθ) .+ angspe.e_SXY[2,:,:] .* (cosθ .- 1) .* cosϕ .* sinϕ .+ angspe.e_SXY[1,:,:] .* sinϕ.^2;
	e_SXY[2,:,:] .= angspe.e_SXY[2,:,:] .* cosϕ.^2 .+ angspe.e_SXY[1,:,:] .* (cosθ .- 1) .* cosϕ .* sinϕ .+ sinϕ .* (-angspe.e_SXY[3,:,:] .* sinθ .+ angspe.e_SXY[2,:,:] .* cosθ .* sinϕ);
	e_SXY[3,:,:] .= angspe.e_SXY[3,:,:] .* cosθ .+ sinθ .* (angspe.e_SXY[1,:,:] .* cosϕ .+ angspe.e_SXY[2,:,:] .* sinϕ);

	k = 2 .* π .* angspe.n ./ angspe.λ;
	e_SXY .= im ./ lens.f(angspe.λ) * k * 2 * π ./ .√(adddims(cosθ, (1,))) .* e_SXY;

	e_SXY .= e_SXY .* (adddims(real.(sinθ), (1,)) .< lens.na);

	x_X = sx_X * lens.f(angspe.λ);
	y_Y = angspe.nsy_Y * lens.f(angspe.λ);

	return (x_X, y_Y, e_SXY)
end

function lightinteractionscalar(lens::Lens, angspe::AbstractFieldAngularSpectrum{T}) where T
	(abs(imag(angspe.n)) < @tol) || error("To apply a lens the medium cannot absorb light")

	sx_X = angspe.nsx_X / real(angspe.n)
	sy_Y = angspe.nsy_Y / real(angspe.n)

	x_X = sx_X * lens.f(angspe.λ)
	y_Y = sy_Y * lens.f(angspe.λ)

	k = 2 * π * angspe.n / angspe.λ
	e_SXY = Array{Complex{T}, 3}(undef, 1, length(sx_X), length(sy_Y))
	@inbounds Threads.@threads for iY in eachindex(sy_Y)
		for iX in eachindex(sx_X)
			cosθ2 = 1 - sx_X[iX]^2 - sy_Y[iY]^2 # cosθ2 is cosθ^2
			if cosθ2 < 0 # Evasnecent waves
				e_SXY[1,iX,iY] = zero(Complex{T})
			else # Plane waves
				(1 - cosθ2 > lens.na^2) && (e_SXY[1,iX,iY] = zero(Complex{T}); continue)
				e_SXY[1,iX,iY] = im / lens.f(angspe.λ) * k * 2π / cosθ2^(1/4) * angspe.e_SXY[1,iX,iY]
			end
		end
	end
	return (x_X, y_Y, e_SXY)
end

function lightinteraction(lens::Lens, space::AbstractFieldSpace{T}) where T
	lensref = space.dir > 0 ? ref1(lens, space.λ) : ref2(lens,space.λ)
	# bestref = planelineintersection(lensref, space.ref)
	spaceref = changereferential(space, lensref);

	if size(spaceref.e_SXY, 1) == 3
		(nsx_X, nsy_Y, e_SXY) = lightinteractionvectorial(lens, spaceref);
	else
		(nsx_X, nsy_Y, e_SXY) = lightinteractionscalar(lens, spaceref);
	end

	lensref2 = space.dir > 0 ? ref2(lens, space.λ) : ref1(lens,space.λ)

	space isa FieldSpaceSymmetric ? (return FieldAngularSpectrumSymmetric{T}(nsx_X, e_SXY, space.λ, space.n, space.dir, ref)) : (return FieldAngularSpectrum{T}(nsx_X, nsy_Y, e_SXY, space.λ, space.n, space.dir, lensref2))
	return angspe;
end

function lightinteractionscalar(lens::Lens, space::AbstractFieldSpace{T}) where T
	(abs(imag(space.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
	f = lens.f(space.λ)
	sx_X = -space.x_X / f
	sy_Y = -space.y_Y / f

	k = 2 * π * space.n / space.λ
	e_SXY = Array{Complex{T},3}(undef, 1, length(space.x_X), length(space.y_Y))
	@inbounds Threads.@threads for iY in eachindex(space.y_Y)
		for iX in eachindex(space.x_X)
			cosθ2 = (f^2 - space.x_X[iX]^2 - space.y_Y[iY]^2)
			if cosθ2 < 0
				e_SXY[1,iX,iY] = zero(Complex{T})
			else
				cosθ2 = √(cosθ2) / f
				(1 - cosθ2 > lens.na^2) || (e_SXY[1,iX,iY] = zero(Complex{T}); continue)
				e_SXY[1,iX,iY] = -im * f / k / 2 / π * √(cosθ2) * space.e_SXY[1,iX,iY];
			end
		end
	end
	sx_X *= real(space.n) # is nsx
	sy_Y *= real(space.n) # is nsy
	return (sx_X, sy_Y, e_SXY)
end

function lightinteractionvectorial(lens::Lens, space::FieldSpace)
	@show "This must be remade. Do not trust results"
	e_SXY = Array{Complex{Float64}, 3}(undef, size(space.e_SXY));
	x_X = space.x_X;
	y_Y = reshape(space.y_Y, 1, :);
	θ = acos.(.√(complex.(lens.f(space.λ).^2 .- x_X.^2 .- y_Y.^2)) ./ lens.f(space.λ));
	ϕ = π .+ atan.(y_Y, x_X);

	cosθ = cos.(θ); sinθ = sin.(θ); cosϕ = cos.(ϕ); sinϕ = sin.(ϕ);

	e_SXY[1,:,:] .= (cosθ .* cosϕ .* cosϕ .+ sinϕ.^2) .* space.e_SXY[1,:,:] .+ cosϕ .* sinθ .* space.e_SXY[3,:,:] .+ (cosθ .- 1) .* cosϕ .* sinϕ .* space.e_SXY[2,:,:];
	e_SXY[2,:,:] .= (cosϕ.^2 .+ cosθ .* sinϕ.^2) .* space.e_SXY[2,:,:] + (cosθ .- 1) .* cosϕ .* sinϕ .* space.e_SXY[1,:,:] .+ sinθ .* sinϕ .* space.e_SXY[3,:,:];
	e_SXY[3,:,:] .= .- sinθ .* cosϕ .* space.e_SXY[1,:,:] .- sinθ .* sinϕ .* space.e_SXY[2,:,:] .+ cosθ .* space.e_SXY[3,:,:];

	k = 2 .* π .* space.n ./ space.λ;
	e_SXY .= -im .* lens.f(space.λ) ./ k ./ 2 ./ π .* .√(adddims(cosθ, (1,))) .* e_SXY;
	e_SXY .= e_SXY .* (adddims(real.(sinθ), (1,)) .< lens.na);

	sx_X = .-space.x_X ./ lens.f(space.λ);
	sy_Y = .-space.y_Y ./ lens.f(space.λ);

	return (sx_X, sy_Y, e_SXY)
end

function ref1(lens::Lens, λ::Real)
	x = lens.ref.x - lens.f(λ) * sin(lens.ref.θ) * cos(lens.ref.ϕ)
	y = lens.ref.y - lens.f(λ) * sin(lens.ref.θ) * sin(lens.ref.ϕ)
	z = lens.ref.z - lens.f(λ) * cos(lens.ref.θ)
	return ReferenceFrame(x, y, z, lens.ref.θ, lens.ref.ϕ)
end

function ref2(lens::Lens, λ::Real)
	x = lens.ref.x + lens.f(λ) * sin(lens.ref.θ) * cos(lens.ref.ϕ)
	y = lens.ref.y + lens.f(λ) * sin(lens.ref.θ) * sin(lens.ref.ϕ)
	z = lens.ref.z + lens.f(λ) * cos(lens.ref.θ)
	return ReferenceFrame(x, y, z, lens.ref.θ, lens.ref.ϕ)
end

mutable struct Lens{T} <: AbstractOpticalComponent{T}
	f::JolabFunction1D{T,T}
	na::T
	ref::ReferenceFrame{T}
	function Lens{A}(f, na, ref) where A
		return new{A}(f, na, ref)
	end
end
Lens(f, na, ref) = Lens{Float64}(f, na, ref)

function ref1(lens::Lens{T}, λ::Real) where T
	x = lens.ref.x - lens.f(λ) * sin(lens.ref.θ) * cos(lens.ref.ϕ)
	y = lens.ref.y - lens.f(λ) * sin(lens.ref.θ) * sin(lens.ref.ϕ)
	z = lens.ref.z - lens.f(λ) * cos(lens.ref.θ)
	return ReferenceFrame{T}(x, y, z, lens.ref.θ, lens.ref.ϕ)
end

function ref2(lens::Lens{T}, λ::Real) where T
	x = lens.ref.x + lens.f(λ) * sin(lens.ref.θ) * cos(lens.ref.ϕ)
	y = lens.ref.y + lens.f(λ) * sin(lens.ref.θ) * sin(lens.ref.ϕ)
	z = lens.ref.z + lens.f(λ) * cos(lens.ref.θ)
	return ReferenceFrame{T}(x, y, z, lens.ref.θ, lens.ref.ϕ)
end

function coefficient_general(lens::Lens{T}, fieldi::FieldAngularSpectrum) where T
	(abs(imag(fieldi.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
	checkorientation(fieldi.ref, lens.ref) || error("Cannot calculate coefficient with change of referential orientation. Use lightinteraction(lens,...) instead")

	ref = (fieldi.dir > 0 ? ref1(lens, fieldi.λ) : ref2(lens, fieldi.λ))
	needProp = checkposition(fieldi.ref, ref)
	needProp && (propM = propagationmatrix(fieldi, ref))
	fieldi.dir > 0 || conj!(propM.diag)

	sizeXY = length(fieldi.nsx_X) * length(fieldi.nsy_Y)
	sx_X = fieldi.nsx_X / real(fieldi.n)
	sy_Y = fieldi.nsy_Y / real(fieldi.n)

	f = lens.f(fieldi.λ)
	x_X = sx_X * f
	y_Y = sy_Y * f

	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r12 = Diagonal(zeros(Complex{T}, sizeXY))
	t21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	k = 2 * π * fieldi.n / fieldi.λ
	i = 1
	@inbounds for iY in eachindex(sy_Y)
		for iX in eachindex(sx_X)
			cosθ² = 1 - sx_X[iX]^2 - sy_Y[iY]^2 # cosθ2 is cosθ^2
			if cosθ² < 0 # Evasnecent waves
				t12.diag[i] = zero(Complex{T})
				t21.diag[i] = zero(Complex{T})
			else # Plane waves
				if (1 - cosθ² > lens.na^2)
					t12.diag[i] = zero(Complex{T})
					t21.diag[i] = zero(Complex{T})
				else
					t12.diag[i] = im / f * k * 2π / cosθ²^(1/4)
					t21.diag[i] = 1 / t12.diag[i]
				end
			end
			i += 1
		end
	end
	if fieldi.dir > 0
		if needProp
			rmul!(t12, propM)
			lmul!(propM, t21)
		end
		fieldl = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, -1, fieldi.ref)
		fieldr = FieldSpace{T}(x_X, y_Y, fieldi.e_SXY, fieldi.λ, fieldi.n, 1, ref2(lens, fieldi.λ))
	else
		aux = t12
		t12 = t21
		t21 = aux
		if needProp
			lmul!(propM, t12)
			rmul!(t21, propM)
		end
		fieldl = FieldSpace{T}(x_X, y_Y, fieldi.e_SXY, fieldi.λ, fieldi.n, -1, ref1(lens, fieldi.λ))
		fieldr = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, 1, fieldi.ref)
	end
	return ScaterringMatrix{T, Diagonal{Complex{T},Vector{Complex{T}}}, typeof(fieldl), typeof(fieldr)}(r12, t12, r12, t21, fieldl, fieldr)
end

function coefficient_general(lens::Lens{T}, field::FieldSpace{A,X}) where {T,A,X}
	(abs(imag(field.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
	checkorientation(field.ref, lens.ref) || error("Cannot calculate coefficient with change of referential orientation. Use lightinteraction(lens,...) instead")

	ref = (field.dir > 0 ? ref1(lens, field.λ) : ref2(lens, field.λ))
	fieldi = changereferential(field, ref)

	sizeXY = length(fieldi.x_X) * length(fieldi.y_Y)
	f = lens.f(fieldi.λ)
	sx_X = -fieldi.x_X / f
	sy_Y = -fieldi.y_Y / f

	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r12 = Diagonal(zeros(Complex{T}, sizeXY))
	t21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	k = 2 * π * fieldi.n / fieldi.λ
	i = 1
	@inbounds for iY in eachindex(sy_Y)
		for iX in eachindex(sx_X)
			cosθ² = 1 - sx_X[iX]^2 - sy_Y[iY]^2 # cosθ2 is cosθ^2
			if cosθ² < 0 # Evasnecent waves
				t12.diag[i] = zero(Complex{T})
				t21.diag[i] = zero(Complex{T})
			else # Plane waves
				if (1 - cosθ² > lens.na^2)
					t12.diag[i] = zero(Complex{T})
					t21.diag[i] = zero(Complex{T})
				else
					t21.diag[i] = im / f * k * 2π / cosθ²^(1/4)
					t12.diag[i] = 1 / t21.diag[i]
				end
			end
			i += 1
		end
	end
	if fieldi.dir > 0
		fieldl = FieldSpace{T}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, -1, field.ref)
		fieldr = FieldAngularSpectrum{T}(real(fieldi.n) * sx_X, real(fieldi.n) * sy_Y, fieldi.e_SXY, fieldi.λ, fieldi.n, 1, ref2(lens, fieldi.λ))
	else
		aux = t12
		t12 = t21
		t21 = aux
		fieldl = FieldAngularSpectrum{T}(real(fieldi.n) * sx_X, real(fieldi.n) * sy_Y, fieldi.e_SXY, fieldi.λ, fieldi.n, -1, ref1(lens, fieldi.λ))
		fieldr = FieldSpace{T}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, 1, fieldi.ref)
	end
	return ScaterringMatrix{T, Diagonal{Complex{T},Vector{Complex{T}}}, typeof(fieldl), typeof(fieldr)}(r12, t12, r12, t21, fieldl, fieldr)
end

@inline coefficient_specific(lens::Lens, field::AbstractFieldMonochromatic) = coefficient_general(lens, field)

# function lightinteractionvectorial(lens::Lens, angspe::FieldAngularSpectrum)
# 	@show("This must be remade. Do not trust result results")
# 	e_SXY = Array{Complex{Float64}, 3}(undef, size(angspe.e_SXY));
# 	(abs(imag(angspe.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
# 	sx_X = angspe.nsx_X;
# 	sy_Y = reshape(angspe.nsy_Y, 1,:)
# 	θ = acos.(.√(complex.(1 .- sx_X.^2 .- sy_Y.^2)))
# 	ϕ = atan.(sx_X, sy_Y)
#
# 	cosθ = cos.(θ); sinθ = sin.(θ); cosϕ = cos.(ϕ); sinϕ = sin.(ϕ);
#
# 	e_SXY[1,:,:] .= cosϕ .* (angspe.e_SXY[1,:,:] .* cosϕ .* cosθ .- angspe.e_SXY[3,:,:] .* sinθ) .+ angspe.e_SXY[2,:,:] .* (cosθ .- 1) .* cosϕ .* sinϕ .+ angspe.e_SXY[1,:,:] .* sinϕ.^2;
# 	e_SXY[2,:,:] .= angspe.e_SXY[2,:,:] .* cosϕ.^2 .+ angspe.e_SXY[1,:,:] .* (cosθ .- 1) .* cosϕ .* sinϕ .+ sinϕ .* (-angspe.e_SXY[3,:,:] .* sinθ .+ angspe.e_SXY[2,:,:] .* cosθ .* sinϕ);
# 	e_SXY[3,:,:] .= angspe.e_SXY[3,:,:] .* cosθ .+ sinθ .* (angspe.e_SXY[1,:,:] .* cosϕ .+ angspe.e_SXY[2,:,:] .* sinϕ);
#
# 	k = 2 .* π .* angspe.n ./ angspe.λ;
# 	e_SXY .= im ./ lens.f(angspe.λ) * k * 2 * π ./ .√(adddims(cosθ, (1,))) .* e_SXY;
#
# 	e_SXY .= e_SXY .* (adddims(real.(sinθ), (1,)) .< lens.na);
#
# 	x_X = sx_X * lens.f(angspe.λ);
# 	y_Y = angspe.nsy_Y * lens.f(angspe.λ);
#
# 	return (x_X, y_Y, e_SXY)
# end
#
# function lightinteractionscalar(lens::Lens, angspe::AbstractFieldAngularSpectrum{T}) where T
# 	(abs(imag(angspe.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
#
# 	sx_X = angspe.nsx_X / real(angspe.n)
# 	sy_Y = angspe.nsy_Y / real(angspe.n)
#
# 	x_X = sx_X * lens.f(angspe.λ)
# 	y_Y = sy_Y * lens.f(angspe.λ)
#
# 	k = 2 * π * angspe.n / angspe.λ
# 	e_SXY = Array{Complex{T}, 3}(undef, 1, length(sx_X), length(sy_Y))
# 	@inbounds Threads.@threads for iY in eachindex(sy_Y)
# 		for iX in eachindex(sx_X)
# 			cosθ2 = 1 - sx_X[iX]^2 - sy_Y[iY]^2 # cosθ2 is cosθ^2
# 			if cosθ2 < 0 # Evasnecent waves
# 				e_SXY[1,iX,iY] = zero(Complex{T})
# 			else # Plane waves
# 				(1 - cosθ2 > lens.na^2) && (e_SXY[1,iX,iY] = zero(Complex{T}); continue)
# 				e_SXY[1,iX,iY] = im / lens.f(angspe.λ) * k * 2π / cosθ2^(1/4) * angspe.e_SXY[1,iX,iY]
# 			end
# 		end
# 	end
# 	return (x_X, y_Y, e_SXY)
# end
#
# function lightinteractionscalar(lens::Lens, space::AbstractFieldSpace{T}) where T
# 	(abs(imag(space.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
# 	f = lens.f(space.λ)
# 	sx_X = -space.x_X / f
# 	sy_Y = -space.y_Y / f
#
# 	k = 2 * π * space.n / space.λ
# 	e_SXY = Array{Complex{T},3}(undef, 1, length(space.x_X), length(space.y_Y))
# 	@inbounds Threads.@threads for iY in eachindex(space.y_Y)
# 		for iX in eachindex(space.x_X)
# 			cosθ2 = (f^2 - space.x_X[iX]^2 - space.y_Y[iY]^2)
# 			if cosθ2 < 0
# 				e_SXY[1,iX,iY] = zero(Complex{T})
# 			else
# 				cosθ2 = √(cosθ2) / f
# 				(1 - cosθ2 > lens.na^2) && (e_SXY[1,iX,iY] = zero(Complex{T}); continue)
# 				e_SXY[1,iX,iY] = -im * f / k / 2 / π * √(cosθ2) * space.e_SXY[1,iX,iY];
# 			end
# 		end
# 	end
# 	sx_X *= real(space.n) # is nsx
# 	sy_Y *= real(space.n) # is nsy
# 	return (sx_X, sy_Y, e_SXY)
# end
#
# function lightinteractionvectorial(lens::Lens, space::FieldSpace)
# 	@show "This must be remade. Do not trust results"
# 	e_SXY = Array{Complex{Float64}, 3}(undef, size(space.e_SXY));
# 	x_X = space.x_X;
# 	y_Y = reshape(space.y_Y, 1, :);
# 	θ = acos.(.√(complex.(lens.f(space.λ).^2 .- x_X.^2 .- y_Y.^2)) ./ lens.f(space.λ));
# 	ϕ = π .+ atan.(y_Y, x_X);
#
# 	cosθ = cos.(θ); sinθ = sin.(θ); cosϕ = cos.(ϕ); sinϕ = sin.(ϕ);
#
# 	e_SXY[1,:,:] .= (cosθ .* cosϕ .* cosϕ .+ sinϕ.^2) .* space.e_SXY[1,:,:] .+ cosϕ .* sinθ .* space.e_SXY[3,:,:] .+ (cosθ .- 1) .* cosϕ .* sinϕ .* space.e_SXY[2,:,:];
# 	e_SXY[2,:,:] .= (cosϕ.^2 .+ cosθ .* sinϕ.^2) .* space.e_SXY[2,:,:] + (cosθ .- 1) .* cosϕ .* sinϕ .* space.e_SXY[1,:,:] .+ sinθ .* sinϕ .* space.e_SXY[3,:,:];
# 	e_SXY[3,:,:] .= .- sinθ .* cosϕ .* space.e_SXY[1,:,:] .- sinθ .* sinϕ .* space.e_SXY[2,:,:] .+ cosθ .* space.e_SXY[3,:,:];
#
# 	k = 2 .* π .* space.n ./ space.λ;
# 	e_SXY .= -im .* lens.f(space.λ) ./ k ./ 2 ./ π .* .√(adddims(cosθ, (1,))) .* e_SXY;
# 	e_SXY .= e_SXY .* (adddims(real.(sinθ), (1,)) .< lens.na);
#
# 	sx_X = .-space.x_X ./ lens.f(space.λ);
# 	sy_Y = .-space.y_Y ./ lens.f(space.λ);
#
# 	return (sx_X, sy_Y, e_SXY)
# end

mutable struct FieldSpace{T,X<:AbstractVector{T}} <: AbstractFieldSpace{T}
	x_X::X
	y_Y::X
	e_SXY::Array{Complex{T}, 3}
	λ::T
	n::Complex{T}
	dir::Int8
	ref::ReferenceFrame{T}
	function FieldSpace{T}(x_X::X, y_Y::Y, e_SXY, λ, n, dir, ref) where {T,X,Y}
		M = promote_type(X,Y)

		length(x_X) != size(e_SXY, 2) ? error("The length of x_X must be the same as the size of e_SXY in the second dimension.") : nothing;
		length(y_Y) != size(e_SXY, 3) ? error("The length of y_Y must be the same as the size of e_SXY in the third dimension.") : nothing;
		size(e_SXY, 1) != 1 && size(e_SXY, 1) != 3 ? error("The size of e_SXY must be 1 (scallar field) or 3 (vectorial field).") : nothing;
		return new{T,M}(x_X, y_Y, e_SXY, λ, n, dir >= 0 ? 1 : -1, ref);
	end
end
# FieldSpace(x_X, y_Y, e_SXY, λ, n, dir, ref) = FieldSpace{Float64}(x_X, y_Y, e_SXY, λ, n, dir, ref)

function FieldSpace_uniform(x_X, y_Y, λ, n, dir, ref)
	T = Float64
	e_SXY = ones(Complex{T}, 1, length(x_X), length(y_Y))
	space = FieldSpace{T}(x_X, y_Y, e_SXY, λ, n, dir, ref);
	int = √(intensity(space))
	space.e_SXY ./= int
	return space
end

function FieldSpace_gaussian(x_X, y_Y, ω, λ, n, dir, ref)
	T = Float64
	norm = 	√(8 / π / ω^2)
	e_SXY = Array{Complex{T}, 3}(undef, 1, length(x_X), length(y_Y))
	@inbounds @simd for iY in eachindex(y_Y)
		for iX in eachindex(x_X)
			e_SXY[1,iX,iY] = norm * exp(- 4 / ω^2 * (x_X[iX]^2 + y_Y[iY]^2));
		end
	end
	return FieldSpace{T}(vec(x_X), vec(y_Y), e_SXY, λ, n, dir, ref);
end

function FieldAngularSpectrum_fromspace(fieldspace::FieldSpace{T}, nsx_X, nsy_Y) where T
	e_SXY = fourriertransform(fieldspace.x_X, fieldspace.y_Y, fieldspace.e_SXY, fieldspace.λ, fieldspace.n, nsx_X, nsy_Y);
	return FieldAngularSpectrum{T}(nsx_X, nsy_Y, e_SXY, fieldspace.λ, fieldspace.n, fieldspace.dir, fieldspace.ref);
end

function FieldAngularSpectrum_fromspacefft(fieldspace::FieldSpace{T,X}; padding=0::Integer) where {T, X<:AbstractRange}
	(nsx_X, nsy_Y, e_SXY) = fourriertransformfft(fieldspace.x_X, fieldspace.y_Y, fieldspace.e_SXY, fieldspace.λ, padding=padding)
	return FieldAngularSpectrum{T}(nsx_X, nsy_Y, e_SXY, fieldspace.λ, fieldspace.n, fieldspace.dir, fieldspace.ref);
end

function changereferential!(fieldspace::FieldSpace, refnew::ReferenceFrame)
	if fieldspace.ref != refnew
		(fieldspace.x_X, fieldspace.y_Y) = translatereferential!(fieldspace.x_X, fieldspace.y_Y, fieldspace.ref, refnew);
		fieldspace.ref = refnew;
	end
end

function changereferential(fieldspace::FieldSpace, refnew::ReferenceFrame)
	tmp_field = deepcopy(fieldspace)
	changereferential!(tmp_field, refnew)
	return tmp_field
end

function translatereferential!(x_X::AbstractVector{<:Real}, y_Y::AbstractVector{<:Real}, refold::ReferenceFrame, refnew::ReferenceFrame)
	(checkorientation(refold, refnew) && checkinplane(refold, refnew)) || error("Both referential must be in plane")
	refΔx = refnew.x - refold.x;
	refΔy = refnew.y - refold.y;
	refΔz = refnew.z - refold.z;

	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, refold.θ, refold.ϕ);
	x_X .+= refΔx
	y_Y .+= refΔy
	return (x_X, y_Y)
end

function translatereferential!(x_X::AbstractRange{<:Real}, y_Y::AbstractRange{<:Real}, refold::ReferenceFrame, refnew::ReferenceFrame)
	(checkorientation(refold, refnew) && checkinplane(refold, refnew)) || error("Both referential must be in plane")
	refΔx = refnew.x - refold.x;
	refΔy = refnew.y - refold.y;
	refΔz = refnew.z - refold.z;

	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, refold.θ, refold.ϕ);
	x_X = x_X .+ refΔx
	y_Y = y_Y .+ refΔy
	return (x_X, y_Y)
end

function intensity(fieldspace::FieldSpace)
	size(fieldspace.e_SXY, 1) == 3 ? (e_SXY = dotdim(fieldspace.e_SXY, fieldspace.e_SXY, 1)) : e_SXY = abs2.(fieldspace.e_SXY[1,:,:])
	return real(fieldspace.n) * ∫∫(e_SXY, fieldspace.x_X, fieldspace.y_Y);
end

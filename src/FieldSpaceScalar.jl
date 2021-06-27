mutable struct FieldSpaceScalar{T,D, X<:AbstractVector{T}, Y<:AbstractVector{Complex{T}}} <: AbstractFieldSpace{T,D}
	x_X::X
	y_Y::X
	e_SXY::Y
	λ::T
	n::Complex{T}
	ref::ReferenceFrame{T}
	function FieldSpaceScalar{T,D,X,Y}(x_X, y_Y, e_SXY::Y, λ, n, ref) where {T,D,X,Y}
		length(x_X) * length(y_Y) == length(e_SXY) || error("length(x_X) * length(y_Y) must be equal to length(e_SXY)")
		return new{T,D,X,Y}(x_X, y_Y, e_SXY, λ, n, ref);
	end
	function FieldSpaceScalar{T,D}(x_X::X1, y_Y::X2, e_SXY::Y, λ, n, ref) where {T,D,X1,X2,Y}
		length(x_X) * length(y_Y) == length(e_SXY) || error("length(x_X) * length(y_Y) must be equal to length(e_SXY)")
		X = promote_type(X1,X2)
		return new{T,D,X,Y}(x_X, y_Y, e_SXY, λ, n, ref);
	end
	function FieldSpaceScalar{T,D,X}(x_X, y_Y, e_SXY::Y, λ, n, ref) where {T,D,X,Y}
		length(x_X) * length(y_Y) == length(e_SXY) || error("length(x_X) * length(y_Y) must be equal to length(e_SXY)")
		return new{T,D,X,Y}(x_X, y_Y, e_SXY, λ, n, ref);
	end
end
# FieldSpace(x_X, y_Y, e_SXY, λ, n, dir, ref) = FieldSpace{Float64}(x_X, y_Y, e_SXY, λ, n, dir, ref)

function FieldSpaceScalar_uniform(::Type{T}, x_X, y_Y, λ, n, dir, ref) where T
	e_SXY = ones(Complex{T}, length(x_X) * length(y_Y))
	space = FieldSpaceScalar{T,dir}(x_X, y_Y, e_SXY, λ, n, ref);
	space.e_SXY ./= √(intensity(space))
	return space
end
FieldSpaceScalar_uniform(arg...) = FieldSpaceScalar_uniform(Float64, arg...)

function FieldSpaceScalar_gaussian(::Type{T}, x_X, y_Y, ω, λ, n, dir, ref) where T
	e_SXY = Vector{Complex{T}}(undef, length(x_X) * length(y_Y))
	space = FieldSpaceScalar{T,dir}(x_X, y_Y, e_SXY, λ, n, ref);
	norm = 	√(T(8 / π / ω^2)) / √T(real(n))
	cart = CartesianIndices(space)
	@inbounds @simd for i in iterator_index(space)
		space.e_SXY[i] = norm * exp(- 4 / T(ω)^2 * (x_X[cart[i][2]]^2 + y_Y[cart[i][3]]^2));
	end
	return space
end
FieldSpaceScalar_gaussian(arg...) = FieldSpaceScalar_gaussian(Float64, arg...)

function changereferenceframe!(fieldspace::AbstractFieldSpace, refnew::ReferenceFrame)
	(checkorientation(fieldspace.ref, refnew) && checkinplane(fieldspace.ref, refnew)) || error("Both referenceframe must be in plane")
	translatereferenceframe!(fieldspace, refnew)
end

function translatereferenceframe!(space::AbstractFieldSpace, refnew::ReferenceFrame)
	(refΔx, refΔy, refΔz) = (refnew.x - space.ref.x, refnew.y - space.ref.y, refnew.z - space.ref.z)
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, space.ref.θ, space.ref.ϕ);
	(space.x_X, space.y_Y) = (space.x_X .+ refΔx, space.y_Y .+ refΔy)
	(space.ref.x, space.ref.y, space.ref.z) = (refnew.x, refnew.y, refnew.z)
end

function intensity(space::FieldSpaceScalar{T}) where T
	int = zero(T)
	cart = CartesianIndices(space)
	@inbounds @simd for i in iterator_index(space)
		Δx = Δvector(space.x_X, cart[i][2])
		Δy = Δvector(space.y_Y, cart[i][3])
		int += abs2(space.e_SXY[i]) * Δx * Δy
	end
	return real(space.n) * int
end

function propagationmatrix(fieldl::L, fieldr::L) where {L <: AbstractFieldSpace}
	error("A field in space cannot be propagated transform to angular spectrum")
end

function samedefinitions(fieldl::L, fieldr::R) where {L<:FieldSpaceScalar, R<:FieldSpaceScalar}
	isapprox(fieldl.x_X, fieldr.x_X, atol = @tol) || error("x_X are different")
	isapprox(fieldl.y_Y, fieldr.y_Y, atol = @tol) || error("y_Y are different")
	isapprox(fieldl.n, fieldr.n, atol = @tol) || error("n are different")
	isapprox(fieldl.λ, fieldr.λ, atol = @tol) || error("λ are different")
	fieldl.ref == fieldr.ref || error("reference frames are different")
	return true
end

function add!(fielda::FieldSpaceScalar{T}, fieldb::FieldSpaceScalar{T}) where T
	samedefinitions(fielda, fieldb) || error("Cannot sum the fields. Different definitions")
	@inbounds @simd for i in iterator_index(fielda)
		fielda.e_SXY[i] += fieldb.e_SXY[i]
	end
end

function Base.:copy(field::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y}
	return FieldSpaceScalar{T,D,X,Y}(deepcopy(field.x_X), deepcopy(field.y_Y), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

function copy_differentD(field::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y}
	return FieldSpaceScalar{T,-D,X,Y}(deepcopy(field.x_X), deepcopy(field.y_Y), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

CartesianIndices(field::FieldSpaceScalar) = Base.CartesianIndices((1, length(field.x_X), length(field.y_Y)))
LinearIndices(field::FieldSpaceScalar) = Base.LinearIndices((1, length(field.x_X), length(field.y_Y)))

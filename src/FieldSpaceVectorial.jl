mutable struct FieldSpaceVectorial{T,D, X<:AbstractVector{T}, Y<:AbstractVector{Complex{T}}} <: AbstractFieldSpace{T,D}
	x_X::X
	y_Y::X
	e_SXY::Y
	λ::T
	n::Complex{T}
	ref::ReferenceFrame{T}
	function FieldSpaceVectorial{T,D,X,Y}(x_X, y_Y, e_SXY::Y, λ, n, ref) where {T,D,X,Y}
		length(x_X) * length(y_Y) * 3 == length(e_SXY) || error("length(x_X) * length(y_Y) * 3 must be equal to length(e_SXY)")
		return new{T,D,X,Y}(x_X, y_Y, e_SXY, λ, n, ref);
	end
	function FieldSpaceVectorial{T,D}(x_X::X1, y_Y::X2, e_SXY::Y, λ, n, ref) where {T,D,X1,X2,Y}
		length(x_X) * length(y_Y) * 3 == length(e_SXY) || error("length(x_X) * length(y_Y) * 3 must be equal to length(e_SXY)")
		X = promote_type(X1,X2)
		return new{T,D,X,Y}(x_X, y_Y, e_SXY, λ, n, ref);
	end
	function FieldSpaceVectorial{T,D,X}(x_X, y_Y, e_SXY::Y, λ, n, ref) where {T,D,X,Y}
		length(x_X) * length(y_Y) * 3 == length(e_SXY) || error("length(x_X) * length(y_Y) * 3 must be equal to length(e_SXY)")
		return new{T,D,X,Y}(x_X, y_Y, e_SXY, λ, n, ref);
	end
end
# FieldSpace(x_X, y_Y, e_SXY, λ, n, dir, ref) = FieldSpace{Float64}(x_X, y_Y, e_SXY, λ, n, dir, ref)

function FieldSpaceVectorial_uniform(::Type{T}, x_X, y_Y, λ, n, dir, ref) where T
	e_SXY = ones(Complex{T}, 3 * length(x_X) * length(y_Y))
	space = FieldSpaceVectorial{T,dir}(x_X, y_Y, e_SXY, λ, n, ref);
	space.e_SXY ./= √(intensity(space))
	return space
end
FieldSpaceVectorial_uniform(arg...) = FieldSpaceVectorial_uniform(Float64, arg...)

function FieldSpaceVectorial_gaussian(::Type{T}, x_X, y_Y, ω, λ, n, dir, ref) where T
	e_SXY = Vector{Complex{T}}(undef, 3 * length(x_X) * length(y_Y))
	space = FieldSpaceVectorial{T,dir}(x_X, y_Y, e_SXY, λ, n, ref);
	norm = 	√(T(8 / π / ω^2))
	cart = CartesianIndices(space)
	@inbounds @simd for i in iterator_index(space)
		space.e_SXY[i] = norm * exp(- 4 / T(ω)^2 * (x_X[cart[i][2]]^2 + y_Y[cart[i][3]]^2));
	end
	return space
end
FieldSpaceVectorial_gaussian(arg...) = FieldSpaceVectorial_gaussian(Float64, arg...)

function intensity(space::FieldSpaceVectorial{T}) where T
	int = zero(T)
	cart = CartesianIndices(space)
	for i in iterator_index(space)
		Δx = Δvector(space.x_X, cart[i][2])
		Δy = Δvector(space.y_Y, cart[i][3])
		int += abs2(space.e_SXY[i]) * Δx * Δy
	end
	return real(space.n) * int
end

function samedefinitions(fieldl::L, fieldr::R) where {L<:FieldSpaceVectorial, R<:FieldSpaceVectorial}
	isapprox(fieldl.x_X, fieldr.x_X, atol = @tol) || error("x_X are different")
	isapprox(fieldl.y_Y, fieldr.y_Y, atol = @tol) || error("y_Y are different")
	isapprox(fieldl.n, fieldr.n, atol = @tol) || error("n are different")
	isapprox(fieldl.λ, fieldr.λ, atol = @tol) || error("λ are different")
	fieldl.ref == fieldr.ref || error("reference frames are different")
	return true
end

function add!(fielda::FieldSpaceVectorial{T}, fieldb::FieldSpaceVectorial{T}) where T
	samedefinitions(fielda, fieldb) || error("Cannot sum the fields. Different definitions")
	@inbounds @simd for i in iterator_index(fielda)
		fielda.e_SXY[i] += fieldb.e_SXY[i]
	end
end

function Base.:copy(field::FieldSpaceVectorial{T,D,X,Y}) where {T,D,X,Y}
	return FieldSpaceVectorial{T,D,X,Y}(deepcopy(field.x_X), deepcopy(field.y_Y), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

function copy_differentD(field::FieldSpaceVectorial{T,D,X,Y}) where {T,D,X,Y}
	return FieldSpaceVectorial{T,-D,X,Y}(deepcopy(field.x_X), deepcopy(field.y_Y), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

CartesianIndices(field::FieldSpaceVectorial) = Base.CartesianIndices((3, length(field.x_X), length(field.y_Y)))
LinearIndices(field::FieldSpaceVectorial) = Base.LinearIndices((3, length(field.x_X), length(field.y_Y)))

mutable struct FieldSpaceScalarRadialSymmetric{T,D,X<:AbstractVector{T},Y<:AbstractVector{Complex{T}}} <: AbstractFieldSpace{T,D}
	r_R::X
	e_SXY::Y
	λ::T
	n::Complex{T}
	ref::ReferenceFrame{T}
	function FieldSpaceScalarRadialSymmetric{T,D,X,Y}(r_R, e_SXY, λ, n, ref) where {T,D,X,Y}
		length(r_R) == length(e_SXY) || error("The length of r_R must be the same as the length of e_SXY.")
		(abs(ref.x) > 1E-10 || abs(ref.y) > 1E-10 || abs(ref.θ) > 1E-5 || abs(ref.ϕ) > 1E-5) && error("The referenceframe of a symetric field must be with x, y, θ and ϕ equal to 0")
		return new{T,D,X,Y}(r_R, e_SXY, λ, n, ref);
	end
	function FieldSpaceScalarRadialSymmetric{T,D}(r_R::X, e_SXY::Y, λ, n, ref) where {T,D,X,Y}
		return  FieldSpaceScalarRadialSymmetric{T,D,X,Y}(r_R, e_SXY, λ, n, ref)
	end
end

function FieldSpaceScalarRadialSymmetric_uniform(::Type{T}, r_R, λ, n, dir, ref) where T
	e_SXY = ones(Complex{T}, length(r_R))
	space = FieldSpaceScalarRadialSymmetric{T,dir}(r_R, e_SXY, λ, n, ref);
	space.e_SXY ./= √(intensity(space))
	return space
end
FieldSpaceScalarRadialSymmetric_uniform(arg...) = FieldSpaceScalarRadialSymmetric_uniform(Float64, arg...)

function FieldSpaceScalarRadialSymmetric_gaussian(::Type{T}, r_R, ω, λ, n, dir, ref) where T
	e_SXY = Vector{Complex{T}}(undef, length(r_R))
	space = FieldSpaceScalarRadialSymmetric{T,dir}(r_R, e_SXY, λ, n, ref);
	norm = 	√(8 / π / T(ω)^2)
	@inbounds @simd for i in iterator_index(space)
		space.e_SXY[i] = norm * exp(- 4 / T(ω)^2 * r_R[i]^2);
	end
	return space
end
FieldSpaceScalarRadialSymmetric_gaussian(arg...) = FieldSpaceScalarRadialSymmetric_gaussian(Float64, arg...)

function changereferenceframe!(fieldspace::FieldSpaceScalarRadialSymmetric, refnew::ReferenceFrame)
	fieldspace.ref == refnew || error("FieldSpaceSpaceRadialSymmetric cannot change the referenceframe.")
end

function intensity(space::FieldSpaceScalarRadialSymmetric{T}) where T
	int = zero(T)
	@inbounds @simd for i in iterator_index(space)
		Δr = Δvector(space.r_R, i)
		int += abs2(space.e_SXY[i]) * Δr * space.r_R[i]
	end
	return real(space.n) * 2π  * int 
end

mutable struct FieldSpaceSymmetric{T<:Real} <: AbstractFieldSpace{T}
	x_X::Vector{Float64}
	y_Y::Vector{Float64}
	e_SXY::Array{Complex{Float64}, 3}
	λ::Float64
	n::Number
	dir::Int64
	ref::ReferenceFrame
	function FieldSpaceSymmetric{T}(x_X, e_SXY, λ, n, dir, ref) where T
		length(x_X) != size(e_SXY, 2) ? error("The length of x_X must be the same as the size of e_SXY in the second dimension.") : nothing;
		size(e_SXY, 3) != 1 ? error("The length of e_SXY must be 1 in the third dimension.") : nothing;
		size(e_SXY, 1) != 1 ? error("The size of e_SXY must be 1 (scallar field).") : nothing;
		abs(ref.x) > 1E-10 || abs(ref.y) > 1E-10 || abs(ref.θ) > 1E-5 || abs(ref.ϕ) > 1E-5 ? error("The referenceframe of a symetric field must be with x, y, θ and ϕ equal to 0") : nothing;
		return new{T}(x_X, [0.], e_SXY, λ, n, dir >= 0 ? 1 : -1, ref);
	end
end

function FieldSpaceSymmetric_gaussian(x_X::AbstractVector{<:Real}, ω::Real, λ::Real, n::Number=1. +0*im, dir::Int=1, ref=ReferenceFrame())::FieldSpaceSymmetric
	norm = 	√(8 / π / ω^2);
	x_X = reshape(x_X, 1, :, 1);
	e_SXY = norm .* exp.(.- 4 ./ ω.^2 .* (x_X.^2));
	x_X = dropdims(x_X, dims=(1, 3));
	return FieldSpaceSymmetric(x_X, e_SXY, λ, n, dir, ref);
end

function changereferenceframe!(fieldspace::FieldSpaceSymmetric, refnew::ReferenceFrame)
	fieldspace.ref != refnew ? error("FieldSpaceSymmetric cannot change the referenceframe.") : nothing
end

function intensity(fieldspace::FieldSpaceSymmetric)::Float64
	e_SXY = abs2.(fieldspace.e_SXY) .* fieldspace.x_X;
	return fieldspace.n * ∫(e_SXY, fieldspace.x_X);
end

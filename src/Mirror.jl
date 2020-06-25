	struct Mirror{T} <: AbstractPropagationComponent{T}
    R::JolabFunction1D{T,T}
    n₁::JolabFunction1D{T,Complex{T}}
    n₂::JolabFunction1D{T,Complex{T}}
    ref::ReferenceFrame{T}
    Mirror{T}(R, n₁, n₂, ref) where T = new{T}(R, n₁, n₂, ref)
end
Mirror(R, n₁, n₂, ref) = Mirror{Float64}(R,n₁,n₂,ref)

function coefficient_geral(mirror::Mirror{T}, field::AbstractFieldSpace) where {T<:Real}
	sizeXY = length(field.x_X) * length(field.y_Y)
    tλ = complex(√(1 - mirror.R(field.λ)));
    rλ = complex(√mirror.R(field.λ));
	r12 = Diagonal(ones(Complex{T}, sizeXY))
	t12 = Diagonal(ones(Complex{T}, sizeXY))
	r21 = Diagonal(ones(Complex{T}, sizeXY))
	t21 = Diagonal(ones(Complex{T}, sizeXY))
	rmul!(r12.diag, rλ)
	rmul!(t12.diag, tλ)
	rmul!(r21.diag, -rλ)
	rmul!(t21.diag, tλ)

	fieldl = FieldSpace{T}(field.x_X, field.y_Y, field.e_SXY, field.λ, mirror.n₁(field.λ), -1, mirror.ref)
	fieldr = FieldSpace{T}(field.x_X, field.y_Y, field.e_SXY, field.λ, mirror.n₂(field.λ), 1, mirror.ref)
	return ScaterringMatrix{T, ScaterringMatrix{T, Diagonal{Complex{T},Vector{Complex{T}}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

function coefficient_geral(mirror::Mirror{T}, field::AbstractFieldAngularSpectrum) where {T<:Real}
	sizeXY = length(field.nsx_X) * length(field.nsy_Y)
    tλ = complex(√(1 - mirror.R(field.λ)));
    rλ = complex(√mirror.R(field.λ));
	r12 = Diagonal(ones(Complex{T}, sizeXY))
	t12 = Diagonal(ones(Complex{T}, sizeXY))
	r21 = Diagonal(ones(Complex{T}, sizeXY))
	t21 = Diagonal(ones(Complex{T}, sizeXY))
	rmul!(r12.diag, -rλ)
	rmul!(t12.diag, tλ)
	rmul!(r21.diag, rλ)
	rmul!(t21.diag, tλ)

	fieldl = FieldAngularSpectrum{T}(field.nsx_X, field.nsy_Y, field.e_SXY, field.λ, mirror.n₁(field.λ), -1, mirror.ref)
	fieldr = FieldAngularSpectrum{T}(field.nsx_X, field.nsy_Y, field.e_SXY, field.λ, mirror.n₂(field.λ), 1, mirror.ref)
	return ScaterringMatrix{T, typeof(r12), typeof(fieldl), typeof(fieldr)}(r12, t12, r21, t21, fieldl, fieldr)
end

@inline coefficient_specific(mirror::Mirror, field::AbstractFieldMonochromatic) = coefficient_geral(mirror, field)

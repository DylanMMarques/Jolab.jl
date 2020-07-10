struct Mirror{T} <: AbstractPropagationComponent{T}
    R::JolabFunction1D{T,T}
    n₁::JolabFunction1D{T,Complex{T}}
    n₂::JolabFunction1D{T,Complex{T}}
    ref::ReferenceFrame{T}
    Mirror{T}(R, n₁, n₂, ref) where T = new{T}(R, n₁, n₂, ref)
end
Mirror(R, n₁, n₂, ref) = Mirror{Float64}(R,n₁,n₂,ref)

function coefficient_general(mirror::Mirror{T}, field::AbstractFieldSpace) where {T<:Real}
	isapprox(field.n, field.dir > 0 ? mirror.n₁(field.λ) : mirror.n₂(field.λ), atol = @tol) || error("Field medium and mirror are different")
	checkinplane(field.ref, mirror.ref) || error("cannot do this")

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
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), Diagonal{Complex{T},Vector{Complex{T}}}, Diagonal{Complex{T},Vector{Complex{T}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

errorToDo() = error("Update in next versions")

function coefficient_general(mirror::Mirror{T}, field::AbstractFieldAngularSpectrum) where {T<:Real}
	isapprox(field.n, field.dir > 0 ? mirror.n₁(field.λ) : mirror.n₂(field.λ), atol = @tol) || error("Field medium and mirror are different")
	checkorientation(field.ref, mirror.ref) || errorToDo()

	sizeXY = length(field.nsx_X) * length(field.nsy_Y)
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

	if !checkposition(field.ref, mirror.ref)
		propM = propagationmatrix(field, mirror.ref)
		if field.dir > 0
			rmul!(r12, propM)
			lmul!(propM, r12)
			rmul!(t12, propM)
			lmul!(propM, t21)
		else
			conj!(propM.diag)
			rmul!(r21, propM)
			lmul!(propM, r21)
			rmul!(t21, propM)
			lmul!(propM, t12)
		end
	end

	if field.dir > 0
		fieldl = FieldAngularSpectrum{T}(field.nsx_X, field.nsy_Y, field.e_SXY, field.λ, mirror.n₁(field.λ), -1, field.ref)
		fieldr = FieldAngularSpectrum{T}(field.nsx_X, field.nsy_Y, field.e_SXY, field.λ, mirror.n₂(field.λ), 1, mirror.ref)
	else
		fieldl = FieldAngularSpectrum{T}(field.nsx_X, field.nsy_Y, field.e_SXY, field.λ, mirror.n₁(field.λ), -1, mirror.ref)
		fieldr = FieldAngularSpectrum{T}(field.nsx_X, field.nsy_Y, field.e_SXY, field.λ, mirror.n₂(field.λ), 1, field.ref)
	end
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), Diagonal{Complex{T},Vector{Complex{T}}}, Diagonal{Complex{T},Vector{Complex{T}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

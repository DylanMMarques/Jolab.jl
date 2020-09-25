struct ScatteringMatrix{T, L, R, X1 <: Union{AbstractMatrix{Complex{T}}, UniformScaling{Complex{T}}}, X2 <: AbstractMatrix{Complex{T}}} <: AbstractCoefficient{T,L,R}
	r₁₂::X1
	t₁₂::X2
	r₂₁::X1
	t₂₁::X2
	fieldl::L
	fieldr::R
end

function checkapplicability(fieldl::L, fieldr::R, scatM::ScatteringMatrix{T, L, R, X}, fieldi::Union{L,R}) where {T,X,L,R}
	samedefinitions(fieldl, scatM.fieldl) || return false
	samedefinitions(fieldr, scatM.fieldr) || return false
	samedefinitions(fieldi, fieldi.dir > 0 ? scatM.fieldl : scatM.fieldr) || return false
	return true
end

function checkapplicability(scatMs::AbstractVector{<:ScatteringMatrix})
	for i in 1:length(scatMs)-1
		samedefinitions(scatMs[i].fieldr, scatMs[i+1].fieldl) || return false
	end
	return true
end

function coefficient_general(coef1::ScatteringMatrix{T,L1,R1,X11}, coef2::ScatteringMatrix{T,L2,R2,X12}) where {T<:Real, L1, L2, R1, R2, X11<:UniformScaling, X12<:UniformScaling}
	r13 = UniformScaling(zero(Complex{T}))
	t31 = coef2.t₂₁ * coef1.t₂₁

	r31 = UniformScaling(zero(Complex{T}))
	t13 = coef1.t₁₂ * coef2.t₁₂

	return ScatteringMatrix{T, L1, R2, typeof(r13), typeof(t13)}(r13, t13, r31, t31, coef1.fieldl, coef2.fieldr)
end

function coefficient_general(coef1::ScatteringMatrix{T,L1,R1,X11}, coef2::ScatteringMatrix{T,L2,R2,X12}) where {T<:Real, L1, R1, L2, R2, X11<:UniformScaling, X12<:AbstractMatrix}
	r13 = coef1.t₁₂ * coef2.r₁₂ * coef1.t₂₁
	t13 = coef1.t₁₂ * coef2.t₁₂

	r31 = coef2.r₂₁
	t31 = coef2.t₂₁ * coef1.t₂₁

	return ScatteringMatrix{T, L1, R2, typeof(r13), typeof(t13)}(r13, t13, r31, t31, coef1.fieldl, coef2.fieldr)
end

function coefficient_general(coef1::ScatteringMatrix{T,L1,R1,X11}, coef2::ScatteringMatrix{T,L2,R2,X12}) where {T<:Real, L1, R1, L2, R2, X11<:AbstractMatrix, X12<:UniformScaling}
	r13 = coef1.r₁₂
	t13 = coef1.t₁₂ * coef2.t₁₂

	r31 = coef2.t₂₁ * coef1.r₂₁ * coef2.t₁₂
	t31 = coef2.t₂₁ * coef1.t₂₁

	return ScatteringMatrix{T, L1, R2, typeof(r13), typeof(t13)}(r13, t13, r31, t31, coef1.fieldl, coef2.fieldr)
end

function coefficient_general(coefs::AbstractVector{<:ScatteringMatrix{T}}) where {T<:Real}
	checkapplicability(coefs) || error("cannot be merged")
	sizeA = length(coefs)
	coefaux = coefs[sizeA]
	for i in length(coefs)-1:-1:1
		coefaux = coefficient_general(coefs[i], coefaux)
	end
	return coefaux
end

function lightinteraction!(fieldl::L, fieldr::R, scatM::ScatteringMatrix{T,L,R,X}, fieldi::Union{L,R}) where {T,X,L,R}
	if fieldi.dir > 0
		mul!(vec(fieldl.e_SXY), scatM.r₁₂, vec(fieldi.e_SXY))
		mul!(vec(fieldr.e_SXY), scatM.t₁₂, vec(fieldi.e_SXY))
	else
		mul!(vec(fieldr.e_SXY), scatM.r₂₁, vec(fieldi.e_SXY))
		mul!(vec(fieldl.e_SXY), scatM.t₂₁, vec(fieldi.e_SXY))
	end
end
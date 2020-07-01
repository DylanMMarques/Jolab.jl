struct ScatteringMatrix{T, L, R, X <: AbstractMatrix{Complex{T}}} <: AbstractCoefficient{T,L,R}
	r₁₂::X
	t₁₂::X
	r₂₁::X
	t₂₁::X
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

function coefficient_general(coefs::AbstractVector{<:ScatteringMatrix{T}}) where {T<:Real}
	checkapplicability(coefs) || error("cannot be merged")

	sizeA = length(coefs)
	r12 = coefs[sizeA].r₁₂
	t12 = coefs[sizeA].t₁₂
	r21 = coefs[sizeA].r₂₁
	t21 = coefs[sizeA].t₂₁

	for i in length(coefs)-1:-1:1
		if coefs[i].fieldr.ref != coefs[i+1].fieldl.ref
			propMatrix = propagationmatrix(coefs[i].fieldr, coefs[i+1].fieldl)
			# (r12, t12, r21, t21) are now (r23, t23, r23, t23) inside the equation
			rmul!(r12, propMatrix)
			lmul!(propMatrix, r12)

			rmul!(t12, propMatrix)
			lmul!(propMatrix, t21)
		end

		aux = inv(I - r12 * coefs[i].r₂₁)
		r13 = coefs[i].r₁₂ + coefs[i].t₂₁ * aux * r12 * coefs[i].t₁₂
		t31 = coefs[i].t₂₁ * aux * t21

		aux = inv(I - coefs[i].r₂₁ * r12)
		r31 = r21 + t12 * aux * coefs[i].r₂₁ * t21
		t13 = t12 * aux * coefs[i].t₁₂
		(r12, t12, r21, t21) = (r13, t13, r31, t31)
	end
	return ScatteringMatrix{T, typeof(coefs[1].fieldl), typeof(coefs[end].fieldr), typeof(r12)}(r12, t12, r21, t21, coefs[1].fieldl, coefs[end].fieldr)
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

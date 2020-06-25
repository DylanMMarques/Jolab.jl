struct ScaterringMatrix{T, X <: AbstractMatrix{Complex{T}}, L <: AbstractFieldMonochromatic{T}, R <: AbstractFieldMonochromatic{T}} <: AbstractCoefficient{T}
	r₁₂::X
	t₁₂::X
	r₂₁::X
	t₂₁::X
	fieldl::L
	fieldr::R
end

function samedefinitions(fieldl::L, fieldr::L) where L <: FieldAngularSpectrum
	isapprox(fieldl.nsx_X, fieldr.nsx_X, atol = @tol) || return false
	isapprox(fieldl.nsy_Y, fieldr.nsx_X, atol = @tol) || return false
	isapprox(fieldl.n, fieldr.n, atol = @tol) || return false
	isapprox(fieldl.λ, fieldr.λ, atol = @tol) || return false
	checkorientation(fieldl.ref, fieldr.ref) || return false
	return true
end

function samedefinitions(fieldl::L, fieldr::L) where L <: FieldSpace
	isapprox(fieldl.x_X, fieldr.x_X, atol = @tol) || return false
	isapprox(fieldl.y_Y, fieldr.x_X, atol = @tol) || return false
	isapprox(fieldl.n, fieldr.n, atol = @tol) || return false
	isapprox(fieldl.λ, fieldr.λ, atol = @tol) || return false
	checkorientation(fieldl.ref, fieldr.ref) || return false
	return true
end

function samedefinitions(fieldl::L, fieldr::L) where L <: FieldModes
	(checkorientation(fieldl.ref, fieldr.ref) && checkinline(fieldl.ref, fieldr.ref)) || return false
	fieldl.modes == fieldr.modes || return false
	return true
end

function checkapplicability(fieldl::L, fieldr::R, scatM::ScaterringMatrix{T,X,L,R}, fieldi::Union{L,R}) where {T,X,L,R}
	samedefinitions(fieldl, scatM.fieldl) || return false
	samedefinitions(fieldr, scatM.fieldr) || return false
	samedefinitions(fieldi, fieldi.dir > 0 ? scatM.fieldl : scatM.fieldr) || return false
	return true
end

function checkapplicability(scatMs::AbstractVector{<:ScaterringMatrix})
	for i in 1:length(scatMs)-1
		samedefinitions(scatMs[i].fieldr, scatMs[i+1].fieldl) || return false
	end
	return true
end

function mergeorientated_propagationcoefficient(coefs::AbstractVector{<:ScaterringMatrix{T}}) where {T<:Real}
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
	return ScaterringMatrix{T, typeof(r12),typeof(coefs[1].fieldl), typeof(coefs[end].fieldr)}(r12, t12, r21, t21, coefs[1].fieldl, coefs[end].fieldr)
end

function lightinteraction!(fieldl::L, fieldr::R, scatM::ScaterringMatrix{T,X,L,R}, fieldi::Union{L,R}) where {T,X,L,R}
	if fieldi.dir > 0
		mul!(vec(fieldl.e_SXY), scatM.r₁₂, vec(fieldi.e_SXY))
		mul!(vec(fieldr.e_SXY), scatM.t₁₂, vec(fieldi.e_SXY))
	else
		mul!(vec(fieldr.e_SXY), scatM.r₂₁, vec(fieldi.e_SXY))
		mul!(vec(fieldl.e_SXY), scatM.t₂₁, vec(fieldi.e_SXY))
	end
end

function lightinteraction(scatM::ScaterringMatrix{T,X,L,R}, fieldi::Union{L,R}) where {T,X,L,R}
	fieldl = deepcopy(scatM.fieldl)
	fieldr = deepcopy(scatM.fieldr)

	samedefinitions(fieldi, fieldi.dir > 0 ? scatM.fieldl : scatM.fieldr)
	changereferential!(fieldi, fieldi.dir > 0 ? scatM.fieldl.ref : scatM.fieldr.ref)
	lightinteraction!(fieldl, fieldr, scatM, fieldi)
	return (fieldl, fieldr)
end

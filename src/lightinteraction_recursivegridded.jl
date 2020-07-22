function intensity_p(field::Union{FieldAngularSpectrum{T}, FieldSpace{T}}) where T
	int = zero(T)
	@inbounds @simd for i in field.e_SXY
		int += abs2(i)
	end
	return int * real(field.n)
end

function lightinteraction_recursivegridded!(fieldl::AbstractFieldMonochromatic{T}, fieldr::AbstractFieldMonochromatic{T}, coefs::AbstractVector{<:AbstractCoefficient{T}}, fieldi::AbstractFieldMonochromatic{T}; rtol = 1E-3one(T)::Real) where {T<:Real}
	# miss check if this can be done
	sizeL = length(coefs) + 1;

	fields_r = Vector{AbstractFieldMonochromatic{T}}(undef, sizeL)
	fields_l = Vector{AbstractFieldMonochromatic{T}}(undef, sizeL)
	int_l = zeros(T, sizeL)
	int_r = zeros(T, sizeL)
	cval_l = ones(T, sizeL)
	cval_r = ones(T, sizeL)

	for i in 1:sizeL-1
		(fields_l[i], tmp) = getfields_lr(coefs[i])
		fields_l[i].dir = -1
		vec(fields_l[i].e_SXY) .= zero(Complex{T})
		fields_r[i] = copy(fields_l[i])
		fields_r[i].dir = 1
	end

	(tmp, fields_l[sizeL]) = getfields_lr(coefs[sizeL-1])
	fields_l[sizeL].dir = -1
	fields_l[sizeL].e_SXY .= zero(Complex{T})
	fields_r[sizeL] = copy(fields_l[sizeL])
	fields_r[sizeL].dir = 1

	fields_aux_r = deepcopy(fields_r)
	fields_aux_l = deepcopy(fields_l)
	rtol = intensity_p(fieldi) * rtol^2

	fieldi.dir > 0 ? vec(fields_r[1].e_SXY) .= vec(fieldi.e_SXY) : vec(fields_l[sizeL].e_SXY) .= vec(fieldi.e_SXY)
	initial_int = intensity_p(fieldi)
	fieldi.dir > 0 ? int_r[1] = initial_int : int_l[sizeL] = initial_int
	fields2_l = deepcopy(fields_l)
	fields2_r = deepcopy(fields_r)

	i = 1
	toSave_l = fields_l
	toSave_r = fields_r
	@inbounds while true
		if i % 2 == 1
			toSave_l = fields_l
			toSave_r = fields_r
			iE_l = fields2_l
			iE_r = fields2_r
		else
			toSave_l = fields2_l
			toSave_r = fields2_r
			iE_l = fields_l
			iE_r = fields_r
		end

		for mls in 1:sizeL
			if mls < sizeL
				if int_r[mls] > rtol * 1E-2
					iE_r[mls].ref == coefs[mls].fieldl.ref || tobedone()
					lightinteraction!(fields_aux_l[mls], fields_aux_r[mls+1], coefs[mls], iE_r[mls])
					add_inplace!(toSave_l[mls], fields_aux_l[mls])
					add_inplace!(toSave_r[mls+1], fields_aux_r[mls+1])
				else
					# add_inplace!(toSave_r[mls], iE_r[mls])
				end
				vec(iE_r[mls].e_SXY) .= zero(Complex{T})
			end
			if mls > 1
				if int_l[mls] > rtol * 1E-2
					iE_l[mls].ref == coefs[mls-1].fieldr.ref || tobedone()
					lightinteraction!(fields_aux_l[mls-1], fields_aux_r[mls], coefs[mls-1], iE_l[mls])
					add_inplace!(toSave_l[mls-1], fields_aux_l[mls-1])
					add_inplace!(toSave_r[mls], fields_aux_r[mls])
				else
					# add_inplace!(toSave_l[mls], iE_l[mls])
				end
				vec(iE_l[mls].e_SXY) .= zero(Complex{T})
			end
		end
		for mls in 1:sizeL
			int_l[mls] = intensity_p(toSave_l[mls])
			int_r[mls] = intensity_p(toSave_r[mls])
		end

		sum(view(int_l,2:sizeL)) + sum(view(int_r, 1:sizeL-1)) < rtol && break
		sum(int_l) + sum(int_r) > 10initial_int && (println("lightinteraction_recursivegridded is not converging. Change cval (increase is recommended)."); break)
		i > 10000 && error("Max number of iterations achieved")
		i += 1
	end
	println("Interactions until convergence: ", i)
	vec(fieldl.e_SXY) .= vec(toSave_l[1].e_SXY)
	vec(fieldr.e_SXY) .= vec(toSave_r[sizeL].e_SXY)
end

function lightinteraction_recursivegridded(coefs::AbstractVector{<:AbstractCoefficient{T}}, field_inc::AbstractFieldMonochromatic{T}; rtol = 1E-3) where T
	(fieldl, tmp) = getfields_lr(coefs[1])
	(tmp, fieldr) = getfields_lr(coefs[end])

	samedefinitions(field_inc, field_inc.dir > 0 ? fieldl : fieldr) || error("Cannot do this")
	lightinteraction_recursivegridded!(fieldl, fieldr, coefs, field_inc, rtol = rtol)
	return (fieldl, fieldr)
end

function lightinteraction_recursivegridded(comps::AbstractVector{<:AbstractOpticalComponent{T}}, fieldi::AbstractFieldMonochromatic{T}; rtol = 1E-3) where T
	coefs = coefficient_specific(comps, fieldi)
	return lightinteraction_recursivegridded(coefs, fieldi, rtol = rtol)
end

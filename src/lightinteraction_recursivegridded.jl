function intensity_p(field::Union{FieldAngularSpectrum{T}, FieldSpace{T}}) where T
	int = zero(T)
	@inbounds @simd for i in field.e_SXY
		int += abs2(i)
	end
	return int * real(field.n)
end

function lightinteraction_recursivegridded!(fieldl::AbstractFieldMonochromatic{T}, fieldr::AbstractFieldMonochromatic{T}, coefs::AbstractVector{<:AbstractCoefficient{T}}, fieldi::AbstractFieldMonochromatic{T}; rtol = 1E-3one(T)::Real, printBool = true) where {T<:Real}
	# miss check if this can be done
	sizeL = length(coefs) + 1;

	fields_r = Vector{AbstractFieldMonochromatic{T,1}}(undef, sizeL)
	fields_l = Vector{AbstractFieldMonochromatic{T,-1}}(undef, sizeL)
	int_l = zeros(T, sizeL)
	int_r = zeros(T, sizeL)
	cval_l = ones(T, sizeL)
	cval_r = ones(T, sizeL)

	for i in 1:sizeL-1
		(fields_l[i], tmp) = getfields_lr(coefs[i])
		vec(fields_l[i].e_SXY) .= zero(Complex{T})
		fields_r[i] = copy_differentD(fields_l[i])
	end

	(tmp, fields_r[sizeL]) = getfields_lr(coefs[sizeL-1])
	fields_r[sizeL].e_SXY .= zero(Complex{T})
	fields_l[sizeL] = copy_differentD(fields_r[sizeL])

	fields_aux_r = deepcopy(fields_r)
	fields_aux_l = deepcopy(fields_l)
	rtol = intensity_p(fieldi) * rtol^2

	dir(fieldi) > 0 ? vec(fields_r[1].e_SXY) .= vec(fieldi.e_SXY) : vec(fields_l[sizeL].e_SXY) .= vec(fieldi.e_SXY)
	initial_int = intensity_p(fieldi)
	dir(fieldi) > 0 ? int_r[1] = initial_int : int_l[sizeL] = initial_int
	fields2_l = deepcopy(fields_l)
	fields2_r = deepcopy(fields_r)
	# fields2_r[end] = fields_r[end]
	# fields2_l[1] = fields_l[1]

	i = 1
	toSave_l = fields_l
	toSave_r = fields_r
	min_int = initial_int
	converge = false
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

		for mls in 1:sizeL-1
			if int_r[mls] > initial_int * @tol
				iE_r[mls].ref == coefs[mls].fieldl.ref || tobedone()
				lightinteraction!(fields_aux_l[mls], fields_aux_r[mls+1], coefs[mls], iE_r[mls])
				add_inplace!(toSave_l[mls], fields_aux_l[mls])
				add_inplace!(toSave_r[mls+1], fields_aux_r[mls+1])
			else
				# add_inplace!(toSave_r[mls], iE_r[mls])
			end
			vec(iE_r[mls].e_SXY) .= zero(Complex{T})
		end
		for mls in 2:sizeL
			if int_l[mls] > initial_int * @tol
				iE_l[mls].ref == coefs[mls-1].fieldr.ref || tobedone()
				lightinteraction!(fields_aux_l[mls-1], fields_aux_r[mls], coefs[mls-1], iE_l[mls])
				add_inplace!(toSave_l[mls-1], fields_aux_l[mls-1])
				add_inplace!(toSave_r[mls], fields_aux_r[mls])
			else
				# add_inplace!(toSave_l[mls], iE_l[mls])
			end
			vec(iE_l[mls].e_SXY) .= zero(Complex{T})
		end
		Threads.@threads for mls in 1:sizeL
			int_l[mls] = intensity_p(toSave_l[mls])
			int_r[mls] = intensity_p(toSave_r[mls])
		end
		now_int = sum(view(int_l,2:sizeL)) + sum(view(int_r, 1:sizeL-1))
		now_int < rtol && (println(""); println("Interactions until convergence: ", i); converge = true; break)

		(now_int < min_int) && (min_int = now_int)

		if sum(int_l) + sum(int_r) > 10 * initial_int || now_int > 10 * min_int
			println("")
			println("Light intensity propagating forward:", int_r)
			println("Light intensity propagating backward:", int_l)
			println("lightinteraction_recursivegridded is not converging. Current number of iterations:", i)
			converge = false
			break
		end

		i > 10000 && (println("Max number of iterations achieved. Current light intensity:", (sum(int_l) + sum(int_r)) / initial_int); converge = false; break)
		i += 1
		if (i % 100 == 99) && printBool
			println("")
			println("Light intensity propagating forward:", int_r)
			println("Light intensity propagating backward:", int_l)
			println("convergence condition: ", sum(view(int_l,2:sizeL)) + sum(view(int_r, 1:sizeL-1)), " < ", rtol)
		end
	end
	vec(fieldl.e_SXY) .= vec(toSave_l[1].e_SXY)
	vec(fieldr.e_SXY) .= vec(toSave_r[sizeL].e_SXY)
	return converge
end

function lightinteraction_recursivegridded(coefs::AbstractVector{<:AbstractCoefficient{T}}, field_inc::AbstractFieldMonochromatic{T}; rtol = 1E-3, printBool = true) where T
	(fieldl, tmp) = getfields_lr(coefs[1])
	(tmp, fieldr) = getfields_lr(coefs[end])

	samedefinitions(field_inc, dir(field_inc) > 0 ? fieldl : fieldr) || error("Cannot do this")
	converge = lightinteraction_recursivegridded!(fieldl, fieldr, coefs, field_inc, rtol = rtol, printBool = printBool)
	return (fieldl, fieldr, converge)
end

function lightinteraction_recursivegridded(comps::AbstractVector{<:AbstractOpticalComponent{T}}, fieldi::AbstractFieldMonochromatic{T}; rtol = 1E-3, printBool = true) where T
	coefs = coefficient_specific(comps, fieldi)
	return lightinteraction_recursivegridded(coefs, fieldi, rtol = rtol, printBool = printBool)
end

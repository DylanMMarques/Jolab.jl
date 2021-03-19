function lightinteraction_recursivegridded!(fieldl::AbstractFieldMonochromatic{T}, fieldr::AbstractFieldMonochromatic{T}, coefs::AbstractVector{<:AbstractCoefficient{T}}, fieldi::AbstractFieldMonochromatic{T}; rtol = 1E-3one(T)::Real, printBool = true) where {T<:Real}
	# miss check if this can be done
	sizeL = length(coefs) + 1;

	fields_r = Vector{AbstractFieldMonochromatic{T,1}}(undef, sizeL)
	fields_l = Vector{AbstractFieldMonochromatic{T,-1}}(undef, sizeL)
	int_l, int_r = zeros(T, sizeL), zeros(T, sizeL)

	for i in 1:sizeL-1
		(fields_l[i], tmp) = getfields_lr(coefs[i])
		fields_l[i].e_SXY .= zero(Complex{T})
		fields_r[i] = copy_differentD(fields_l[i])
	end

	(tmp, fields_r[sizeL]) = getfields_lr(coefs[sizeL-1])
	fields_r[sizeL].e_SXY .= zero(Complex{T})
	fields_l[sizeL] = copy_differentD(fields_r[sizeL])
	fields_aux_r = deepcopy(fields_r)
	fields_aux_l = deepcopy(fields_l)
	fields2_l = deepcopy(fields_l)
	fields2_r = deepcopy(fields_r)
	fields_r[sizeL] = fields2_r[sizeL]
	fields_l[1] = fields2_l[1]

	rtol = intensity(fieldi) * rtol^2

	dir(fieldi) > 0 ? fields_r[1].e_SXY .= fieldi.e_SXY : fields_l[sizeL].e_SXY .= fieldi.e_SXY
	initial_int = intensity(fieldi)
	dir(fieldi) > 0 ? int_r[1] = initial_int : int_l[sizeL] = initial_int

	i = 1
	toSave_l, toSave_r = fields_l, fields_r
	min_int = initial_int
	converge = false

	# Create status bar and relation between time and intensity
	m_prog = 100 / (log(rtol) - log(initial_int))
	b_prog = - m_prog * log(initial_int)
	p = Progress(100)

	@inbounds while true
		if i % 2 == 1
			(toSave_l, toSave_r) = (fields_l, fields_r)
			(iE_l, iE_r) = (fields2_l, fields2_r)
		else
			(toSave_l, toSave_r) = (fields2_l, fields2_r)
			(iE_l, iE_r) = (fields_l, fields_r)
		end

		Threads.@threads for mls in 1:sizeL-1
			if int_r[mls] > rtol * @tol
				iE_r[mls].ref == coefs[mls].fieldl.ref || tobedone()
				lightinteraction!(fields_aux_l[mls], fields_aux_r[mls+1], coefs[mls], iE_r[mls])
				add!(toSave_l[mls], fields_aux_l[mls])
				add!(toSave_r[mls+1], fields_aux_r[mls+1])
			end
			setzeros!(iE_r[mls])
		end
		Threads.@threads for mls in 2:sizeL
			if int_l[mls] > rtol * @tol
				iE_l[mls].ref == coefs[mls-1].fieldr.ref || tobedone()
				lightinteraction!(fields_aux_l[mls-1], fields_aux_r[mls], coefs[mls-1], iE_l[mls])
				add!(toSave_l[mls-1], fields_aux_l[mls-1])
				add!(toSave_r[mls], fields_aux_r[mls])
			end
			setzeros!(iE_l[mls])
		end
		Threads.@threads for mls in 1:sizeL
			int_l[mls] = intensity(toSave_l[mls])
			int_r[mls] = intensity(toSave_r[mls])
		end
		now_int = sum(view(int_l,2:sizeL)) + sum(view(int_r, 1:sizeL-1))

		# Update status bar
		status = round(Int, log(now_int) * m_prog + b_prog)
		status < 0 && (status = 0)
		update!(p, status)

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
		# if (i % 100 == 99) && printBool
		# 	println("")
		# 	println("Light intensity propagating forward:", int_r)
		# 	println("Light intensity propagating backward:", int_l)
		# 	println("convergence condition: ", sum(view(int_l,2:sizeL)) + sum(view(int_r, 1:sizeL-1)), " < ", rtol)
		# end
	end
	fieldl.e_SXY .= fields_l[1].e_SXY
	fieldr.e_SXY .= fields_r[sizeL].e_SXY
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

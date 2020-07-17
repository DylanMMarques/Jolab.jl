function intensity_p(field::Union{FieldAngularSpectrum{T}, FieldSpace{T}}) where T
	int = zero(T)
	@inbounds @simd for i in field.e_SXY
		int += abs2(i)
	end
	return int * real(field.n)
end

function lightinteraction_recursivegridded!(fieldl::AbstractFieldMonochromatic{T}, fieldr::AbstractFieldMonochromatic{T}, coefs::AbstractVector{<:AbstractCoefficient{T}}, fieldi::AbstractFieldMonochromatic{T}; rtol = 1E-3one(T)::Real, cval = one(T)::Real) where {T<:Real}
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

	i = 1
	arg_l = 1
	arg_r = 1
	while true
		# Select the next field to consider
		max_r = zero(T)
		max_l = zero(T)
		for a in 1:sizeL-1
			if max_r < cval_r[a] * int_r[a]
				max_r = cval_r[a] * int_r[a]
				arg_r = a
			end
			if max_l < cval_l[a+1] * int_l[a+1]
				max_l = cval_l[a+1] * int_l[a+1]
				arg_l = a
			end
		end

		if max_r > max_l # Field is propagating forward
			max_r < rtol && break
			mls = arg_r
			fields_r[mls].ref == coefs[mls].fieldl.ref || tobedone()
			lightinteraction!(fields_aux_l[mls], fields_aux_r[mls+1], coefs[mls], fields_r[mls])
			add_inplace!(fields_l[mls], fields_aux_l[mls])
			add_inplace!(fields_r[mls+1], fields_aux_r[mls+1])
			int_l[mls] = intensity_p(fields_l[mls])
			int_r[mls+1] = intensity_p(fields_r[mls+1])
			vec(fields_r[mls].e_SXY) .= zero(Complex{T})
			int_r[mls] = zero(T)
			cval_r[mls] = one(T)
		else # field is going backward
			max_l < rtol && break
			mls = arg_l + 1 # need to add +1 because length start from 2
			fields_l[mls].ref == coefs[mls-1].fieldr.ref || tobedone()
			lightinteraction!(fields_aux_l[mls-1], fields_aux_r[mls], coefs[mls-1], fields_l[mls])
			add_inplace!(fields_l[mls-1], fields_aux_l[mls-1])
			add_inplace!(fields_r[mls], fields_aux_r[mls])
			int_l[mls-1] = intensity_p(fields_l[mls-1])
			int_r[mls] = intensity_p(fields_r[mls])
			vec(fields_l[mls].e_SXY) .= zero(Complex{T})
			int_l[mls] = zero(T)
			cval_l[mls] = one(T)
		end
		sum(int_l) + sum(int_r) > 10initial_int && (println("lightinteraction_recursivegridded is not converging. Change cval (increase is recommended)."); break)
		i > 100000 && error("Max number of iterations achieved")
		i += 1

		cval_l .*= cval
		cval_r .*= cval

		# @show cval_l
		# @show int_l
		# @show cval_r
		# @show int_r
		# @show initial_int

	end
	# @show i
	vec(fieldl.e_SXY) .= vec(fields_l[1].e_SXY)
	vec(fieldr.e_SXY) .= vec(fields_r[sizeL].e_SXY)
end

function lightinteraction_recursivegridded(coefs::AbstractVector{<:AbstractCoefficient{T}}, field_inc::AbstractFieldMonochromatic{T}; rtol = 1E-3, cval = one(T)) where T
	(fieldl, tmp) = getfields_lr(coefs[1])
	(tmp, fieldr) = getfields_lr(coefs[end])

	samedefinitions(field_inc, field_inc.dir > 0 ? fieldl : fieldr) || error("Cannot do this")
	lightinteraction_recursivegridded!(fieldl, fieldr, coefs, field_inc, rtol = rtol, cval = cval)
	return (fieldl, fieldr)
end

function lightinteraction_recursivegridded(comps::AbstractVector{<:AbstractOpticalComponent{T}}, fieldi::AbstractFieldMonochromatic{T}; rtol = 1E-3, cval = one(T)) where T
	coefs = coefficient_specific(comps, fieldi)
	return lightinteraction_recursivegridded(coefs, fieldi, rtol = rtol, cval = cval)
end

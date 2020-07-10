function intensity_p(field::Union{FieldAngularSpectrum{T}, FieldSpace{T}}) where T
	int = zero(T)
	@inbounds @simd for i in field.e_SXY
		int += abs2(i)
	end
	return int * real(field.n)
end

function lightinteraction_recursivegridded!(fieldl::AbstractFieldMonochromatic{T}, fieldr::AbstractFieldMonochromatic{T}, coefs::AbstractVector{<:AbstractCoefficient{T}}, fieldi::AbstractFieldMonochromatic{T}; rtol = convert(T,1E-4)::Real, cval = convert(T,2)::Real) where {T<:Real}
	# miss check if this can be done
	sizeL = length(coefs) + 1;

	fields_r = Vector{AbstractFieldMonochromatic{T}}(undef, sizeL)
	fields_l = Vector{AbstractFieldMonochromatic{T}}(undef, sizeL)
	int_l = zeros(T, sizeL)
	int_r = zeros(T, sizeL)

	for i in 1:sizeL-1
		(fields_l[i], tmp) = getfields_lr(coefs[i])
		fields_l[i].dir = -1
		vec(fields_l[i].e_SXY) .= zero(Complex{T})
		fields_r[i] = deepcopy(fields_l[i])
		fields_r[i].dir = 1
	end

	(tmp, fields_l[sizeL]) = getfields_lr(coefs[sizeL-1])
	fields_l[sizeL].dir = -1
	fields_l[sizeL].e_SXY .= zero(Complex{T})
	fields_r[sizeL] = deepcopy(fields_l[sizeL])
	fields_r[sizeL].dir = 1

	fields_aux_r = deepcopy(fields_r)
	fields_aux_l = deepcopy(fields_l)
	rtol = intensity_p(fieldi) * rtol^2

	fieldi.dir > 0 ? vec(fields_r[1].e_SXY) .= vec(fieldi.e_SXY) : vec(fields_l[sizeL].e_SXY) .= vec(fieldi.e_SXY)
	initial_int = intensity_p(fieldi)
	fieldi.dir > 0 ? int_r[1] = initial_int : int_l[sizeL] = initial_int

	i = 1
	while true
		# Select the next field to consider
		(max_l, arg_l) = findmax(view(int_l,2:sizeL)) # Removes the first as it is the field going out
		(max_r, arg_r) = findmax(view(int_r,1:sizeL-1)) # Removes the last as it is the field going out

		view(int_l, 2:sizeL) .*= cval
		view(int_r, 1:sizeL-1) .*= cval
		if max_r > max_l # Field is propagating forward
			max_r < rtol && break
			mls = arg_r
			fields_r[mls].ref == coefs[mls].fieldl.ref || tobedone()
			lightinteraction!(fields_aux_l[mls], fields_aux_r[mls+1], coefs[mls], fields_r[mls])
			vec(fields_l[mls].e_SXY) .+= vec(fields_aux_l[mls].e_SXY)
			vec(fields_r[mls+1].e_SXY) .+= vec(fields_aux_r[mls+1].e_SXY)
			int_l[mls] = intensity_p(fields_l[mls])
			int_r[mls+1] = intensity_p(fields_r[mls+1])
			int_r[mls] = zero(T)
			vec(fields_r[mls].e_SXY) .= zero(Complex{T})
		else # field is going backward
			max_l < rtol && break
			mls = arg_l + 1 # need to add +1 because length start from 2
			fields_l[mls].ref == coefs[mls-1].fieldr.ref || tobedone()
			lightinteraction!(fields_aux_l[mls-1], fields_aux_r[mls], coefs[mls-1], fields_l[mls])
			vec(fields_l[mls-1].e_SXY) .+= vec(fields_aux_l[mls-1].e_SXY)
			vec(fields_r[mls].e_SXY) .+= vec(fields_aux_r[mls].e_SXY)
			int_l[mls-1] = intensity_p(fields_l[mls-1])
			int_r[mls] = intensity_p(fields_r[mls])
			int_l[mls] = zero(T)
			vec(fields_l[mls].e_SXY) .= zero(Complex{T})
		end
		sum(int_l) + sum(int_r) > 10initial_int && error("lightinteraction_recursivegridded is not converging. Change cval (increase is recommended).")
		i += 1
	end
	vec(fieldl.e_SXY) .= vec(fields_l[1].e_SXY)
	vec(fieldr.e_SXY) .= vec(fields_r[sizeL].e_SXY)
end

function lightinteraction_recursivegridded(coefs::AbstractVector{<:AbstractCoefficient{T}}, field_inc::AbstractFieldMonochromatic{T}; rtol = 1E-4, cval = 2one(T)) where T
	(fieldl, tmp) = getfields_lr(coefs[1])
	(tmp, fieldr) = getfields_lr(coefs[end])

	samedefinitions(field_inc, field_inc.dir > 0 ? fieldl : fieldr) || error("Cannot do this")
	lightinteraction_recursivegridded!(fieldl, fieldr, coefs, field_inc, rtol = rtol, cval = cval)
	return (fieldl, fieldr)
end

function lightinteraction_recursivegridded(comps::AbstractVector{<:AbstractOpticalComponent{T}}, fieldi::AbstractFieldMonochromatic{T}; rtol = 1E-4, cval = 2one(T)) where T
	coefs = coefficient_specific(comps, fieldi)
	return lightinteraction_recursivegridded(coefs, fieldi, rtol = rtol, cval = cval)
end

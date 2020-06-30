function lightinteraction_recursivegridded!(fieldl::AbstractFieldMonochromatic{T}, fieldr::AbstractFieldMonochromatic{T}, coefs::AbstractVector{<:AbstractCoefficient{T}}, fieldi::AbstractFieldMonochromatic{T} , thresold::Real) where {T<:Real}
	# miss check if this can be donne
	sizeL = length(coefs) + 1;

	fields_r = Vector{AbstractFieldMonochromatic{T}}(undef, sizeL)
	fields_l = Vector{AbstractFieldMonochromatic{T}}(undef, sizeL)
	int_l = zeros(T, sizeL)
	int_r = zeros(T, sizeL)

	for i in 1:sizeL-1
		(fields_l[i], tmp) = getfields_lr(coefs[i])
		fields_l.dir = -1
		fields_l.e_SXY .= zero(Complex{T})
		fields_r[i] = deepcopy(fields_l[i])
		fields_r.dir = 1
		fields_r.e_SXY .= zero(Complex{T})
	end
	(tmp, fields_l[sizeL]) = getfields_lr(coefs[sizeL-1])
	fields_l.dir = -1
	fields_l.e_SXY .= zero(Complex{T})
	fields_r[sizeL] = deepcopy(fields_l[sizeL])
	fields_r.dir = 1
	fields_r.e_SXY .= zero(Complex{T})

	rtol = intensity(fieldi) * convert(T, thresold^2)

	fieldi.dir > 0 ? fields_r[1] = deepcopy(fieldi) : fields_l[sizeL] = deepcopy(fieldi)

	i = 1
	mls = 1
	isForward = true

	while true
		# Select the next field to consider
		(maxBack, argBack) = findmax(view(intBackward_L, 2:sizeL)) # Removes the first as it is the field going out
		(maxForw, argForw) = findmax(view(intForward_L, 1:sizeL-1)) # Removes the last as it is the field going out
		# @show isForward
		# @show intBackward_L
		# @show intForward_L
		# @show (sum(intBackward_L) + sum(intForward_L)) / rtol * thresold^2
		# readline()
		# (maxForw + maxBack) / rtol * thresold^2 > 10 && error("Energy is being created.")
		isForward = maxForw > maxBack
		if isForward
			maxForw < rtol && break
			mls = argForw
			field_inc = fields_r[mls]
		else
			maxBack < rtol && break
			mls = argBack + 1 # need to add +1 because length start from 2
			field_inc = fields_l[mls]
		end

		if isForward
			# need to check this
			lightinteraction!(fields_l[mls], fields_r[mls], coefs[mls], field_inc)
			intBackward_L[mls] = intensity(fields_l[mls])
			intForward_L[mls+1] = intensity(fields_r[mls])
			intForward_L[mls] = zero(T)
		else
			#need to check this
			lightinteraction!(fields_l[mls], fields_r[mls], coefs[mls], field_inc)
			intBackward_L[mls] = intensity(fields_l[mls])
			intForward_L[mls+1] = intensity(fields_r[mls])
			intForward_L[mls] = zero(T)
		end
		e_SXY .= zero(Complex{T})
		i = i + 1;
		i > 5000 && break
	end
	ebackward_SXY .= view(eBackwardLayers_SXYL,:,:,:,1)
	eforward_SXY .= view(eForwardLayers_SXYL,:,:,:,sizeL)
	return (ebackward_SXY, eforward_SXY)
end

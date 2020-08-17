struct PhantomLayer{T} <: AbstractOpticalComponent{T}
	ref::ReferenceFrame{T}
end
PhantomLayer(ref) = PhantomLayer{Float64}(ref)

function coefficient_general(phantom::PhantomLayer{T}, field::AbstractFieldSpace) where {T<:Real}
	checkinplane(field.ref, phatom.ref) || error("cannot do this")

	sizeXY = length(field.x_X) * length(field.y_Y)
	r = spzeros(Complex{T}, sizeXY, sizeXY)
	t = Diagonal(ones(Complex{T}, sizeXY))

	fieldl = FieldSpace{T}(field.x_X, field.y_Y, field.e_SXY, field.λ, field.n, -1, phantom.ref)
	fieldr = FieldSpace{T}(field.x_X, field.y_Y, field.e_SXY, field.λ, field.n, 1, phantom.ref)
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), typeof(r), Diagonal{Complex{T},Vector{Complex{T}}}}(r, t, r, t, fieldl, fieldr)
end

function coefficient_general(phantom::PhantomLayer{T}, field::AbstractFieldAngularSpectrum) where {T<:Real}
	checkorientation(field.ref, phantom.ref) || errorToDo()

	sizeXY = length(field.nsx_X) * length(field.nsy_Y)

	r = spzeros(Complex{T}, sizeXY, sizeXY)
	t = Diagonal(ones(Complex{T}, sizeXY))

	if !checkposition(field.ref, phantom.ref)
		propM = propagationmatrix(field, phantom.ref)
		if field.dir > 0
			rmul!(t, propM)
		else
			conj!(propM.diag)
			rmul!(t, propM)
		end
	end

	if field.dir > 0
		fieldl = FieldAngularSpectrum{T}(field.nsx_X, field.nsy_Y, field.e_SXY, field.λ, field.n, -1, field.ref)
		fieldr = FieldAngularSpectrum{T}(field.nsx_X, field.nsy_Y, field.e_SXY, field.λ, field.n, 1, phantom.ref)
	else
		fieldl = FieldAngularSpectrum{T}(field.nsx_X, field.nsy_Y, field.e_SXY, field.λ, field.n, -1, phantom.ref)
		fieldr = FieldAngularSpectrum{T}(field.nsx_X, field.nsy_Y, field.e_SXY, field.λ, field.n, 1, field.ref)
	end
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), typeof(r), Diagonal{Complex{T},Vector{Complex{T}}}}(r, t, r, t, fieldl, fieldr)
end

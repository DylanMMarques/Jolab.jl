mutable struct SpatialLightModulator{T} <: AbstractOpticalComponent{T}
	r::JolabFunction2D{T}
	t::JolabFunction2D{T}
    ref::ReferenceFrame{T}
end

SpatialLightModulator_slm(t, ref) = return SpatialLightModulator(0, t, ref)
SpatialLightModulator_mask(t, ref) = return SpatialLightModulator(1 - t, t, ref)
SpatialLightModulator_reflectivemask(r, ref) = return SpatialLightModulator(r, 1 - r, ref)

function coefficient_general(slm::SpatialLightModulator, field::FieldSpaceScalar{T}) where {T<:Real}
	checkorientation(field.ref, slm.ref) || errorToDo()
    fieldi = changereferenceframe(field, slm.ref)

	sizeI = length(field.e_SXY)
	cart = CartesianIndices(field)
	scat = get_scatteringmatrixtype(slm, field)
	@inbounds @simd for i in 1:sizeI
		scat.t₁₂.diag[i] = slm.t(field.x_X[cart[i][2]], field.y_Y[cart[i][3]])
		scat.r₁₂.diag[i] = slm.r(field.x_X[cart[i][2]], field.y_Y[cart[i][3]])
	end
end

function get_scatteringmatrixtype(slm::SpatialLightModulator, field::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y}
	sizeI = length(field.e_SXY)
	r = Diagonal(Vector{Complex{T}}(undef, sizeI))
	t = Diagonal(Vector{Complex{T}}(undef, sizeI))
	(fieldl, fieldr) = getfields_lr(slm, field)
	return ScatteringMatrix{T, FieldSpaceScalar{T,-1,X,Y}, FieldSpaceScalar{T,-1,X,Y}, Diagonal{Complex{T},Vector{Complex{T}}}, Diagonal{Complex{T},Vector{Complex{T}}}}(r, t, r, t, fieldl, fieldr)
end

function getfields_lr(slm::SpatialLightModulator, field::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y}
	fieldl = FieldSpaceScalar{T,-1,X,Y}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, -1, fieldi.ref)
	fieldr = FieldSpaceScalar{T,1,X,Y}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, 1, fieldi.ref)
	return (fieldl, fieldr)
end

mutable struct SpatialLightModulator{T, A<:JolabFunction2D{T}, B<:JolabFunction2D{T}} <: AbstractOpticalComponent{T}
	r::A
	t::B
    ref::ReferenceFrame{T}
	function SpatialLightModulator(r, t, ref::ReferenceFrame{T}) where {T} 
		r = convert(JolabFunction2D{T}, r)
		t = convert(JolabFunction2D{T}, t)
		new{T, typeof(r), typeof(t)}(r, t, ref)
	end
end


SpatialLightModulator_slm(t, ref) = return SpatialLightModulator(ComplexF64(0), t, ref)
SpatialLightModulator_mask(t, ref) = begin
	@show "needs to be done"
	 return SpatialLightModulator(1 - t, t, ref)
end
SpatialLightModulator_reflectivemask(r, ref) = begin 
	@show "needs to be done"
	return SpatialLightModulator(r, ComplexF64(0), ref)
end

function coefficient_general(slm::SpatialLightModulator, field::FieldSpaceScalar{T}) where {T<:Real}
	checkorientation(field.ref, slm.ref) || errorToDo()
    fieldi = changereferenceframe(field, slm.ref)

	sizeI = length(field.e_SXY)
	cart = CartesianIndices(field)
	scat = get_scatteringmatrixtype(slm, field)
	for i in 1:sizeI
		scat.t₁₂.diag[i] = slm.t(field.x_X[cart[i][2]], field.y_Y[cart[i][3]])
		scat.r₁₂.diag[i] = slm.r(field.x_X[cart[i][2]], field.y_Y[cart[i][3]])
	end
	return scat
end

function get_scatteringmatrixtype(slm::SpatialLightModulator, field::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y}
	sizeI = length(field.e_SXY)
	r = Diagonal(Vector{Complex{T}}(undef, sizeI))
	t = Diagonal(Vector{Complex{T}}(undef, sizeI))
	(fieldl, fieldr) = getfields_lr(slm, field)
	return ScatteringMatrix{T, FieldSpaceScalar{T,-1,X,Y}, FieldSpaceScalar{T,1,X,Y}, Diagonal{Complex{T},Vector{Complex{T}}}, Diagonal{Complex{T},Vector{Complex{T}}}}(r, t, r, t, fieldl, fieldr)
end

function getfields_lr(slm::SpatialLightModulator, field::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y}
	fieldl = FieldSpaceScalar{T,-1,X,Y}(copy(field.x_X), copy(field.y_Y), field.e_SXY, field.λ, field.n, field.ref)
	fieldr = FieldSpaceScalar{T,1,X,Y}(copy(field.x_X), copy(field.y_Y), field.e_SXY, field.λ, field.n, field.ref)
	return (fieldl, fieldr)
end

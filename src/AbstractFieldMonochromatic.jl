changereferenceframe!(field::AbstractFieldMonochromatic, ref::ReferenceFrame) = error("Cannot do this")
samedefinitions(field1::AbstractFieldMonochromatic, field2::AbstractFieldMonochromatic) = return false
propagationmatrix(field::AbstractFieldMonochromatic, ref::ReferenceFrame) = error("Cannot do this")

function changereferenceframe(field::AbstractFieldMonochromatic, refnew::ReferenceFrame)
	tmp_field = deepcopy(field)
	changereferenceframe!(tmp_field, refnew)
	return tmp_field
end

function Base.:(+)(fielda::A, fieldb::A) where {A <: AbstractFieldMonochromatic}
	fieldc = copy(fielda)
	return add_inplace(fieldc, fieldb)
end

Base.:copy(field::AbstractFieldMonochromatic) = deepcopy(field)

dir(field::AbstractFieldMonochromatic{T,-1}) where T = -1
dir(field::AbstractFieldMonochromatic{T,1}) where T = 1 

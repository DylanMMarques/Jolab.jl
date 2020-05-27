mutable struct FieldPolichromatic
	fields::Vector{AbstractField}
end

function lightinteraction(component::AbstractOpticalComponent, field::FieldPolichromatic)
	lightinteraction.(Ref(component), field.e)
end

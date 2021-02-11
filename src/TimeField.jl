mutable struct FieldTime{T,D,F<:AbstractFieldMonochromatic{T,D}} <: AbstractField{T,D}
	fields::Vector{F}
end


function FieldTime_frompolychromatic(poly::FieldPolychromatic{T,D,F}, t)


end

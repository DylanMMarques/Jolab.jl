Base.eltype(coef::AbstractCoefficient{T}) where T = T

lightinteraction!(fieldl::AbstractFieldMonochromatic, fieldr::AbstractFieldMonochromatic, coef::AbstractCoefficient, fieldi::AbstractFieldMonochromatic) = error("Cannot do this")

getfields_lr(coef::AbstractCoefficient) = return (deepcopy(coef.fieldl), deepcopy(coef.fieldr))

function lightinteraction(coef::AbstractCoefficient, fieldi::AbstractFieldMonochromatic)
	(fieldl, fieldr) = getfields_lr(coef)

	samedefinitions(fieldi, fieldi.dir > 0 ? fieldl : fieldr)
	changereferenceframe!(fieldi, fieldi.dir > 0 ? fieldl.ref : fieldr.ref)
	lightinteraction!(fieldl, fieldr, coef, fieldi, cache)
	return (fieldl, fieldr)
end

function lightinteraction(coef::Vector{<:AbstractCoefficient}, fieldi::AbstractFieldMonochromatic)
	aux = coefficient_general(coef)
	return lightinteraction(aux, fieldi)
end

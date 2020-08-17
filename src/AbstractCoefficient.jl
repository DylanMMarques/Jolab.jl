Base.eltype(coef::AbstractCoefficient{T}) where T = T

getfields_lr(coef::AbstractCoefficient) = return (deepcopy(coef.fieldl), deepcopy(coef.fieldr))

function lightinteraction(coef::AbstractCoefficient, fieldi::AbstractFieldMonochromatic)
	(fieldl, fieldr) = getfields_lr(coef)

	changereferenceframe!(fieldi, fieldi.dir > 0 ? fieldl.ref : fieldr.ref)
	samedefinitions(fieldi, fieldi.dir > 0 ? fieldl : fieldr)
	lightinteraction!(fieldl, fieldr, coef, fieldi)
	return (fieldl, fieldr)
end

function lightinteraction(coefs::Vector{<:AbstractCoefficient}, fieldi::AbstractFieldMonochromatic)
	aux = coefficient_general(coefs)
	return lightinteraction(aux, fieldi)
end

coefficient_specific(comp::AbstractOpticalComponent, field::AbstractFieldMonochromatic) = coefficient_general(comp, field)

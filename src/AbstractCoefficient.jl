lightinteraction!(fieldl::AbstractFieldMonochromatic, fieldr::AbstractFieldMonochromatic, coef::AbstractCoefficient, fieldi::AbstractFieldMonochromatic, cache::AbstractCoefficientCache) = error("Cannot do this")

getcache(coef::AbstractCoefficient) = error("Cannot do this")
getfields_lr(coef::AbstractCoefficient) = error("Cannot do this")

function lightinteraction(coef::AbstractCoefficient, fieldi::AbstractFieldMonochromatic)
	(fieldl, fieldr) = getfields_lr(coef)
	cache = getcache(coef)

	samedefinitions(fieldi, fieldi.dir > 0 ? fieldl : fieldr)
	changereferenceframe!(fieldi, fieldi.dir > 0 ? fieldl.ref : fieldr.ref)
	lightinteraction!(fieldl, fieldr, coef, fieldi, cache)
	return (fieldl, fieldr)
end

function lightinteraction(coef::Vector{<:AbstractCoefficient}, fieldi::AbstractFieldMonochromatic)
	aux = coefficient_general(coef)
	return lightinteraction(aux, fieldi)
end

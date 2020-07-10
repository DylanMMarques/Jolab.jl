Base.eltype(coef::AbstractCoefficient{T}) where T = T

getfields_lr(coef::AbstractCoefficient) = return (deepcopy(coef.fieldl), deepcopy(coef.fieldr))

function lightinteraction(coef::AbstractCoefficient, fieldi::AbstractFieldMonochromatic)
	(fieldl, fieldr) = getfields_lr(coef)

	samedefinitions(fieldi, fieldi.dir > 0 ? fieldl : fieldr)
	changereferenceframe!(fieldi, fieldi.dir > 0 ? fieldl.ref : fieldr.ref)
	lightinteraction!(fieldl, fieldr, coef, fieldi)
	return (fieldl, fieldr)
end

function lightinteraction(coef::Vector{<:AbstractCoefficient}, fieldi::AbstractFieldMonochromatic)
	aux = coefficient_general(coef)
	return lightinteraction(aux, fieldi)
end

coefficient_specific(comp::AbstractOpticalComponent, field::AbstractFieldMonochromatic) = coefficient_general(comp, field)

function coefficient_general(comps::AbstractVector{<:AbstractOpticalComponent{T}}, fieldi::AbstractFieldMonochromatic) where T
    coef = Vector{ScatteringMatrix{T}}(undef, length(comps))
    fieldaux = fieldi
    if fieldi.dir > 0
        for i in eachindex(comps)
            coef[i] = coefficient_general(comps[i], fieldaux)
            (tmp, fieldaux) = getfields_lr(coef[i])
        end
    else
        sizeA = length(comps)
        for i in eachindex(comps)
            coef[sizeA - i + 1] = coefficient_general(comps[sizeA-i + 1], fieldaux)
            (fieldaux, tmp) = getfields_lr(coef[sizeA-i+1])
        end
    end
    return coefficient_general(coef)
end

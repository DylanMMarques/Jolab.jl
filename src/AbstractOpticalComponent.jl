function lightinteraction(comp::AbstractOpticalComponent{T}, fieldi::AbstractFieldMonochromatic) where T
    coef = coefficient_specific(comp, fieldi)
    return lightinteraction(coef, fieldi)
end

function lightinteraction(comps::AbstractVector{<:AbstractOpticalComponent{T}}, fieldi::AbstractFieldMonochromatic) where T
    coef = Vector{ScaterringMatrix{T}}(undef, length(comps))
    fieldaux = fieldi
    if fieldi.dir > 0
        for i in eachindex(comps)
            coef[i] = coefficient_geral(comps[i], fieldaux)
            fieldaux = coef[i].fieldr
        end
    else
        sizeA = length(comps)
        for i in eachindex(comps)
            coef[sizeA - i] = coefficient_geral(comps[sizeA-i], fieldaux)
            fieldaux = coef[sizeA-i].fieldl
        end
    end
    coef = mergeorientated_propagationcoefficient(coef)
    return lightinteraction(coef, fieldi)
end

function lightinteraction(comp::AbstractOpticalComponent{T}, fieldi::AbstractFieldMonochromatic) where T
    coef = coefficient_specific(comp, fieldi)
    return lightinteraction(coef, fieldi)
end

function coefficient_general(comps::AbstractVector{<:AbstractOpticalComponent{T}}, fieldi::AbstractFieldMonochromatic) where T
    coef = Vector{ScaterringMatrix{T}}(undef, length(comps))
    fieldaux = fieldi
    if fieldi.dir > 0
        for i in eachindex(comps)
            coef[i] = coefficient_general(comps[i], fieldaux)
            fieldaux = coef[i].fieldr
        end
    else
        sizeA = length(comps)
        for i in eachindex(comps)
            coef[sizeA - i + 1] = coefficient_general(comps[sizeA-i + 1], fieldaux)
            fieldaux = coef[sizeA-i+1].fieldl
        end
    end
    coef = mergeorientated_propagationcoefficient(coef)
end

function lightinteraction(comps::AbstractVector{<:AbstractOpticalComponent{T}}, fieldi::AbstractFieldMonochromatic) where T
    coef = coefficient_general(comps, fieldi)
    return lightinteraction(coef, fieldi)
end

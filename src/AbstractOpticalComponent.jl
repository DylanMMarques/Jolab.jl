function lightinteraction(comp::AbstractOpticalComponent{T}, fieldi::AbstractFieldMonochromatic) where T
    coef = coefficient_specific(comp, fieldi)
    return lightinteraction(coef, fieldi)
end

coefficient_specific(comp::AbstractOpticalComponent, field::AbstractFieldMonochromatic) = coefficient_general(comp, field)

function coefficient_general(comps::AbstractVector{<:AbstractOpticalComponent{T}}, fieldi::AbstractFieldMonochromatic) where T
    coef = Vector{ScatteringMatrix{T}}(undef, length(comps))
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
    coef = coefficient_general(coef)
end

function lightinteraction(comps::AbstractVector{<:AbstractOpticalComponent{T}}, fieldi::AbstractFieldMonochromatic) where T
    coef = coefficient_general(comps, fieldi)
    return lightinteraction(coef, fieldi)
end

rotatestructure!(comp::AbstractOpticalComponent, ref_ref::ReferenceFrame, θ::Real, ϕ::Real) = rotatereferenceframe!(ref_ref, comp.ref, θ, ϕ)
function rotatestructure(comp::AbstractOpticalComponent, ref_ref::ReferenceFrame, θ::Real, ϕ::Real)
	aux = deepcopy(comp)
	rotatestructure!(aux, ref_ref, θ, ϕ)
	return aux
end

translatestructure!(comp::AbstractOpticalComponent, x::Real, y::Real, z::Real) = translatereferenceframe!(comp.ref, x, y, z)
function translatestructure(comp::AbstractOpticalComponent, x::Real, y::Real, z::Real)
	aux = deepcopy(comp)
	translatestructure!(aux, x, y, z)
	return aux
end

function rotatestructure(comps::AbstractVector{T}, ref_ref::ReferenceFrame, θ::Real, ϕ::Real) where {T<:AbstractOpticalComponent}
	out = Vector{T}(undef, length(comps))
	@inbounds for i in eachindex(comps)
		out[i] = rotatestructure(comps[i], ref_ref, θ, ϕ)
	end
	return out
end

function translatestructure(comps::AbstractVector{T}, x::Real, y::Real, z::Real) where {T<:AbstractOpticalComponent}
	out = Vector{T}(undef, length(comps))
	@inbounds for i in eachindex(comps)
		out[i] = translatestructure(comps[i], x, y, z)
	end
	return out
end

ref1(comp::AbstractOpticalComponent) = comp.ref
ref2(comp::AbstractOpticalComponent) = comp.ref

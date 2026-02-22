struct Beam{T, D, M<:StructArray{<:AbstractFieldMode{T,D}}} <: AbstractField{T,D}
    modes::M
    function Beam(modes::M) where {M<:StructArray{<:AbstractFieldMode{T,D}}} where {T,D}
        new{T, D, M}(modes)
    end
end

function Base.isapprox(a::Beam{T1,D1,<:StructArray{M1}}, b::Beam{T2,D2,<:StructArray{M2}}; kwargs...) where {T1, T2, D1, D2, M1, M2}
    D1 == D2 || return false
    same_mode_type(M1, M2) || return false
    fields = fieldnames(M1)
    for field in fields
        all(isapprox.(component(a.modes, field), component(b.modes, field); kwargs...)) || return false
    end
    return true
end

function intensity(beam::Beam)
    sum(abs2, beam.modes.e)
end

function fill_zeros!(beam::Beam)
    fill!(beam.modes.e, 0)
end

function _unchecked_add!(beam::Beam, other::Beam)
    beam.modes.e .+= other.modes.e
end

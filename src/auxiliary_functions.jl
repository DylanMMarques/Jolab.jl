reverse_if_backward(::Type{Forward}, x) = x
reverse_if_backward(::Type{Backward}, x) = reverse(x)

import Base.!
Base.length(beam::Beam) = length(beam.modes)
Base.size(beam::Beam) = size(beam.modes)
(!)(::Type{Forward}) = Backward
(!)(::Type{Backward}) = Forward

function check_same_definition(mode1::M1, mode2::M2) where {M2<:AbstractFieldMode{<:Any,D2}, M1<:AbstractFieldMode{<:Any,D1}} where {D1, D2}
    same_mode_type(M1, M2) || return false
    D1 == D2 || return false
    fields = fieldnames(M1)
    for field in Iterators.filter(!=(:e), fields)
        getfield(mode1, field) ≈ getfield(mode2, field) || return false
    end
    return true
end

function check_same_definition(beam1::Beam{T1,D1,<:StructArray{M1}}, beam2::Beam{T2,D2,<:StructArray{M2}}) where {T1,T2,D1,D2,M1,M2}
    same_mode_type(M1, M2) || return false
    D1 == D2 || return false 
    fields = fieldnames(M1)
    for field in Iterators.filter(!=(:e), fields)
        all(component(beam1.modes, field) .≈ component(beam2.modes, field)) || return false
    end
    return true
end

same_mode_type(::Type{<:PlaneWaveScalar}, ::Type{<:PlaneWaveScalar}) = true
same_mode_type(::Type{<:PlaneWaveVectorial}, ::Type{<:PlaneWaveVectorial}) = true
same_mode_type(::Type{<:AbstractFieldMode}, ::Type{<:AbstractFieldMode}) = false

function find_local_minima(::Type{T}, f, x, initial_size)
    local_minima = Vector{T}(undef, initial_size)
    ind_zeros = 1
    value = f(x[1])
    next_value = f(x[2])
    @inbounds for ix in view(x, 3:length(x))
        prev_value = value
        value = next_value
        next_value = f(ix)
        if value < prev_value && value < next_value
            local_minima[ind_zeros] = ix - step(x)
            ind_zeros += 1
        end
    end
    local_minima[1:(ind_zeros-1)]
end
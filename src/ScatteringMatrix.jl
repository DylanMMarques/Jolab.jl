struct ScatteringMatrix{T, F2<:Beam{T, Backward}, F3<:Beam{T, Forward}, R1<:AbstractMatrix{RealOrComplex{T}}, T1<:AbstractMatrix{RealOrComplex{T}}, F1<:Beam{T}}
    field_b::F2
    field_f::F3
    mat_itob::R1
    mat_itof::T1
    field_i::F1
    function ScatteringMatrix(::Type{T}, field_b::F2, field_f::F3, i_to_b::R1, i_to_f::T1, fieldi::F1) where {T, R1, T1, F1<:Beam{T}, F2, F3}
        number_modes_field_i = number_modes(fieldi)
        number_modes_field_b = number_modes(field_b)
        number_modes_field_f = number_modes(field_f)
        @argcheck size(i_to_f) == (number_modes_field_f, number_modes_field_i) DimensionMismatch
        @argcheck size(i_to_b) == (number_modes_field_b, number_modes_field_i) DimensionMismatch
        new{T, R1, T1, F1, F2, F3}(field_b, field_f, i_to_b, i_to_f, fieldi)
    end
end

ScatteringMatrix(field_b, field_t, i_to_b, i_to_f, field_i) = ScatteringMatrix(Float64, field_b, field_t, i_to_b, i_to_f, field_i)

function light_interaction!(backwardfield::Beam, forwardfield::Beam, sm::ScatteringMatrix, field_i::Beam{T, D}) where {T, D}
    @argcheck check_same_definition(field_i, sm.field_i) ArgumentError
    (smfield_b, smfield_f) = reverse_if_backward(D, (sm.field_r, field_t))
        
    @argcheck check_same_definition(forwardfield, smfield_f) ArgumentError
    @argcheck check_same_definition(backwardfield, smfield_b) ArgumentError

    (rfield, tfield) = reverse_if_backward(D, (backwardfield, forwardfield))
    mul!(vec(rfield.modes.e), sm.r, vec(field_i.modes.e))
    mul!(vec(tfield.modes.e), sm.t, vec(field_i.modes.e))
    return (backwardfield, forwardfield)
end

function light_interaction(sm::ScatteringMatrix, field_i::Beam{T, D}) where {T, D}
    rfield = deepcopy(sm.field_r)
    tfield = deepcopy(sm.field_t)

    (backwardfield, forwardfield) = reverse_if_backward(D, (rfield, tfield))

    light_interaction!(backwardfield, forwardfield, sm, field_i)
end

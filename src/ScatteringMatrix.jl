struct ScatteringMatrix{T, F2<:Beam{T, Backward}, F3<:Beam{T, Forward}, R1<:AbstractMatrix{<:RealOrComplex{T}}, T1<:AbstractMatrix{<:RealOrComplex{T}}, F1<:Beam{T}}
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
        new{T, F2, F3, R1, T1, F1}(field_b, field_f, i_to_b, i_to_f, fieldi)
    end
end

ScatteringMatrix(field_b, field_t, i_to_b, i_to_f, field_i) = ScatteringMatrix(Float64, field_b, field_t, i_to_b, i_to_f, field_i)

function _light_interaction!(backwardfield::Beam, forwardfield::Beam, sm::ScatteringMatrix, field_i::Beam{T, D}) where {T, D}
    mul!(vec(backwardfield.modes.e), sm.mat_itob, vec(field_i.modes.e))
    mul!(vec(forwardfield.modes.e), sm.mat_itof, vec(field_i.modes.e))
    return (backwardfield, forwardfield)
end

function light_interaction(sm::ScatteringMatrix, field_i)
    @argcheck check_same_definition(sm.field_i, field_i) ArgumentError
    _light_interaction!(deepcopy(sm.field_b), deepcopy(sm.field_f), sm, field_i)
end

function ScatteringMatrix(field_b, field_f, comp, field_i)
    @argcheck check_input_field(comp, field_i) ArgumentError
    @argcheck check_output_fields(field_b, field_f, comp, field_i) ArgumentError
    _ScatteringMatrix(field_b, field_f, comp, field_i)
end

function ScatteringMatrix(comp, field_i)
    (field_b, field_f) = forward_backward_field(comp, field_i)
    _ScatteringMatrix(field_b, field_f, comp, field_i)
end

import Base: *

(*)(sm::ScatteringMatrix, field_i) = light_interaction(sm, field_i)
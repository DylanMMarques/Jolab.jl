struct ScatteringMatrix{T, F2<:AbstractField{T, Backward}, F3<:AbstractField{T, Forward}, R1<:AbstractMatrix{<:RealOrComplex{T}}, T1<:AbstractMatrix{<:RealOrComplex{T}}, F1<:AbstractField{T}}
    field_b::F2
    field_f::F3
    mat_itob::R1
    mat_itof::T1
    field_i::F1
    function ScatteringMatrix(::Type{T}, field_b::F2, field_f::F3, i_to_b::R1, i_to_f::T1, fieldi::F1) where {T, R1, T1, F1<:AbstractField{T}, F2, F3}
        number_modes_field_i = length(fieldi.e)
        number_modes_field_b = length(field_b.e)
        number_modes_field_f = length(field_f.e)
        @argcheck size(i_to_f) == (number_modes_field_f, number_modes_field_i) DimensionMismatch
        @argcheck size(i_to_b) == (number_modes_field_b, number_modes_field_i) DimensionMismatch
        new{T, F2, F3, R1, T1, F1}(field_b, field_f, i_to_b, i_to_f, fieldi)
    end
end

ScatteringMatrix(field_b, field_t, i_to_b, i_to_f, field_i) = ScatteringMatrix(Float64, field_b, field_t, i_to_b, i_to_f, field_i)

function _light_interaction!(backwardfield::MeshedBeam, forwardfield::MeshedBeam, sm::ScatteringMatrix, field_i::MeshedBeam{T, D}) where {T, D}
    mul!(vec(backwardfield.e), sm.mat_itob, vec(field_i.e))
    mul!(vec(forwardfield.e), sm.mat_itof, vec(field_i.e))
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
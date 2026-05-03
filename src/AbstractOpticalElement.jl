function light_interaction(
    comps::NTuple{N,AbstractOpticalElement},
    field::MeshedBeam{T,D},
) where {T,D,N}
    scat = ScatteringMatrix(comps, field)
    ind = D == Forward ? 2 : 1
    light_interaction(scat, field)
end

function ScatteringMatrix(
    comps::NTuple{N,AbstractOpticalElement{T}},
    field_forward::MeshedBeam{T,D},
) where {N,D,T}
    function fg(scat12, comp)
        prop = Propagation(
            T,
            (scat12[2].field_f.frame, first(comp.frames)),
            scat12[2].field_f.medium,
        )
        scat23 = ScatteringMatrix(prop, scat12[2].field_f)
        scat32 = ScatteringMatrix(prop, reverse_direction(scat23.field_f))
        scat13 = recursive_light_interaction_inverse(scat12, (scat32, scat23))
        scat34 = ScatteringMatrix(comp, scat13[2].field_f)
        scat43 = ScatteringMatrix(comp, reverse_direction(scat34.field_f))
        recursive_light_interaction_inverse(scat13, (scat43, scat34))
    end
    scat12 = ScatteringMatrix(comps[1], field_forward)
    scat21 = ScatteringMatrix(comps[1], reverse_direction(scat12.field_f))
    res = mapreduce(identity, fg, comps[2:end], init = (scat21, scat12))
    return res[D == Forward ? 2 : 1]
end

function reverse_direction(
    field::MeshedBeam{T,D,C,P},
) where {T,D<:Union{Forward,Backward},C,P}
    return MeshedBeam{T,!D,C,P}(field.mesh, field.e, field.medium, field.frame)
end

function check_mergable_scattering_matrices(mat_1, mat_2)
    code = zero(UInt64)
    mat_12, mat_21 = mat_1[2], mat_1[1]
    mat_23, mat_32 = mat_2[2], mat_2[1]
    check_same_definition(mat_12.field_f, mat_23.field_i) ||
        (code |= 1 << INVALID_BEAM_SHAPE)
    check_same_definition(mat_21.field_i, mat_32.field_b) ||
        (code |= 1 << INVALID_BEAM_SHAPE)
    code
end

function recursive_light_interaction_inverse(
    mat_1::Tuple{ScatteringMatrix{T1},ScatteringMatrix{T2}},
    mat_2::Tuple{ScatteringMatrix{T3},ScatteringMatrix{T4}},
) where {T1,T2,T3,T4}
    code = check_mergable_scattering_matrices(mat_1, mat_2)
    code == 0 || throw_error_msg(code)

    mat_1_21, mat_1_12 = mat_1
    mat_2_21, mat_2_12 = mat_2

    r12 = mat_1_12.mat_itob
    t12 = mat_1_12.mat_itof
    r21 = mat_1_21.mat_itof
    t21 = mat_1_21.mat_itob

    r23 = mat_2_12.mat_itob
    t23 = mat_2_12.mat_itof
    r32 = mat_2_21.mat_itof
    t32 = mat_2_21.mat_itob

    aux_1 = inv(I - r23 * r21)
    r_13 = r12 + t21 * aux_1 * r23 * t12
    t_31 = t21 * aux_1 * t32

    aux_2 = inv(I - r21 * r23)
    r_31 = r32 + t23 * aux_2 * r21 * t32
    t_13 = t23 * aux_2 * t12

    T = promote_type(T1, T2, T3, T4)
    scat_13 = ScatteringMatrix(
        T,
        mat_1_12.field_b,
        mat_2_12.field_f,
        r_13,
        t_13,
        mat_1_12.field_i,
    )
    scat_31 = ScatteringMatrix(
        T,
        mat_1_12.field_b,
        mat_2_12.field_f,
        t_31,
        r_31,
        mat_2_21.field_i,
    )
    return (scat_31, scat_13)
end

function solver(comp::AbstractOpticalElement, field_i)
    comp
end

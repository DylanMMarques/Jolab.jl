function light_interaction(comps::NTuple{N,AbstractOpticalElement}, field::MeshedBeam{T,D}) where {T,D,N}
    scat = ScatteringMatrix(comps, field)
    ind = D == Forward ? 2 : 1
    light_interaction(scat[ind], field)
end

function ScatteringMatrix(comps::NTuple{N,AbstractOpticalElement{T}}, field_forward::MeshedBeam{T,D}) where {N,D,T}
    function f(comp, field_i)
        scat_forward = ScatteringMatrix(comp, field_i)
        scat_backward = ScatteringMatrix(comp, reverse_direction(scat_forward.field_f))
        (scat_backward, scat_forward)
    end
    scats = map(comp_i -> f(comp_i, field_forward), comps)
    
    function g(scat_1, scat_2)
        frame_f = scat_1[2].field_f.frame 
        frame_f_i = scat_2[2].field_i.frame
        prop = Propagation(T, (frame_f, frame_f_i), scat_1[2].field_f.medium)
        scat_f = ScatteringMatrix(prop, scat_1[2].field_f)
        
        frame_b = scat_2[1].field_b.frame
        frame_b_i = scat_1[1].field_i.frame
        @show frame_b == frame_f_i
        @show frame_b_i == frame_f
        scat_b = ScatteringMatrix(prop, scat_2[1].field_b)
        (scat_b, scat_f)
    end
    scat_props = g.(scats[1:end-1], scats[2:end])
    return scats
    scats_2 = map(recursive_light_interaction_inverse, scat_props, scats[2:end])
    res = reduce(recursive_light_interaction_inverse, scats_2, init = first(scats))
    res[D == Forward ? 2 : 1]
end

function reverse_direction(field::MeshedBeam{T,D,C}) where {T,D<:Union{Forward, Backward},C}
    return MeshedBeam{T,!D,C}(field.mesh, field.e, field.medium, field.frame)
end

function check_mergable_scattering_matrices(mat_1, mat_2)
    code = zero(UInt64)
    isapprox(mat_1[2].field_f, mat_2[2].field_i) || (error(); code |= 1 << INVALID_BEAM_SHAPE)
    isapprox(mat_2[1].field_b, mat_1[1].field_i) || (error(); code |= 1 << INVALID_BEAM_SHAPE)
    code
end

function recursive_light_interaction_inverse(mat_1::Tuple{ScatteringMatrix{T1}, ScatteringMatrix{T2}}, mat_2::Tuple{ScatteringMatrix{T3}, ScatteringMatrix{T4}}) where {T1, T2, T3, T4}
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
    scat_13 = ScatteringMatrix(T, mat_1_12.field_b, mat_2_12.field_f, r_13, t_13, mat_1_12.field_i)
    scat_31 = ScatteringMatrix(T, mat_1_12.field_b, mat_2_12.field_f, t_31, r_31, mat_2_21.field_i)
    return (scat_31, scat_13)
end



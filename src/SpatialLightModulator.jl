struct SpatialLightModulator{T,C,R,F} <: AbstractOpticalElement{T}
    r::R
    t::F
    frame::ReferenceFrame{T}
    media::Tuple{Medium{T,T}, Medium{T,T}}
    function SpatialLightModulator{T,C}(r::R, t::F, frame, media) where {T,C,R,F} 
        new{T,C,R,F}(r, t, frame, media)
    end
end

function SpatialLightModulator_aperture(::Type{T}, radius, media, frame) where T
    slm_r = (r, theta, lambda) -> false
    slm_t = (r, theta, lambda) -> radius >= r
    SpatialLightModulator{T,R_θ_λ}(slm_r, slm_t, frame, media)
end
function SpatialLightModulator_diskaperture(::Type{T}, internal_radius, outer_radius, media, frame) where T
    slm_r = (r, theta, lambda) -> false
    slm_t = (r, theta, lambda) -> internal_radius <= r <= outer_radius
    SpatialLightModulator{T,R_θ_λ}(slm_r, slm_t, frame, media)
end

function reflection_transmission_function(slm::SpatialLightModulator{<:Any, C}, ::Type{C}) where {C}
    slm.r, slm.t
end,
function reflection_transmission_function(slm::SpatialLightModulator{<:Any, C}, ::Type{C1}) where {C,C1}
    coord_trans(c1,c2,c3) = C(C1(c1,c2,c3))
    new_t(c1, c2, c3) = begin
        new_coord = coord_trans(c1, c2, c3)
        slm.t(new_coord[1], new_coord[2], new_coord[3])
    end
    new_r(c1, c2, c3) = begin
        new_coord = coord_trans(c1, c2, c3)
        slm.r(new_coord[1], new_coord[2], new_coord[3])
    end
    (new_r, new_t)
end


function _light_interaction!(field_b::MeshedBeam{<:Any,Backward}, field_f::MeshedBeam{<:Any,Forward}, slm::SpatialLightModulator, field_i::MeshedBeam{T,D,C,P}) where {T,D,C<:SpatialCoords,P}
    (field_r, field_t) = reverse_if_backward(D, (field_b, field_f))
    
    slm_r, slm_t = reflection_transmission_function(slm, C)

    tmp_vec = StructVector{Tuple{Complex{T}, Complex{T}}}((vec(field_r.e), vec(field_t.e)))
    broadcast!(tmp_vec, eachindex(field_i.e)) do i
        coord = centroid(field_i.mesh, i)
        (slm_r(coord[1], coord[2], coord[3]) * field_i.e[i], slm_t(coord[1], coord[2], coord[3]) * field_i.e[i])
    end

    (field_b, field_f)
end

function check_input_field(slm::SpatialLightModulator, field::MeshedBeam{T, D, C, P}) where {T,D,C<:SpatialCoords,P}
    code = zero(UInt64)
    slm.frame ≈ field.frame || (code |= 1 << INVALID_FRAME)
    medium = D == Forward ? first(slm.media) : last(slm.media)
    medium ≈ field.medium || (code |= 1 << INVALID_MEDIUM)
    code
end

function forward_backward_field(slm::SpatialLightModulator, field_i::MeshedBeam{T,D,C,P}) where {T,D,C<:SpatialCoords,P}
    e_b = similar(field_i.e, Complex{T}, size(field_i.mesh))
    e_f = similar(field_i.e, Complex{T}, size(field_i.mesh))

    field_f = MeshedBeam{T, Forward, C, P}(field_i.mesh, e_f, field_i.medium, field_i.frame)
    field_b = MeshedBeam{T, Backward, C, P}(field_i.mesh, e_b, field_i.medium, field_i.frame)
    reverse_if_backward(D, (field_b, field_f))
end

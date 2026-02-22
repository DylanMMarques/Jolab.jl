struct Propagation{T, R<:ReferenceFrame{T}, M<:Medium{T}} <: AbstractOpticalElement{T}
    frames::NTuple{2, R}
    medium::M
    function Propagation(::Type{T}, frames::NTuple{2, R}, medium::M) where {T, R<:ReferenceFrame{T}, M<:Medium{T}}
        new{T, R, M}(frames, medium)
    end
end
Propagation(frames, medium) = Propagation(Float64, frames, medium)


function _light_interaction!(field_b, field_f, prop::Propagation, field_i::MeshedAngularSpectrum{T,D,C}) where {T,D,C}
    frame = D == Forward ? prop.frames[2] : prop.frames[1]

    field_r, field_t = reverse_if_backward(D, (field_b, field_f))
    fill_zeros!(field_r)

    rpos = delta_pos_referenceframe(field_i.frame, frame.origin)

    f(ind) = t(prop, D, prop.medium.n, rpos, C(centroid(field_i.mesh, ind)))
    vec(field_t.e) .= f.(eachindex(field_i.e)) .* vec(field_i.e)
  
    (field_b, field_f)
end

"""
Calculates the cartesians distance between the origin of a reference frame and a new origin. The new origin is given in the coordinate system of the input reference frame.
"""
function delta_pos_referenceframe(frame::ReferenceFrame{T}, new_origin::X_Y_Z) where T
    Δpos = new_origin - frame.origin
    rot = _RotXYZ(frame.direction)
    X_Y_Z(inv(rot) * Δpos)
end

function t(prop::Propagation{T}, ::Type{D}, n, rΔpos::X_Y_Z, coord::NSX_NSY_λ) where {D,T}
    (nsx, nsy, λ) = coord
    positive_nsz = √(complex(n^2 - nsx^2 - nsy^2))
    nsz = D == Forward ? positive_nsz : -positive_nsz
    exp(im * 2T(π) / λ * dot(rΔpos, (nsx, nsy, nsz)))
end
t(prop::Propagation, ::Type{D}, medium, rΔpos::X_Y_Z, coord::NSR_NSθ_λ) where D = t(prop, D, medium, rΔpos, NSX_NSY_λ(coord))

function _ScatteringMatrix(field_b, field_f, prop::Propagation, field_i::MeshedAngularSpectrum{T,D,C}) where {T,D,C<:AngularSpectrumCoords}
    t_vec = similar(field_i.e, Complex{T}, length(field_i.e))
    r = Zeros(T, (length(field_i.e), length(field_i.e)))
    frame = D == Forward ? prop.frames[2] : prop.frames[1]
    rpos = delta_pos_referenceframe(field_i.frame, frame.origin)
    f(ind) = t(prop, D, prop.medium.n, rpos, C(centroid(field_i.mesh, ind)))
    t_vec .= f.(eachindex(field_i.e))
    
    (mat_i_to_b, mat_i_to_f) = reverse_if_backward(D, (r, Diagonal(vec(t_vec))))
    ScatteringMatrix(T, field_b, field_f, mat_i_to_b, mat_i_to_f, field_i)
end

function check_input_field(prop::Propagation, field::MeshedAngularSpectrum)
    code = zero(UInt64)
    field.medium ≈ prop.medium || (code |= 1 << INVALID_MEDIUM)
    code
end

function check_input_field(prop::Propagation, field::MeshedBeam{<:Any, D}) where D
    code = zero(UInt64)
    field.medium ≈ prop.medium || (code |= 1 << INVALID_MEDIUM)
    frame = D == Forward ? prop.frames[2] : prop.frames[1]
    field.frame ≈ frame || (code |= 1 << INVALID_CANNOT_BE_TRANSLATED)
    code
end

function forward_backward_field(prop::Propagation, field_i::MeshedAngularSpectrum{T,D,C,P}) where {T,D,C,P}
    frame = D == Forward ? prop.frames[2] : prop.frames[1]
    field_r = MeshedAngularSpectrum{T, !D, C,P}(field_i.mesh, Zeros(T, size(field_i.e)), field_i.medium, field_i.frame)
    field_t = MeshedAngularSpectrum{T, D, C,P}(field_i.mesh, similar(field_i.e, Complex{T}), field_i.medium, frame)
    reverse_if_backward(D, (field_r, field_t))
end

function _ScatteringMatrix(field_b, field_f, prop::Propagation, field_i::MeshedBeam{T, D}) where {T,D}
    # Only called if same reference frame. Defined for inverse scattering matrix solving
    t_vec = Ones(T, length(field_i.e))
    r = Zeros(T, (length(field_i.e), length(field_i.e)))
    (mat_i_to_b, mat_i_to_f) = reverse_if_backward(D, (r, Diagonal(vec(t_vec))))
    ScatteringMatrix(T, field_b, field_f, mat_i_to_b, mat_i_to_f, field_i)
end

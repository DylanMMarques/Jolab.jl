struct Lens{T, F, N, M1<:Medium{T}, M2<:Medium{T}} <: AbstractOpticalElement{T}
    focal_length::F
    numerical_aperture::N
    mat::Tuple{M1, M2}
    frames::Tuple{ReferenceFrame{T}, ReferenceFrame{T}}
    function Lens(::Type{T}, focal_length::F, numerical_aperture::N, media::Tuple{M1,M2}, frame) where {T,F,N,M1,M2}

        any(is_complex_medium, media) && throw(ArgumentError("A lens cannot be in a medium with absortion. Use media with real refractive index."))

        frame_1 = frame.origin + _RotXYZ(frame.direction) * X_Y_Z(zero(T), zero(T), -focal_length)
        frame_2 = frame.origin + _RotXYZ(frame.direction) * X_Y_Z(zero(T), zero(T), focal_length)
        new{T,F,N,M1,M2}(focal_length, numerical_aperture, media, map(o -> ReferenceFrame(o, frame.direction), (frame_1, frame_2)))
    end
end
export Lens

Lens(focal_length, numerical_aperture, media, frame) = Lens(Float64, focal_length, numerical_aperture, media, frame)

@inbounds function _light_interaction!(field_b::MeshedBeam{<:Any,Backward}, field_f::MeshedBeam{<:Any,Forward}, lens::Lens, field_i::MeshedBeam{T,D,C,PolarizationScalar}) where {T,D,C}
    (field_r, field_t) = reverse_if_backward(D, (field_b, field_f))
    fill!(field_r.e, 0)

    function t_aux(ind)
        t(lens, C(centroid(field_i.mesh, ind))) * field_i.e[ind]
    end

    vec(field_t.e) .= t_aux.(eachindex(field_i.e))

    (field_b, field_f)
end

@inbounds function _light_interaction!(field_b::MeshedBeam{<:Any,Backward}, field_f::MeshedBeam{<:Any,Forward}, lens::Lens, field_i::MeshedBeam{T,D,C, PolarizationXYZ}) where {T,D,C<:SpatialCoords}
    (field_r, field_t) = reverse_if_backward(D, (field_b, field_f))
    fill!(field_r.e, 0)
    
    field_i_x = view(field_i.e, :, :, 1)
    field_i_y = view(field_i.e, :, :, 2)
    field_i_z = view(field_i.e, :, :, 3)
    field_t_s = view(field_t.e, :, :, 1)
    field_t_p = view(field_t.e, :, :, 2)
    function t_aux(ind)

        r_pupil = R_θ_λ(centroid(field_i.mesh, ind))
        t_coef = t(lens, r_pupil)

        n_s = unit_vector(R_θ_Z(r_pupil[1], r_pupil[2], zero(T)), Val(:θ))
        n_p = unit_vector(R_θ_Z(r_pupil[1], r_pupil[2], zero(T)), Val(:R))
    
        e_i_vector = SVector(field_i_x[ind], field_i_y[ind], field_i_z[ind])
        comp_s = dot(n_s, e_i_vector)
        comp_p = dot(n_p, e_i_vector)
        field_t_s[ind] = comp_s * t_coef
        field_t_p[ind] = comp_p * t_coef
    end
    map(t_aux, eachindex(field_i.mesh))
    (field_b, field_f)
end

function t(lens::Lens{T}, coord::NSR_NSθ_λ) where T<:AbstractFloat
    (nsr, θ, λ) = coord
    cosθ² = 1 - (nsr / real(lens.mat[1].n))^2
    if cosθ² < 0 # Evasnecent waves
    	zero(Complex{T})
    else # Plane waves
    	if (1 - cosθ² > lens.numerical_aperture^2) # above the lens NA
    	    zero(Complex{T})
    	else
            1 / (T(lens.focal_length * cosθ²^(1/4) / (2π / λ * 2π)) / im) # might be wrong
        end
    end
end,
function t(lens::Lens{T}, coord::R_θ_λ) where T<:AbstractFloat
    (r, θ, λ) = coord
    cosθ² = 1 - (r / lens.focal_length)^2
    if cosθ² < 0 # Evasnecent waves
    	zero(Complex{T})
    else # Plane waves
    	if (1 - cosθ² > lens.numerical_aperture^2) # above the lens NA
    	    zero(Complex{T})
    	else
    	    T(lens.focal_length * cosθ²^(1/4) / (2π / λ * 2π)) / im # Needs checking
    	end
    end
end,
function t(lens::Lens, coord::X_Y_λ)
    t(lens, R_θ_λ(coord))
end,
function t(lens::Lens, coord::NSX_NSY_λ)
    t(lens, NSR_NSθ_λ(coord))
end,

function _ScatteringMatrix(field_b::MeshedBeam, field_f::MeshedBeam, comp::Lens, field_i::MeshedBeam{T,D,C,P}) where {T,D, C,P<:PolarizationScalar}
    (field_r, field_t) = reverse_if_backward(D, (field_b, field_f))
    r = Zeros(length(field_r.e), length(field_i.e))

    t_dia = similar(field_i.e, Complex{T}, length(field_i.e))

    t_aux(ind) = t(comp, C(centroid(field_i.mesh, ind)))

    vec(t_dia) .= t_aux.(eachindex(field_i.mesh))
    
    (mat_i_to_b, mat_i_to_f) = reverse_if_backward(D, (r, Diagonal(vec(t_dia))))
    ScatteringMatrix(T, field_b, field_f, mat_i_to_b, mat_i_to_f, field_i)
end

function forward_backward_field(lens::Lens, field_i::MeshedBeam{T,D,C,P}) where {T,D,C,P}
    P_out = get_polarization_type(Lens, P)

    e_t = similar(field_i.e, Complex{T}, (size(field_i.mesh)[1:2]..., number_components(P_out)))
    e_r = Zeros(T, size(field_i.e))

    f = lens.focal_length
    n = D == Forward ? first(lens.mat).n : last(lens.mat).n 
    
    mesh = if C <: AngularSpectrumCoords
        if C <: NSR_NSθ_λ
            Scale((T(f / n), T(1), T(1)))(change_coordinate_type(R_θ_λ, field_i.mesh))
        else
            Scale((T(f / n), T(f / n), T(1)))(change_coordinate_type(X_Y_λ, field_i.mesh))
        end
    else
        if C <: R_θ_λ
            Scale((T(n / f), T(1), T(1)))(change_coordinate_type(NSR_NSθ_λ, field_i.mesh))
        else
            Scale((-T(n / f), -T(n / f), T(1)))(change_coordinate_type(NSX_NSY_λ, field_i.mesh))
        end
    end

    frame_t = D == Forward ? last(lens.frames) : first(lens.frames)
    (dir_r, dir_t) = reverse_if_backward(D, (Backward, Forward))
    field_r = MeshedBeam{T, dir_r, C, P}(field_i.mesh, e_r, field_i.medium, field_i.frame)
    P_out = get_polarization_type(Lens, P)
    field_t = MeshedBeam{T, dir_t, get_transmitted_coord_type(Lens, C), P_out}(mesh, e_t, lens.mat[2], frame_t)
    reverse_if_backward(D, (field_r, field_t))
end

get_polarization_type(::Type{Lens}, P::Type{PolarizationXYZ}) = PolarizationSP
get_polarization_type(::Type{Lens}, P::Type{PolarizationSP}) = PolarizationXYZ
get_polarization_type(::Type{Lens}, P::Type{PolarizationScalar}) = PolarizationScalar

get_transmitted_coord_type(::Type{Lens}, ::Type{NSX_NSY_λ}) = X_Y_λ
get_transmitted_coord_type(::Type{Lens}, ::Type{NSR_NSθ_λ}) = R_θ_λ
get_transmitted_coord_type(::Type{Lens}, ::Type{X_Y_λ}) = NSX_NSY_λ
get_transmitted_coord_type(::Type{Lens}, ::Type{R_θ_λ}) = NSR_NSθ_λ

function check_output_fields(field_b, field_f, comp::Lens, field_i::MeshedAngularSpectrum)
    error("to be done")
    (field_r, field_t) = reserve_if_backward(D, (field_b, field_f))
    (n_r, n_t) = reverse_if_backward(D, comp.media)
    msg = ""
    (field_r.mesh == field_i.mesh) || (msg *= "Meshes of the reflected field must be the same as the input field.\n")
    (field_r.medium == field_i.medium) || (msg *= "Medium of the reflected field must be the same as the input field.\n")
    (n_r == field_r.medium) || (msg *= "Medium of the reflected field must be the same as the lens medium which the field is incident upon.\n")
    (n_t == field_t.medium) || (msg *= "Medium of the transmitted field must be the same as the lens medium which the field is transmitted to.\n")
    f = comp.focal_length
    field_t.mesh == Scale(f / n_t, f / n_t, 1)(field_i.mesh) || (msg *= "Meshes of the transmitted field must be scaled by the focal length and refractive index.\n")

    (ref_r, ref_t) = reverse_if_backward(D, comp.frames)
    (field_r.frame == ref_r) || (msg *= "Reference frame of the reflected field must be the same as the lens reference frame.\n")
    (field_t.frame == ref_t) || (msg *= "Reference frame of the transmitted field must be the same as the lens reference frame.\n")
    return isempty(msg) ? (true, msg) : (false, msg)
end

 function check_input_field(lens::Lens, field_i::MeshedBeam{T,D}) where {D,T}
    ref = (D == Forward ? first : last)(lens.frames)
    code = 0 
    (field_i.frame ≈ ref) || (code |= 1 << INVALID_FRAME)
    (field_i.medium ≈ (D == Forward ? first : last)(lens.mat)) || (code |= 1 << INVALID_MEDIUM)
    code
end

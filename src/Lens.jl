struct Lens{T, F, N, M1<:Medium{T}, M2<:Medium{T}} <: AbstractOpticalElement{T}
    focal_length::F
    numerical_aperture::N
    mat::Tuple{M1, M2}
    frames::Tuple{ReferenceFrame{T}, ReferenceFrame{T}}
    function Lens(::Type{T}, focal_length::F, numerical_aperture::N, media::Tuple{M1,M2}, frame) where {T,F,N,M1,M2}

        any(is_complex_medium, media) && throw(ArgumentError("A lens cannot be in a medium with absortion. Use media with real refractive index."))

        frame_1 = frame.origin + RotXYZ(frame.direction) * (Vec3(0, 0, -focal_length))
        frame_2 = frame.origin + RotXYZ(frame.direction) * (Vec3(0, 0, focal_length))
        new{T,F,N,M1,M2}(focal_length, numerical_aperture, media, map(o -> ReferenceFrame(o, frame.direction), (frame_1, frame_2)))
    end
end
export Lens

Lens(focal_length, numerical_aperture, media, frame) = Lens(Float64, focal_length, numerical_aperture, media, frame)

function _light_interaction!(field_b::MeshedBeam{<:Any,Backward}, field_f::MeshedBeam{<:Any,Forward}, lens::Lens, field_i::MeshedBeam{T,D,C}) where {T,D,C}
    (field_r, field_t) = reverse_if_backward(D, (field_b, field_f))
    field_r.e .= 0
    function t_aux(ind)
        t(lens, C(centroid(field_i.mesh, ind))) * field_i.e[ind]
    end
    vec(field_t.e) .= t_aux.(eachindex_nonzeros(field_i.e))
    (field_b, field_f)
end

function t(lens::Lens, coord::NSR_NSθ_λ{T}) where T<:AbstractFloat
    (nsr, θ, λ) = coord.coords
    cosθ² = 1 - (nsr / real(lens.mat[1].n))^2
    return if cosθ² < 0 # Evasnecent waves
    	zero(Complex{T})
    else # Plane waves
    	if (1 - cosθ² > lens.numerical_aperture^2) # above the lens NA
    		zero(Complex{T})
    	else
            cte = 1 / √(16π^4 / λ^2 * real(lens.mat[1].n))
    		cte / (T(lens.focal_length * cosθ²^(1/4) / (2π / λ * 2π)) / im) # might be wrong
        end
    end
end,
function t(lens::Lens, coord::NSX_NSY_λ)
    t(lens, NSR_NSθ_λ(coord))
end,
function t(lens::Lens, coord::R_θ_λ{T}) where T<:AbstractFloat
    (r, θ, λ) = coord.coords
    cosθ² = 1 - (r / lens.focal_length)^2
	if cosθ² < 0 # Evasnecent waves
		return zero(Complex{T})
	else # Plane waves
		if (1 - cosθ² > lens.numerical_aperture^2) # above the lens NA
			return zero(Complex{T})
		else
            cte = √(16π^4 / λ^2 * real(lens.mat[1].n))
			cte * (T(lens.focal_length * cosθ²^(1/4) / (2π / λ * 2π)) / im) # Needs checking
		end
	end
end,
function t(lens::Lens, coord::X_Y_λ)
    t(lens, R_θ_λ(coord))
end

function forward_backward_field(lens::Lens, field_i::MeshedBeam{T,D,C}) where {T,D,C}
    e_t = similar(field_i.e, Complex{T})
    e_r = Zeros(T, size(field_i.e))

    f = lens.focal_length
    n = D == Forward ? first(lens.mat).n : last(lens.mat).n 
    
    mesh = if C <: AngularSpectrumCoords
        Scale(T(f / n), T(f / n), T(1))(field_i.mesh)
    else
        # Scale(T(-n / f), T(-n / f), T(1))(field_i.mesh) # Miss a negative factor here due to Meshes limitations
        Affine(SA[T(-n / f) 0 0; 0 T(-n / f) 0; 0 0 1], SA[0, 0.0, 0])(field_i.mesh) # Should be using Scale but it is not implemented in Meshes
    end

    frame_t = D == Forward ? last(lens.frames) : first(lens.frames)
    (dir_r, dir_t) = reverse_if_backward(D, (Backward, Forward))
    field_r = MeshedBeam{T, dir_r, C}(field_i.mesh, e_r, field_i.medium, field_i.frame)
    field_t = MeshedBeam{T, dir_t, get_transmitted_coord_type(Lens, C)}(mesh, e_t, lens.mat[2], frame_t)
    reverse_if_backward(D, (field_r, field_t))
end

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

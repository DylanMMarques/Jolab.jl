export DielectricStack, Mirror

struct Mirror{T, R, T2, M1<:Medium{T}, M2<:Medium{T}} <: AbstractOpticalElement{T}
    reflection_coefficient::R
    transmission_coefficient::T2
    mat::Tuple{M1, M2}
    frames::Tuple{ReferenceFrame{T}, ReferenceFrame{T}}
    function Mirror(::Type{T}, media::Tuple{M1,M2}, frame; reflectivity::R) where {T,R,M1,M2}
        r = (reflectivity)^(1/2)
        t = (1 - reflectivity)^(1/2) * (real(media[1].n) / real(media[2].n))^(1/2)
        new{T,R,R,M1,M2}(r, t, (media), (frame, frame))
    end
    function Mirror(::Type{T}, r::R, t::T2, media::Tuple{M1,M2}, frame) where {T,T2,R,M1,M2}
        new{T,R,T2,M1,M2}(r, t, (media), frame)
    end
end

Mirror(media, frame; reflectivity) = Mirror(Float64, media, frame; reflectivity = reflectivity)
Mirror(r, t, media, frame) = Mirror(Float64, r, t, media, frame)

reflection_coefficient(mirror::Mirror{<:Any, <:Number}, λ, nsx, nsy) = mirror.reflection_coefficient
transmission_coefficient(mirror::Mirror{<:Any, <:Any, <:Number}, λ, nsx, nsy) = mirror.transmission_coefficient
reflection_coefficient(mirror::Mirror, λ, nsx, nsy) = mirror.reflection_coefficient(λ, nsx, nsy)
transmission_coefficient(mirror::Mirror, λ, nsx, nsy) = mirror.transmission_coefficient(λ, nsx, nsy)

function rtss(mirror::Mirror, ::Type{Forward}, nsx, nsy, λ)
    (reflection_coefficient(mirror, λ, nsx, nsy), transmission_coefficient(mirror, λ, nsx, nsy))
end

function Base.reverse(mirror::Mirror{T}) where T 
    Mirror(T, -mirror.reflection_coefficient, mirror.transmission_coefficient, reverse(mirror.mat), reverse(mirror.frames))
end


## Dielectric Stack

struct DielectricStack{T, N<:AbstractVector, H<:AbstractVector} <: AbstractOpticalElement{T}
    mat::N
    h::H
    frames::Tuple{ReferenceFrame{T}, ReferenceFrame{T}}
    function DielectricStack{T,N,H}(n::N, h::H, frame::ReferenceFrame) where {N,H,T}
        @argcheck length(n) == length(h) + 2
        last_frame = ReferenceFrame(frame.origin + RotXYZ(frame.direction) * Point3D(0,0,sum(h)), frame.direction)
        new{T,N,H}(n, h, (frame, last_frame))
    end
end

function DielectricStack(::Type{T}, n::N, h::H, frame) where {T<:Real,N<:AbstractVector{<:Medium{T,<:Union{T,Complex{T}}}}, H<:AbstractVector{T}}
    DielectricStack{T,N,H}(n, h, frame)
end

function DielectricStack(::Type{T}, n::N, h::H, frame) where {T,N<:AbstractVector{<:Medium{M1,M2}}, H<:AbstractVector{<:Real}} where {M1<:Real, M2<:Number}
    M = M2 <: Complex ? Complex{T} : T
    DielectricStack(T, convert.(Medium{T, M}, n), convert.(T, h), frame)
end
DielectricStack(n, h, frame) = DielectricStack(Float64, n, h, frame)

## Plane wave calculations

reflectioncoefficient_interfacep(n1, sz1, n2, sz2) = (n2 * sz1 - n1 * sz2) / (n2 * sz1 + n1 * sz2);

reflectioncoefficient_interfaces(n1, sz1, n2, sz2) = (n1 * sz1 - n2 * sz2) / (n1 * sz1 + n2 * sz2);

transmissioncoefficient_interfaces(n1, sz1, n2, sz2) = (2 * n1 * sz1) / (n1 * sz1 + n2 * sz2);

transmissioncoefficient_interfacep(n1, sz1, n2, sz2) = (2 * n1 * sz1) / (n2 * sz1 + n1 * sz2);

function rtss(stack::DielectricStack{<:Real, N}, ::Type{Forward}, nsx::T, nsy::T, λ) where {T, N<:AbstractVector{<:DefinedMedium}}
    nsr = √(nsx^2 + nsy^2)
    sizeA = length(stack.mat)
	sz2 = complex(1 - (nsr / stack.mat[sizeA].n)^2)^T(1/2) # Cannot use sqrt for now as Enzyme not working with sqrt of complex numbers
	sz1 = complex(1 - (nsr / stack.mat[sizeA-1].n)^2)^T(1/2)
	ri = reflectioncoefficient_interfaces(stack.mat[sizeA-1].n, sz1, stack.mat[sizeA].n, sz2)
	ti = transmissioncoefficient_interfaces(stack.mat[sizeA-1].n, sz1, stack.mat[sizeA].n, sz2)
 	imk = im * T(2π) / λ
	@inbounds for iA in (sizeA-2):-1:1
		sz2 = sz1
		sz1 = complex(1 - (nsr / stack.mat[iA].n)^2)^T(1/2)
		propagationTerm = exp(imk * stack.mat[iA+1].n * sz2 * stack.h[iA])
		rinterface = reflectioncoefficient_interfaces(stack.mat[iA].n, sz1, stack.mat[iA+1].n, sz2)
		tinterface = transmissioncoefficient_interfaces(stack.mat[iA].n, sz1, stack.mat[iA+1].n, sz2)
		ti = tinterface * ti * propagationTerm / (1 + rinterface * ri * propagationTerm^2)
		ri = (rinterface + ri * propagationTerm^2) / (1 + rinterface * ri * propagationTerm^2)
	end
	return (ri, ti)
end
@inline function rtss(stack::Union{DielectricStack, Mirror}, ::Type{Backward}, nsx, nsy, λ)
    rtss(reverse(stack), Forward, nsx, nsy, λ)
end

function Base.reverse(stack::DielectricStack{T}) where T
    mat_rev = reverse(stack.mat)
    h_rev = reverse(stack.h)
    DielectricStack(T, mat_rev, h_rev, first(stack.frames))
end

function _light_interaction(stack::Union{DielectricStack, Mirror}, pw::PlaneWaveScalar{T,D}) where {T,D}
    # change coordinate frame
    (r,t) = rtss(stack, D, pw.nsx, pw.nsy, pw.wavelength)
    r_pw = PlaneWaveScalar(T, D, pw.nsx, pw.nsy, r * pw.e, pw.wavelength, first(stack.mat), pw.frame, pw.dA)
    t_pw = PlaneWaveScalar(T, D, pw.nsx, pw.nsy, t * pw.e, pw.wavelength, last(stack.mat), pw.frame, pw.dA)
    
    return reverse_if_backward(D, (r_pw, t_pw))
end


## Beam calculations
function _light_interaction!(field_b, field_f, comp::Union{DielectricStack, Mirror}, beam::ScalarAngularSpectrumBeam{T,D}) where {T,D}
    ipw = beam.modes
    aux_struct = StructArray{Point2D{Complex{T}}}(reverse_if_backward(D, (field_b.modes.e, field_f.modes.e)))
    function f(nsx, nsy, ipwe, wavelength) 
        (a, b) = rtss(comp, D, nsx, nsy, wavelength)
        Point2D(a * ipwe, b * ipwe)
    end
    broadcast!(f, aux_struct, ipw.nsx, ipw.nsy, ipw.e, ipw.wavelength)
    (field_b, field_f)
end

function _ScatteringMatrix(field_b::Beam, field_f::Beam, comp::Union{DielectricStack{<:Any, <:AbstractVector{M2}}, Mirror{<:Any,<:Any}}, field_i::AngularSpectrumBeam{T,D}) where {T,D, M2}
    r = similar(field_i.modes.e, Complex{T})
    t = similar(r)

    aux_struct = StructArray{Point2D{Complex{T}}}((r, t))
    function f(nsx, nsy, wavelength) 
        (a, b) = rtss(comp, D, nsx, nsy, wavelength)
        Point2D(a, b)
    end
    broadcast!(f, aux_struct, field_i.modes.nsx, field_i.modes.nsy, field_i.modes.wavelength)
    
    (mat_i_to_b, mat_i_to_f) = reverse_if_backward(D, (r, t))
    @show size(mat_i_to_b)
    ScatteringMatrix(T, field_b, field_f, Diagonal(vec(mat_i_to_b)), Diagonal(vec(mat_i_to_f)), field_i)
end

function check_output_fields(field_b, field_f, comp, field_i::AngularSpectrumBeam)
    all(field_b.nsx .≈ field_f.nsx .≈ field_i.nsx) || return false
    all(field_b.nsy .≈ field_f.nsy .≈ field_i.nsy) || return false
    all(field_b.wavelength .≈ field_f.wavelength .≈ field_i.wavelength) || return false
    all(field_b.medium .≈ first.(comp.mat)) || return false
    all(field_f.medium .≈ last.(comp.mat)) || return false
    all(field_f.dA .≈ field_i.dA .≈ field_b.dA) || return false
    all(field_b.frame .≈ first(comp.frames)) || return false
    all(field_f.frame .≈ last(comp.frames)) || return false
    return true
end

function check_input_field(comp::Union{DielectricStack, Mirror}, field_i::AngularSpectrumBeam{T,D}) where {D,T}
    all(field_i.modes.frame .≈ Ref((D == Forward ? first : last)(comp.frames))) || return false
    all(field_i.modes.medium .≈ Ref((D == Forward ? first : last)(comp.mat))) || return false
    return true
end

function check_input_field(comp::Union{DielectricStack, Mirror}, field_i::PlaneWaveScalar{T,D}) where {T,D}
    field_i.frame ≈ ((D == Forward ? first : last)(comp.frames)) || return false
    field_i.medium ≈ ((D == Forward ? first : last)(comp.mat)) || return false
    return true
end

function forward_backward_field(comp::Union{DielectricStack{<:Any, <:AbstractVector{M2}}, Mirror{<:Any,<:Any,<:Any,M2}}, field_i::AngularSpectrumBeam{T,D}) where {T,D, M2}
    ipw = field_i.modes
    frame_b = Fill(first(comp.frames), size(ipw))
    frame_f = Fill(last(comp.frames), size(ipw))
    medium_b = Fill(first(comp.mat), size(ipw))
    medium_f = Fill(last(comp.mat), size(ipw))
    e_b = similar(ipw.e, Complex{T})
    e_f = similar(ipw.e, Complex{T})

    field_b = Beam(StructArray{PlaneWaveScalar{T,Backward,Complex{T},M2}}((ipw.nsx, ipw.nsy, e_b, ipw.wavelength, medium_b, frame_b, ipw.dA)))
    field_f = Beam(StructArray{PlaneWaveScalar{T,Forward,Complex{T},M2}}((ipw.nsx, ipw.nsy, e_f, ipw.wavelength, medium_f, frame_f, ipw.dA)))
    (field_b, field_f)
end

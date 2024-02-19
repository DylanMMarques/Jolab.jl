export DielectricStack, Mirror

struct Mirror{T, R, T2, M1<:Medium{T}, M2<:Medium{T}} <: AbstractOpticalElement{T}
    reflection_coefficient::R
    transmission_coefficient::T2
    mat::Tuple{M1, M2}
    frame::ReferenceFrame{T}
    function Mirror(::Type{T}, media::Tuple{M1,M2}, frame; reflectivity::R) where {T,R,M1,M2}
        r = (reflectivity)^(1/2)
        t = (1 - reflectivity)^(1/2) * (real(media[1].n) / real(media[2].n))^(1/2)
        new{T,R,R,M1,M2}(r, t, (media), frame)
    end
    function Mirror(::Type{T}, r::R, t::T2, media::Tuple{M1,M2}, frame) where {T,T2,R,M1,M2}
        new{T,R,T2,M1,M2}(r, t, (media), frame)
    end
end

Mirror(reflectivity, media, frame) = Mirror(Float64, reflectivity, media, frame)

reflection_coefficient(mirror::Mirror{<:Any, <:Number}, λ, nsx, nsy) = mirror.reflection_coefficient
transmission_coefficient(mirror::Mirror{<:Any, <:Any, <:Number}, λ, nsx, nsy) = mirror.transmission_coefficient
reflection_coefficient(mirror::Mirror, λ, nsx, nsy) = mirror.reflection_coefficient(λ, nsx, nsy)
transmission_coefficient(mirror::Mirror, λ, nsx, nsy) = mirror.transmission_coefficient(λ, nsx, nsy)

function rtss(mirror::Mirror, ::Type{Forward}, nsx, nsy, λ)
    (reflection_coefficient(mirror, λ, nsx, nsy), transmission_coefficient(mirror, λ, nsx, nsy))
end

function Base.reverse(mirror::Mirror{T}) where T 
    Mirror(T, -mirror.reflection_coefficient, mirror.transmission_coefficient, reverse(mirror.mat), mirror.frame)
end


## Dielectric Stack

struct DielectricStack{T, N<:AbstractVector, H<:AbstractVector} <: AbstractOpticalElement{T}
    mat::N
    h::H
    frame::ReferenceFrame{T}
    function DielectricStack{T,N,H}(n::N, h::H, frame::ReferenceFrame) where {N,H,T}
        @argcheck length(n) == length(h) + 2
        new{T,N,H}(n, h, frame)
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
 	imk= im * T(2π) / λ
	@inbounds for iA in (sizeA-2):-1:1
		sz2 = sz1
		sz1 = complex(1 - (nsr / stack.mat[iA].n)^2)^T(1/2)
		propagationTerm = exp(imsqrtk * stack.mat[iA+1].n * sz2 * stack.h[iA])
		rinterface = reflectioncoefficient_interfaces(stack.mat[iA].n, sz1, stack.mat[iA+1].n, sz2)
		tinterface = transmissioncoefficient_interfaces(stack.mat[iA].n, sz1, stack.mat[iA+1].n, sz2)
		ti = tinterface * ti * propagationTerm / (1 + rinterface * ri * propagationTerm^2)
		ri = (rinterface + ri * propagationTerm^2) / (1 + rinterface * ri * propagationTerm^2)
	end
	return (ri, ti)
end
@inline function rtss(stack::Union{DielectricStack, Mirror}, dir::Type{Backward}, nsx, nsy, λ)
    rtss(reverse(stack), Forward, nsx, nsy, λ)
end

function Base.reverse(stack::DielectricStack{T}) where T
    mat_rev = reverse(stack.mat)
    h_rev = reverse(stack.h)
    total_h = sum(stack.h)
    origin_rev = stack.frame.origin + RotXYZ(stack.frame.direction) * Point3D(0,0,total_h)
    frame_rev = ReferenceFrame(origin_rev, stack.frame.direction)
    DielectricStack(T, mat_rev, h_rev, frame_rev)
end


function notchecked_light_interaction(stack::Union{DielectricStack, Mirror}, pw::PlaneWaveScalar{T,D}) where {T,D}
    # change coordinate frame
    (r,t) = rtss(stack, D, pw.nsx, pw.nsy, pw.wavelength)
    r_pw = PlaneWaveScalar(T, D, pw.nsx, pw.nsy, r * pw.e, pw.wavelength, first(stack.mat), pw.frame, pw.dA)
    t_pw = PlaneWaveScalar(T, D, pw.nsx, pw.nsy, t * pw.e, pw.wavelength, last(stack.mat), pw.frame, pw.dA)
    
    return reverse_if_backward(D, (r_pw, t_pw))
end

function light_interaction(stack::Union{DielectricStack, Mirror}, pw::PlaneWaveScalar)
    # check if medium is wavelength defined
    @argcheck first(stack.mat) ≈ pw.medium ArgumentError
    @argcheck stack.frame ≈ pw.frame ArgumentError
    notchecked_light_interaction(stack, pw)
end

## Beam calculations
function light_interaction(comp::Union{DielectricStack{<:Any, <:AbstractVector{M2}}, Mirror{<:Any, <:Any, <:Any, <:Any, M2}}, beam::ScalarAngularSpectrumBeam{T,D,E,M}) where {M2,D,T,E,M}
    ipw = beam.modes
    @argcheck all(ipw.frame .≈ Ref(comp.frame)) ArgumentError
    @argcheck all(ipw.medium .≈ Ref(first(comp.mat))) ArgumentError
    
    r = similar(ipw.e, Complex{T})
    t = similar(ipw.e, Complex{T})
    aux_struct = StructArray{Point2D{Complex{T}}}((r, t))
    
    function f(nsx, nsy, ipwe, wavelength) 
        (a, b) = rtss(comp, D, nsx, nsy, wavelength)
        Point2D(a * ipwe, b * ipwe)
    end
    broadcast!(f, aux_struct, ipw.nsx, ipw.nsy, ipw.e, ipw.wavelength)

    mat_t = Fill((D == Forward ?  last : first)(comp.mat), size(beam)) 
    r_beam = Beam(StructArray{PlaneWaveScalar{T,!D,Complex{T},M}}((ipw.nsx, ipw.nsy, r, ipw.wavelength, ipw.medium, ipw.frame, ipw.dA)))
    t_beam = Beam(StructArray{PlaneWaveScalar{T,D,Complex{T},M2}}((ipw.nsx, ipw.nsy, t, ipw.wavelength, mat_t, ipw.frame, ipw.dA)))

    return reverse_if_backward(D, (r_beam, t_beam))
end


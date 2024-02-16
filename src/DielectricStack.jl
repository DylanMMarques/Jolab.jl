export DielectricStack

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

function rtss(stack::DielectricStack{<:Real, N}, nsx::T, nsy::T, λ) where {T, N<:AbstractVector{<:DefinedMedium}}
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
		propagationTerm = exp(imk * stack.mat[iA+1].n * sz2 * stack.h[iA])
		rinterface = reflectioncoefficient_interfaces(stack.mat[iA].n, sz1, stack.mat[iA+1].n, sz2)
		tinterface = transmissioncoefficient_interfaces(stack.mat[iA].n, sz1, stack.mat[iA+1].n, sz2)
		ti = tinterface * ti * propagationTerm / (1 + rinterface * ri * propagationTerm^2)
		ri = (rinterface + ri * propagationTerm^2) / (1 + rinterface * ri * propagationTerm^2)
	end
	return (ri, ti)
end

function notchecked_light_interaction(stack::DielectricStack, pw::PlaneWaveScalar{T}) where T
    # change coordinate frame
    (r,t) = rtss(stack, pw.nsx, pw.nsy, pw.wavelength)
    r_pw = PlaneWaveScalar(T, pw.nsx, pw.nsy, r * pw.e, pw.wavelength, first(stack.mat), pw.frame, pw.dA)
    t_pw = PlaneWaveScalar(T, pw.nsx, pw.nsy, t * pw.e, pw.wavelength, last(stack.mat), pw.frame, pw.dA)
    return (r_pw, t_pw)
end


function light_interaction(stack::DielectricStack, pw::PlaneWaveScalar)
    # check if medium is wavelength defined
    @argcheck first(stack.mat) ≈ pw.medium ArgumentError
    @argcheck stack.frame ≈ pw.frame ArgumentError
    notchecked_light_interaction(stack, pw)
end

## Beam calculations


const ScalarAngularSpectrumBeam{T,E,M,D,V,C} = Beam{T, StructArray{PlaneWaveScalar{T,E,M},D,V,C}} 
function light_interaction(comp::DielectricStack{<:Any, <:AbstractVector{M2}}, beam::ScalarAngularSpectrumBeam{T,E,M}) where {M2,T,E,M}
    ipw = beam.modes
    @argcheck all(ipw.frame .≈ Ref(comp.frame)) ArgumentError
    @argcheck all(ipw.medium .≈ Ref(first(comp.mat))) ArgumentError
    
    r = similar(ipw.e, Complex{T})
    t = similar(ipw.e, Complex{T})
    aux_struct = StructArray{Point2D{Complex{T}}}((r, t))
    
    function f(nsx, nsy, ipwe, wavelength) 
        (a, b) = rtss(comp, nsx, nsy, wavelength)
        Point2D(a * ipwe, b * ipwe)
    end
    broadcast!(f, aux_struct, ipw.nsx, ipw.nsy, ipw.e, ipw.wavelength)
    r_modes = StructArray{PlaneWaveScalar{T,Complex{T},M}}((ipw.nsx, ipw.nsy, r, ipw.wavelength, ipw.medium, ipw.frame, ipw.dA))
    t_modes = StructArray{PlaneWaveScalar{T,Complex{T},M2}}((ipw.nsx, ipw.nsy, t, ipw.wavelength, Fill(last(comp.mat), size(ipw.nsx)), ipw.frame, ipw.dA))

    return (Beam(r_modes), Beam(t_modes))
end

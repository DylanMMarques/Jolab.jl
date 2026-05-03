export DielectricStack, Mirror

struct Mirror{T,R,T2,M1<:Medium{T},M2<:Medium{T}} <: AbstractOpticalElement{T}
    reflection_coefficient::R
    transmission_coefficient::T2
    mat::Tuple{M1,M2}
    frames::Tuple{ReferenceFrame{T},ReferenceFrame{T}}
    function Mirror(
        ::Type{T},
        media::Tuple{M1,M2},
        frame;
        reflectivity::R,
    ) where {T,R,M1,M2}
        r = √reflectivity
        t = √(1 - reflectivity)
        new{T,R,R,M1,M2}(r, t, (media), (frame, frame))
    end
    function Mirror(::Type{T}, r::R, t::T2, media::Tuple{M1,M2}, frame) where {T,T2,R,M1,M2}
        new{T,R,T2,M1,M2}(r, t, (media), frame)
    end
end

Mirror(media, frame; reflectivity) =
    Mirror(Float64, media, frame; reflectivity = reflectivity)
Mirror(r, t, media, frame) = Mirror(Float64, r, t, media, frame)

reflection_coefficient(mirror::Mirror{<:Any,<:Number}) = mirror.reflection_coefficient
transmission_coefficient(mirror::Mirror{<:Any,<:Any,<:Number}) =
    mirror.transmission_coefficient

function rtss(mirror::Mirror, ::Type{Forward}, coords::NSR_NSθ_λ)
    r = reflection_coefficient(mirror)
    t = transmission_coefficient(mirror)

    # nsz1 = √(first(mirror.mat).n^2 - coords.coords[1]^2)
    # nsz2 = √(last(mirror.mat).n^2 - coords.coords[1]^2)
    return (r, t)
end,
function rtss(mirror::Mirror, ::Type{Forward}, coords::NSX_NSY_λ)
    rtss(mirror, Forward, NSR_NSθ_λ(coords))
end

function Base.reverse(mirror::Mirror{T}) where {T}
    Mirror(
        T,
        -mirror.reflection_coefficient,
        mirror.transmission_coefficient,
        reverse(mirror.mat),
        reverse(mirror.frames),
    )
end


## Dielectric Stack

struct DielectricStack{T,N<:AbstractVector,H<:AbstractVector} <: AbstractOpticalElement{T}
    mat::N
    h::H
    frames::Tuple{ReferenceFrame{T},ReferenceFrame{T}}
    function DielectricStack{T,N,H}(n::N, h::H, frame::ReferenceFrame) where {N,H,T}
        @argcheck length(n) == length(h) + 2
        last_frame = ReferenceFrame(
            frame.origin + X_Y_Z(_RotXYZ(frame.direction) * SVector((T(0), T(0), sum(h)))),
            frame.direction,
        )
        new{T,N,H}(n, h, (frame, last_frame))
    end
end


function DielectricStack(
    ::Type{T},
    n::N,
    h::H,
    frame,
) where {T<:Real,N<:AbstractVector{<:Medium{T,<:Union{T,Complex{T}}}},H<:AbstractVector{T}}
    DielectricStack{T,N,H}(n, h, frame)
end

function DielectricStack(
    ::Type{T},
    n::N,
    h::H,
    frame,
) where {
    T,
    N<:AbstractVector{<:Medium{M1,M2}},
    H<:AbstractVector{<:Real},
} where {M1<:Real,M2<:Number}
    M = M2 <: Complex ? Complex{T} : T
    DielectricStack(T, convert.(Medium{T,M}, n), convert.(T, h), frame)
end
DielectricStack(n, h, frame) = DielectricStack(Float64, n, h, frame)

## Plane wave calculations

reflectioncoefficient_interfacep(n1, sz1, n2, sz2) =
    (n2 * sz1 - n1 * sz2) / (n2 * sz1 + n1 * sz2);

reflectioncoefficient_interfaces(n1, sz1, n2, sz2) =
    (n1 * sz1 - n2 * sz2) / (n1 * sz1 + n2 * sz2);

transmissioncoefficient_interfaces(n1, sz1, n2, sz2) =
    (2 * n1 * sz1) / (n1 * sz1 + n2 * sz2);

transmissioncoefficient_interfacep(n1, sz1, n2, sz2) =
    (2 * n1 * sz1) / (n2 * sz1 + n1 * sz2);

function recursive_stack(
    f_r,
    f_t,
    stack::DielectricStack{<:Real,N},
    nsr::T,
    λ,
) where {T,N<:AbstractVector{<:DefinedMedium}}
    sizeA = length(stack.mat)
    sz2 = √(complex(1 - (nsr / stack.mat[sizeA].n)^2))
    sz1 = √(complex(1 - (nsr / stack.mat[sizeA-1].n)^2))
    ri = f_r(stack.mat[sizeA-1].n, sz1, stack.mat[sizeA].n, sz2)
    ti = f_t(stack.mat[sizeA-1].n, sz1, stack.mat[sizeA].n, sz2)
    imk = im * T(2π) / λ
    @inbounds for iA = (sizeA-2):-1:1
        sz2 = sz1
        sz1 = √(complex(1 - (nsr / stack.mat[iA].n)^2))
        propagationTerm = exp(imk * stack.mat[iA+1].n * sz2 * stack.h[iA])
        rinterface = f_r(stack.mat[iA].n, sz1, stack.mat[iA+1].n, sz2)
        tinterface = f_t(stack.mat[iA].n, sz1, stack.mat[iA+1].n, sz2)
        ti = tinterface * ti * propagationTerm / (1 + rinterface * ri * propagationTerm^2)
        ri =
            (rinterface + ri * propagationTerm^2) /
            (1 + rinterface * ri * propagationTerm^2)
    end
    return (ri, ti)
end

for i in (:s, :p)
    f = Symbol(:rt, i, i)
    @eval begin
        function $f(
            stack::DielectricStack{<:Real,N},
            ::Type{Forward},
            nsr::T,
            λ,
        ) where {T,N<:AbstractVector{<:DefinedMedium}}
            recursive_stack(
                $(Symbol(:reflectioncoefficient_interface, i)),
                $(Symbol(:transmissioncoefficient_interface, i)),
                stack,
                nsr,
                λ,
            )
        end

        $f(stack::DielectricStack, ::Type{Forward}, coords::NSX_NSY_λ) =
            $f(stack, Forward, √(coords[1]^2 + coords[2]^2), coords[3])
        $f(stack::DielectricStack, ::Type{Forward}, coords::NSR_NSθ_λ) =
            $f(stack, Forward, coords[1], coords[3])
    end
end

@inline function rtss(stack::Union{DielectricStack,Mirror}, ::Type{Backward}, coords)
    rtss(reverse(stack), Forward, coords)
end

function Base.reverse(stack::DielectricStack{T}) where {T}
    mat_rev = reverse(stack.mat)
    h_rev = reverse(stack.h)
    DielectricStack(T, mat_rev, h_rev, first(stack.frames))
end

## Beam calculations
function _light_interaction!(
    field_b,
    field_f,
    comp::Union{DielectricStack,Mirror},
    beam::MeshedAngularSpectrum{T,D,C,P},
) where {T,D,C<:AngularSpectrumCoords,P}
    (field_r, field_t) = reverse_if_backward(D, (field_b.e, field_f.e))

    polarization_components = number_components(P)
    for iP = 1:polarization_components
        field_i_e = view(beam.e,:,:,:,iP)
        field_r_e = view(field_r,:,:,:,iP)
        field_t_e = view(field_t,:,:,:,iP)

        f = if isone(iP)
            (index, e) -> begin
                (r, t) = rtss(comp, D, C(centroid(beam.mesh, index)))
                (r * e, t * e)
            end
        else
            (index, e) -> begin
                (r, t) = rtpp(comp, D, C(centroid(beam.mesh, index)))
                (r * e, t * e)
            end
        end
        tmp = StructArray{Tuple{Complex{T},Complex{T}}}((vec(field_r_e), vec(field_t_e)))
        tmp .= f.(eachindex(field_i_e), vec(field_i_e))
    end
    (field_b, field_f)
end

# Need inline or bug with enzyme
@inline function _ScatteringMatrix(
    field_b::MeshedAngularSpectrum,
    field_f::MeshedAngularSpectrum,
    comp::Union{DielectricStack{<:Any,<:AbstractVector{M2}},Mirror{<:Any,<:Any}},
    field_i::MeshedAngularSpectrum{T,D,C},
) where {T,D,M2,C}
    r = similar(field_i.e, Complex{T}, length(field_i.e))
    t = similar(r)

    f(index) = rtss(comp, D, C(centroid(field_i.mesh, index)))
    tmp = StructArray{Tuple{Complex{T},Complex{T}}}((vec(r), vec(t)))
    tmp .= f.(eachindex(field_i.e))

    (mat_i_to_b, mat_i_to_f) = reverse_if_backward(D, (r, t))
    ScatteringMatrix(
        T,
        field_b,
        field_f,
        Diagonal(vec(mat_i_to_b)),
        Diagonal(vec(mat_i_to_f)),
        field_i,
    )
end

function check_input_field(
    comp::Union{DielectricStack,Mirror},
    field_i::MeshedAngularSpectrum{T,D},
) where {D,T}
    msg_code = zero(UInt64)
    field_i.frame ≈ (D == Forward ? first : last)(comp.frames) ||
        (msg_code |= 1 << INVALID_FRAME)
    field_i.medium ≈ (D == Forward ? first : last)(comp.mat) ||
        (msg_code |= 1 << INVALID_MEDIUM)
    msg_code
end

function forward_backward_field(
    comp::Union{DielectricStack{<:Any,<:AbstractVector{M2}},Mirror{<:Any,<:Any,<:Any,M2}},
    field_i::MeshedAngularSpectrum{T,D,C,P},
) where {T,D,M2,C,P}
    frame_b = first(comp.frames)
    frame_f = last(comp.frames)
    medium_b = first(comp.mat)
    medium_f = last(comp.mat)
    e_b = similar(field_i.e, Complex{T})
    e_f = similar(field_i.e, Complex{T})

    field_b = MeshedBeam{T,Backward,C,P}(field_i.mesh, e_b, medium_b, frame_b)
    field_t = MeshedBeam{T,Forward,C,P}(field_i.mesh, e_f, medium_f, frame_f)
    (field_b, field_t)
end

export PlaneWaveScalar, PlaneWaveVectorial

struct PlaneWaveScalar{T, E<:Union{T, Complex{T}}, M<:Medium{T, <:Union{T, Complex{T}}}} <: AbstractPlaneWave{T}
    nsx::T
    nsy::T
    e::E
    wavelength::T
    medium::M
    frame::ReferenceFrame{T}
    dA::T
end

const FieldVector2D3D{T} = Union{FieldVector{2, T}, FieldVector{3, T}}

struct PlaneWaveVectorial{T, E<:FieldVector2D3D{<:Union{Complex{T}, T}}, M<:Medium}
    nsx::T
    nsy::T
    e::E
    wavelength::T
    medium::M
    frame::ReferenceFrame{T}
    dA::T
end

function PlaneWaveScalar(::Type{T}, nsx, nsy, e::E1, wavelength, medium::Medium{M1, M2}, frame, dA = one(T)) where {T <: AbstractFloat, E1, M1, M2<:Number}
    E = E1 <: Complex ? Complex{T} : T
    M = M2 <: Complex ? Complex{T} : T
    PlaneWaveScalar{T, E, Medium{T, M}}(nsx, nsy, e, wavelength, medium, frame, dA)
end

PlaneWaveScalar(arg...) = PlaneWaveScalar(Float64, arg...)
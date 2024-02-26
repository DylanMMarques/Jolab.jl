export PlaneWaveScalar, PlaneWaveVectorial

struct PlaneWaveScalar{T, D, E<:Union{T, Complex{T}}, M<:Medium{T, <:Union{T, Complex{T}}}} <: AbstractPlaneWave{T,D}
    nsx::T
    nsy::T
    e::E
    wavelength::T
    medium::M
    frame::ReferenceFrame{T}
    dA::T
end

const FieldVector2D3D{T} = Union{FieldVector{2, T}, FieldVector{3, T}}

struct PlaneWaveVectorial{T, D, E<:FieldVector2D3D{<:Union{Complex{T}, T}}, M<:Medium} <: AbstractPlaneWave{T,D}
    nsx::T
    nsy::T
    e::E
    wavelength::T
    medium::M
    frame::ReferenceFrame{T}
    dA::T
end

function PlaneWaveScalar(::Type{T}, ::Type{D}, nsx, nsy, e::E1, wavelength, medium::Medium{M1, M2}, frame, dA = one(T)) where {T <: AbstractFloat, D<:AbstractDirection, E1, M1, M2<:Number}
    E = E1 <: Complex ? Complex{T} : T
    M = M2 <: Complex ? Complex{T} : T
    PlaneWaveScalar{T, D, E, Medium{T, M}}(nsx, nsy, e, wavelength, medium, frame, dA)
end
PlaneWaveScalar(::Type{D}, nsx, nsy, e, wavelength, medium, frame, dA = one(Float64)) where D = PlaneWaveScalar(Float64, D, nsx, nsy, e, wavelength, medium, frame, dA)

intensity(pw::PlaneWaveScalar) = abs2(pw.e) * pw.dA * pw.medium.n
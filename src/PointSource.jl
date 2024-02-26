export ScalarPointSource, VectorialPointSource

struct ScalarPointSource{T, D, E<:RealOrComplex{T}, M<:Medium{T, <:RealOrComplex{T}}} <: AbstractPointSource{T,D}
    x::T
    y::T
    e::E
    wavelength::T
    medium::M
    frame::ReferenceFrame{T}
    dA::T
end

struct VectorialPointSource{T, D, E<:FieldVector{3, RealOrComplex{T}}, M<:Medium{T, <:RealOrComplex}} <: AbstractPointSource{T,D}
    x::T
    y::T
    e::E
    wavelength::T
    medium::M
    frame::ReferenceFrame{T}
    dA::T
end

function ScalarPointSource(::Type{T}, ::Type{D}, x, y, e::E1, wavelength, medium::Medium{M1, M2}, frame, dA = one(T)) where {T <: AbstractFloat, D<:AbstractDirection, E1, M1, M2<:Number}
    E = E1 <: Complex ? Complex{T} : T
    M = M2 <: Complex ? Complex{T} : T
    ScalarPointSource{T, D, E, Medium{T, M}}(x, y, e, wavelength, medium, frame, dA)
end
ScalarPointSource(::Type{D}, x, y, e, wavelength, medium, frame, dA = one(Float64)) where D = ScalarPointSource(Float64, D, x, y, e, wavelength, medium, frame, dA)
ScalarPointSource(::Type{T}, ::Type{D}, x, y, e, wavelength, medium, frame, dA = one(T)) where {T, D} = ScalarPointSource(T, D, x, y, e, wavelength, medium, frame, dA)

intensity(pw::ScalarPointSource) = abs2(pw.e) * pw.dA
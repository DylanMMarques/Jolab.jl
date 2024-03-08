struct Beam{T, D, M<:StructArray{<:AbstractFieldMode{T,D}}} <: AbstractField{T,D}
    modes::M
    function Beam(modes::M) where {M<:StructArray{<:AbstractFieldMode{T,D}}} where {T,D}
        new{T, D, M}(modes)
    end
end

const ScalarAngularSpectrumBeam{T,D,E,M,N,V,C} = Beam{T, D, StructArray{PlaneWaveScalar{T,D,E,M},N,V,C}} 
const VectorialAngularSpectrumBeam{T,D,E,M,N,V,C} = Beam{T, D, StructArray{PlaneWaveVectorial{T,D,E,M},N,V,C}} 
const AngularSpectrumBeam{T,D,E,M,N,V,C} = Union{ScalarAngularSpectrumBeam{T,D,E,M,N,V,C}, VectorialAngularSpectrumBeam{T,D,E,M,N,V,C}}

const ScalarSpatialBeam{T,D,E,M,N,V,C} = Beam{T, D, StructArray{ScalarPointSource{T,D,E,M},N,V,C}} 
const VectorialSpatialBeam{T,D,E,M,N,V,C} = Beam{T, D, StructArray{VectorialPointSource{T,D,E,M},N,V,C}} 
const SpatialBeam{T,D,E,M,N,V,C} = Union{ScalarSpatialBeam{T,D,E,M,N,V,C}, VectorialSpatialBeam{T,D,E,M,N,V,C}}


function monochromatic_angularspectrum(::Type{T}, ::Type{D}, nsx, nsy, e::AbstractArray{E}, wavelength, medium::M, frame::ReferenceFrame{T}) where {T, M<:Medium, E,D}
    E2 = E <: Complex ? Complex{T} : T

    number_modes = length(e)
    @argcheck all(i -> length(i) == number_modes, (nsx, nsy)) DimensionMismatch

    modes = StructArray{PlaneWaveScalar{T,D,E2,M}}((vec(nsx), vec(nsy), vec(e), Fill(wavelength, number_modes), Fill(medium, number_modes), Fill(frame, number_modes), Fill((nsx[2] - nsx[1]) * (nsy[1,2] - nsy[1]), number_modes)))
    Beam(modes)
end
monochromatic_angularspectrum(::Type{D}, nsx, nsy, e, wavelength, medium, frame) where D = monochromatic_angularspectrum(Float64, D, nsx, nsy, e, wavelength, medium, frame)

function monochromatic_fieldspace(::Type{T}, ::Type{D}, x, y, e::AbstractArray{E}, wavelength, medium::M, frame::ReferenceFrame{T}) where {T, M<:Medium, E,D}
    E2 = E <: Complex ? Complex{T} : T

    number_modes = size(e)
    @argcheck all(i -> size(i) == number_modes, (x, y)) DimensionMismatch
    modes = StructArray{ScalarPointSource{T,D,E2,M}}((x, y, e, Fill(wavelength, number_modes), Fill(medium, number_modes), Fill(frame, number_modes), Fill((x[2] - x[1]) * (y[1,2] - y[1]), number_modes)))
    Beam(modes)
end
monochromatic_fieldspace(::Type{D}, x, y, e, wavelength, medium, frame) where D = monochromatic_fieldspace(Float64, D, x, y, e, wavelength, medium, frame)

function intensity(beam::ScalarAngularSpectrumBeam) 
    f(e, dA, mat) = abs2(e) * (dA * mat.n)
    mapreduce(f, +, beam.modes.e, beam.modes.dA, beam.modes.medium)
end

function Base.isapprox(a::Beam{T1,D1,<:StructArray{M1}}, b::Beam{T2,D2,<:StructArray{M2}}; kwargs...) where {T1, T2, D1, D2, M1, M2}
    D1 == D2 || return false
    same_mode_type(M1, M2) || return false
    fields = fieldnames(M1)
    for field in fields
        all(isapprox.(component(a.modes, field), component(b.modes, field); kwargs...)) || return false
    end
    return true
end
coord_list = (:X_Y_λ,
    :X_Y_t,
    :X_NSY_λ,
    :X_NSY_t,
    :NSX_Y_λ,
    :NSX_Y_t,
    :NSX_NSY_λ,
    :NSX_NSY_t,
    :R_θ_λ,
    :R_θ_t,
    :NSR_NSθ_λ,
    :NSR_NSθ_t
)
for coordi in coord_list
    eval(quote
        struct $coordi{T}
            coords::Vec{3, T}
        end
        Base.getindex(c::$coordi, i::Integer) = c.coords[i]
        $coordi(coords::Point3) = $coordi(coords.coords)
    end)
end

struct MeshedBeam{T, D, C, V, E<:AbstractArray{<:RealOrComplex{T},3}, M<:Medium{T,<:RealOrComplex{T}}, T2<:RealOrComplex{T}} <: AbstractField{T, D}
    mesh::V
    e::E
    medium::M
    frame::ReferenceFrame{T}
    function MeshedBeam{T,D,C}(mesh::V, e::E, medium::Medium{<:Any, M}, frame::ReferenceFrame) where {T,V,E <: AbstractArray{T2},D,C,M} where T2
        @argcheck size(e) == size(mesh) DimensionMismatch
        M1 = M <: Complex ? Complex{T} : T
        T3 = T2 <: Complex ? Complex{T} : T
        new{T,D,C,V,E,Medium{T, M1},T3}(mesh, e, medium, frame)
    end
end

function MonochromaticAngularSpectrum(::Type{T}, ::Type{D}, nsx::AbstractRange, nsy::AbstractRange, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {T, D}
    mesh = CartesianGrid((length(nsx), length(nsy), 1),
        Point3(first(nsx), first(nsy), λ - eps(T) / 2), 
        (step(nsx), step(nsy), eps(T)))
    MeshedBeam{T, D, NSX_NSY_λ}(mesh, reshape(e, size(e)..., 1), medium, frame)
end
MonochromaticAngularSpectrum(::Type{D}, nsx::AbstractRange, nsy::AbstractRange, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {D} = MonochromaticAngularSpectrum(Float64, D, nsx, nsy, e, λ, medium, frame)

function MonochromaticSpatialBeam(::Type{T}, ::Type{D}, x::AbstractVector, y::AbstractVector, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {T, D}
    @argcheck size(e) == (length(x), length(y)) DimensionMismatch
    mesh = CartesianGrid((length(x), length(y), 1),
        Point3(first(x), first(y), λ - eps(T) / 2),  # The -0.5 is to center the point on the face
        (step(x), step(y), eps(T)))
    MeshedBeam{T, D, X_Y_λ}(mesh, reshape(e, size(e)..., 1), medium, frame)
end
MonochromaticSpatialBeam(::Type{D}, x::AbstractRange, y::AbstractRange, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {D} = MonochromaticSpatialBeam(Float64, D, x, y, e, λ, medium, frame)

function intensity(beam::MeshedBeam)
    f(e, ind) = abs2(e) * (volume(beam.mesh, ind) * beam.medium.n)
    mapreduce(f, +, vec(beam.e), eachindex(beam.e))
end

const AngularSpectrumCoords = Union{NSX_NSY_λ, NSR_NSθ_λ}
const MeshedAngularSpectrum{T,D} = MeshedBeam{T,D,<:AngularSpectrumCoords}

function MeshedPlaneWaveScalar(::Type{T}, ::Type{D}, nsx, nsy, e::T2, λ, medium, frame) where {T,D,T2}
    mesh = CartesianGrid((1, 1, 1),
        Point{3,T}(nsx - eps(T) / 2, nsy - eps(T) / 2, λ - eps(T) / 2), 
        (eps(T), eps(T), eps(T))
        )
    TE = T2 <: Complex ? Complex{T} : T
    MeshedBeam{T, D, NSX_NSY_λ}(mesh, (@SArray [TE(e);;;]), medium, frame)
end
MeshedPlaneWaveScalar(::Type{D}, nsx, nsy, e, λ, medium, frame) where D = MeshedPlaneWaveScalar(Float64, D, nsx, nsy, e, λ, medium, frame)

function Base.isapprox(beam1::MeshedBeam{T1,D,C}, beam2::MeshedBeam{T2,D,C}; kwargs...) where {T1, T2, D, C}
    isapprox(beam1.mesh, beam2.mesh; kwargs...) || return false
    isapprox(beam1.frame, beam2.frame; kwargs...) || return false
    isapprox(beam1.medium, beam2.medium; kwargs...) || return false
    return true
end

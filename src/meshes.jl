abstract type Domain{C, Dim, T} end
abstract type AbstractPolarization end

struct PolarizationSP <: AbstractPolarization end 
struct PolarizationXYZ <: AbstractPolarization end
struct PolarizationScalar <: AbstractPolarization end

number_components(::Type{PolarizationSP}) = 2
number_components(::Type{PolarizationXYZ}) = 3
number_components(::Type{PolarizationScalar}) = 1

struct X_Y_Z{T} <: FieldVector{3, T}
    x::T
    y::T
    z::T
end

struct X_Y_λ{T} <: FieldVector{3, T}
    x::T
    y::T
    λ::T
end

struct R_θ_Z{T} <: FieldVector{3, T}
    r::T
    θ::T
    z::T
end

struct R_θ_ϕ{T} <: FieldVector{3, T}
    r::T
    θ::T
    ϕ::T
end

struct R_θ_λ{T} <: FieldVector{3, T}
    r::T
    θ::T
    λ::T
end

struct NSX_NSY_λ{T} <: FieldVector{3, T}
    nsx::T
    nsy::T
    λ::T
end

struct NSR_NSθ_λ{T} <: FieldVector{3, T}
    nsr::T
    nsθ::T
    λ::T
end

function unit_vector(coord::X_Y_Z{T}, ::Val{:x}) where T
    X_Y_Z(one(T), zero(T), zero(T))
end,
function unit_vector(coord::X_Y_Z{T}, ::Val{:y}) where T
    X_Y_Z(zero(T), one(T), zero(T))
end,
function unit_vector(coord::X_Y_Z{T}, ::Val{:Z}) where T
    X_Y_Z(zero(T), zero(T), one(T))
end,
function unit_vector(coord::R_θ_Z{T}, ::Val{:R}) where T
    X_Y_Z(cos(coord.θ), sin(coord.θ), zero(T))
end,
function unit_vector(coord::R_θ_Z{T}, ::Val{:θ}) where T
    X_Y_Z(-sin(coord.θ), cos(coord.θ), zero(T))
end,
function unit_vector(coord::R_θ_Z{T}, ::Val{:Z}) where T
    X_Y_Z(zero(T), zero(T), one(T))
end,
function unit_vector(coord::R_θ_ϕ{T}, ::Val{:R}) where T
    X_Y_Z(sin(coord.θ)*cos(coord.ϕ), sin(coord.θ)*sin(coord.ϕ), cos(coord.θ))
end,
function unit_vector(coord::R_θ_ϕ{T}, ::Val{:θ}) where T
    X_Y_Z(cos(coord.θ)*cos(coord.ϕ), cos(coord.θ)*sin(coord.ϕ), -sin(coord.θ))
end,
function unit_vector(coord::R_θ_ϕ{T}, ::Val{:ϕ}) where T
    X_Y_Z(-sin(coord.ϕ), cos(coord.ϕ), zero(T))
end

for (cart, polar) in zip((:NSX_NSY_λ, :X_Y_λ, :X_Y_t), (:NSR_NSθ_λ, :R_θ_λ))
    eval(quote
        function $cart(coords::$polar{T}) where T
            car = CartesianFromPolar()(Polar(coords[1], coords[2]))
            $cart(car.x, car.y, coords[3])
        end
        function $polar(coords::$cart{T}) where T
            car = PolarFromCartesian()(@view coords[1:2])
            $polar(car.r, car.θ, coords[3])
        end
    end)
end

struct CartesianGrid{C, Dim, T} <: Domain{C, Dim,T}
    spacing::NTuple{Dim, T}
    origin::C
    lengths::NTuple{Dim, Int}
    function CartesianGrid(lengths::NTuple{Dim, Int}, origin::C, spacing::NTuple{Dim, T}) where {Dim, T, C}
        new{C, Dim, T}(spacing, origin, lengths)
    end
end

function CartesianGrid(::Type{C}, x::AbstractRange, y::AbstractRange, z::AbstractRange) where C
    xyz = (x, y, z)
    lengths = length.(xyz)
    spacing = step.(xyz)
    start = first.(xyz)
    CartesianGrid(lengths, C(start...), spacing)
end


area(grid::CartesianGrid{C, 2}, ind::Union{Integer, CartesianIndex{2}}) where C = prod(grid.spacing)
volume(grid::CartesianGrid{C, 3}, ind::Union{Integer, CartesianIndex{3}}) where C = prod(grid.spacing)
function centroid(grid::CartesianGrid{C, Dim}, ind::CartesianIndex{Dim}) where {C, Dim}
    @boundscheck checkbounds(grid, ind)
    C(grid.origin .+ grid.spacing .* (ind.I .- 1))
end

function Base.checkbounds(grid::Domain{C, Dim}, ind::CartesianIndex{Dim}) where {C, Dim}
    all(i -> 1 <= i[1] <= i[2], Iterators.zip(ind.I, grid.lengths)) || throw(BoundsError(grid, ind))
end,
function Base.checkbounds(grid::Domain{C, Dim}, ind::Integer) where {C, Dim}
    1 <= ind <= nelements(grid) || throw(BoundsError(grid, ind))
end

struct CylindricalGrid{C, Dim,T} <: Domain{C, Dim, T}
    spacing::NTuple{Dim, T}
    origin::C
    lengths::NTuple{Dim, Int}
    function CylindricalGrid(lengths::NTuple{Dim, Int}, origin::C, spacing::NTuple{Dim, T}) where {C, Dim, T}
        θ_start = origin[2]
        θ_stop = θ_start + spacing[2] * (lengths[2] - 1)
        0 <= θ_start <= 2π || throw(ArgumentError("CylindricalGrid requires the angular origin to be within [0, 2π]."))
        0 <= θ_stop <= 2π || throw(ArgumentError("CylindricalGrid requires the angular range to stay within [0, 2π]."))
        new{C, Dim, T}(spacing, origin, lengths)
    end
end
function area(grid::CylindricalGrid{C,2}, ind::Integer) where C
    pos = centroid(grid, ind)
    grid.spacing[1] * grid.spacing[2] * pos[1]
end

function volume(grid::CylindricalGrid{C,3}, ind::Union{Integer, CartesianIndex}) where C
    pos = centroid(grid, ind)
    grid.spacing[1] * grid.spacing[2] * pos[1] * grid.spacing[3]
end

function centroid(grid::CylindricalGrid{C,3}, ind::CartesianIndex{3}) where C
    @boundscheck checkbounds(grid, ind)
    vec3 = grid.origin .+ grid.spacing .* (ind.I .- 1)
    C(vec3)
end

function change_coordinate_type(::Type{C_New}, grid::CartesianGrid) where C_New
    CartesianGrid(grid.lengths, C_New(grid.origin), grid.spacing) 
end,
function change_coordinate_type(::Type{C_New}, grid::CylindricalGrid) where C_New
    CylindricalGrid(grid.lengths, C_New(grid.origin), grid.spacing) 
end

integration_space(grid::CylindricalGrid{C, 2}, ind) where C = area(grid, ind)
integration_space(grid::CylindricalGrid{C, 3}, ind) where C = volume(grid, ind)
integration_space(grid::CartesianGrid{C, 2}, ind) where C = area(grid, ind)
integration_space(grid::CartesianGrid{C, 3}, ind) where C = volume(grid, ind)

centroid(grid::Domain, ind::Integer)= centroid(grid, CartesianIndices(grid.lengths)[ind])
nelements(grid::Domain) = prod(grid.lengths)

Base.eachindex(grid::Domain) = Base.OneTo(nelements(grid))
Base.axes(grid::Domain, ind) = Base.OneTo(size(grid, ind))

Base.size(grid::Domain) = grid.lengths
Base.size(grid::Domain, ind::Integer) = grid.lengths[ind]

function Base.isapprox(grid::CylindricalGrid{C, 3}, grid2::CylindricalGrid{C, 3}; kwargs...) where C
    all(isapprox.(grid.spacing, grid2.spacing, kwargs...)) &&
    isapprox(grid.origin, grid2.origin, kwargs...) &&
    grid.lengths == grid2.lengths
end

function Base.isapprox(grid::CartesianGrid{C, 3}, grid2::CartesianGrid{C, 3}; kwargs...) where C
    all(isapprox.(grid.spacing, grid2.spacing; kwargs...)) &&
    isapprox(grid.origin, grid2.origin; kwargs...) &&
    grid.lengths == grid2.lengths
end

function Base.getindex(grid::CylindricalGrid{C, 3}, ::Colon, ::Colon, ind::Int) where C
    @boundscheck 1 <= ind <= grid.lengths[3] 
    new_origin = C(grid.origin .+ (0, 0, grid.spacing[3] * (ind - 1)))
    CylindricalGrid(grid.lengths[1:2], new_origin, grid.spacing)
end

function get_mesh(::Type{<:CartesianGrid}, x::AbstractRange, y::AbstractRange, z::AbstractRange)
    xyz = (x, y, z)
    lengths = length.(xyz)
    spacing = step.(xyz)
    start = first.(xyz)
    CartesianGrid(lengths, Point(start), spacing)
end

function get_mesh(::Type{<:CylindricalGrid}, r::AbstractRange, θ::AbstractRange, z::AbstractRange)
    rθz = (r, θ, z)
    lengths = length.(rθz)
    spacing = step.(rθz)
    start = first.(rθz)
    CylindricalGrid(lengths, Point(start), spacing)
end

struct Scale{Dim, T}
    factors::NTuple{Dim, T}
    function Scale(factors::NTuple{Dim, T}) where {Dim, T}
        new{Dim, T}(factors)
    end
end
# Scale(factors...) = Scale(factors)

function (scale::Scale{N})(mesh::CylindricalGrid{C, N}) where {N,C}
    CylindricalGrid(mesh.lengths, C(mesh.origin .* scale.factors), mesh.spacing .* scale.factors)
end

function (scale::Scale{N})(mesh::CartesianGrid{C, N}) where {C,N} 
    CartesianGrid(mesh.lengths, C(mesh.origin .* scale.factors), mesh.spacing .* scale.factors)
end

for domain in (:CartesianGrid, :CylindricalGrid)
    eval(quote
        function Base.view(grid::$domain{C,3,T}, ::Colon, ::Colon, ind::UnitRange) where {C,T}
            # @boundscheck 1 <= ind <= grid.lengths[3] 
            # Needs to do bound check
            $domain((grid.lengths[1:2]..., length(ind)), 
                    C(grid.origin[1:2]..., grid.origin[3] + grid.spacing[3] * (first(ind) - 1)),
                    (grid.spacing[1:2]..., grid.spacing[3]))
        end
        function Base.getindex(grid::$domain{C,3,T}, x, y, z) where {C,T}
            @view grid[grid, x, y, z]
        end
    end)
end

volume(grid::CartesianGrid{3}, ind) = prod(grid.spacing)
volume(grid::Domain, ind) = Meshes.volume(element(grid, ind))

measure(grid::CartesianGrid, ind) = prod(grid.spacing)
measure(grid::Domain, ind) = Meshes.measure(element(grid, ind))

struct CylindricalGrid{Dim,T}
    spacing::NTuple{Dim, T}
    origin::Point{Dim, T}
    lengths::NTuple{Dim, Int}
    function CylindricalGrid(lengths::NTuple{Dim, Int}, origin::Point{Dim, T}, spacing::NTuple{Dim, T}) where {Dim, T}
        new{Dim, T}(spacing, origin, lengths)
    end
end

function volume(grid::CylindricalGrid{3}, ind)
    pos = centroid(grid, ind)
    prod(grid.spacing) * pos.coords[1]
end

function Meshes.centroid(grid::CylindricalGrid{3}, ind::CartesianIndex{3})
    vec3 = grid.origin.coords .+ grid.spacing .* (ind.I .- 1 ./ 2)
    Point(vec3.coords)
end
Meshes.centroid(grid::CylindricalGrid{3}, ind::Integer) = centroid(grid, CartesianIndices(grid.lengths)[ind])

Base.eachindex(grid::CylindricalGrid{3}) = 1:prod(grid.lengths)
Base.size(grid::CylindricalGrid{3}) = grid.lengths

function Base.isapprox(grid::CylindricalGrid{3}, grid2::CylindricalGrid{3}; kwargs...)
    all(isapprox.(grid.spacing, grid2.spacing, kwargs...)) &&
    isapprox(grid.origin, grid2.origin, kwargs...) &&
    grid.lengths == grid2.lengths
end

function Base.getindex(grid::CylindricalGrid{3}, ::Colon, ::Colon, ind::Int)
    @boundscheck 1 <= ind <= grid.lengths[3] 
    new_origin = Point(grid.origin.coords .+ (0, 0, grid.spacing[3] * (ind - 1)))
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
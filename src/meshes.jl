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
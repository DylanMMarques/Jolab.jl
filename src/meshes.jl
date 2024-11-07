abstract type Domain{C, Dim, T} end

struct Point{C, Dim, T}
    coords::SVector{Dim, T}
    function Point(::Type{C}, coords::SVector{Dim, T}) where {Dim, T, C}
        new{C, Dim, T}(coords)
    end
end
function Base.getindex(p::Point, i)
    @boundscheck checkbounds(p.coords, i)
    @inbounds p.coords[i]
end
Point(::Type{C}, coords::NTuple) where C = Point(C, SVector(coords...))

Base.isapprox(p1::Point, p2::Point; kwargs...) = isapprox(p1.coords, p2.coords, kwargs...)

abstract type CoordType end
struct X_Y_λ <: CoordType end
struct X_Y_t <: CoordType end
struct X_NSY_λ <: CoordType end
struct X_NSY_t <: CoordType end
struct NSX_Y_λ <: CoordType end
struct NSX_Y_t <: CoordType end
struct NSX_NSY_λ <: CoordType end
struct NSX_NSY_t <: CoordType end
struct R_θ_λ <: CoordType end
struct R_θ_t <: CoordType end
struct NSR_NSθ_λ <: CoordType end
struct NSR_NSθ_t <: CoordType end

for (cart, polar) in zip((:NSX_NSY_λ, :NSX_NSY_t, :X_Y_λ, :X_Y_t), (:NSR_NSθ_λ, :NSR_NSθ_t, :R_θ_λ, :R_θ_t))
    eval(quote
        function $cart(coords::Point{$polar, 3})
            car = CartesianFromPolar()(Polar(coords[1], coords[2]))
            Point($cart, (car.x, car.y, coords[3]))
        end
        function $polar(coords::Point{$cart, 3})
            car = PolarFromCartesian()(coords[1:2])
            Point($polar, (car.r, car.θ, coords[3]))
        end
    end)
end

for i in (:X_Y_λ, :X_Y_t, :X_NSY_λ, :X_NSY_t, :NSX_Y_λ, :NSX_Y_t, :NSX_NSY_λ, :NSX_NSY_t, :R_θ_λ, :R_θ_t, :NSR_NSθ_λ, :NSR_NSθ_t)
    eval(quote
        $i(coords::Point{$i, 3}) = coords
    end)
end

struct CartesianGrid{C, Dim, T} <: Domain{C, Dim,T}
    spacing::NTuple{Dim, T}
    origin::Point{C, Dim, T}
    lengths::NTuple{Dim, Int}
    function CartesianGrid(lengths::NTuple{Dim, Int}, origin::Point{C, Dim, T}, spacing::NTuple{Dim, T}) where {Dim, T, C}
        new{C, Dim, T}(spacing, origin, lengths)
    end
end

area(grid::CartesianGrid{C, 2}, ind::Integer) where C = grid.spacing[1] * grid.spacing[2]
volume(grid::CartesianGrid{C, 3}, ind::Integer) where C = grid.spacing[1] * grid.spacing[2] * grid.spacing[3]
centroid(grid::CartesianGrid{C, Dim}, ind::CartesianIndex{Dim}) where {C, Dim} = Point(C, grid.origin.coords .+ grid.spacing .* (ind.I .- 1 ./ 2))

struct CylindricalGrid{C, Dim,T} <: Domain{C, Dim, T}
    spacing::NTuple{Dim, T}
    origin::Point{C, Dim, T}
    lengths::NTuple{Dim, Int}
    function CylindricalGrid(lengths::NTuple{Dim, Int}, origin::Point{C, Dim, T}, spacing::NTuple{Dim, T}) where {C, Dim, T}
        new{C, Dim, T}(spacing, origin, lengths)
    end
end
function area(grid::CylindricalGrid{C,2}, ind::Integer) where C
    pos = centroid(grid, ind)
    grid.spacing[1] * grid.spacing[2] * pos.coords[1]
end

function volume(grid::CylindricalGrid{C,3}, ind::Integer) where C
    pos = centroid(grid, ind)
    grid.spacing[1] * grid.spacing[2] * pos.coords[1] * grid.spacing[3]
end

function centroid(grid::CylindricalGrid{C,3}, ind::CartesianIndex{3}) where C
    vec3 = grid.origin.coords .+ grid.spacing .* (ind.I .- 1 ./ 2)
    Point(C, vec3)
end

centroid(grid::Domain, ind::Integer)= centroid(grid, CartesianIndices(grid.lengths)[ind])
nelements(grid::Domain) = prod(grid.lengths)

Base.eachindex(grid::Domain) = Base.OneTo(nelements(grid))
Base.size(grid::Domain) = grid.lengths
Base.size(grid::Domain, ind::Integer) = grid.lengths[ind]

function Base.isapprox(grid::CylindricalGrid{C, 3}, grid2::CylindricalGrid{C, 3}; kwargs...) where C
    all(isapprox.(grid.spacing, grid2.spacing, kwargs...)) &&
    isapprox(grid.origin, grid2.origin, kwargs...) &&
    grid.lengths == grid2.lengths
end

function Base.isapprox(grid::CartesianGrid{C, 3}, grid2::CartesianGrid{C, 3}; kwargs...) where C
    all(isapprox.(grid.spacing, grid2.spacing, kwargs...)) &&
    isapprox(grid.origin, grid2.origin, kwargs...) &&
    grid.lengths == grid2.lengths
end

function Base.getindex(grid::CylindricalGrid{C, 3}, ::Colon, ::Colon, ind::Int) where C
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

struct Scale{Dim, T}
    factors::NTuple{Dim, T}
    function Scale(factors::NTuple{Dim, T}) where {Dim, T}
        new{Dim, T}(factors)
    end
end
# Scale(factors...) = Scale(factors)

function (scale::Scale{N})(mesh::CylindricalGrid{C, N}) where {N,C}
    CylindricalGrid(mesh.lengths, Point(C, mesh.origin.coords .* scale.factors), mesh.spacing .* scale.factors)
end

function (scale::Scale{N})(mesh::CartesianGrid{C, N}) where {C,N} 
    CartesianGrid(mesh.lengths, Point(C, mesh.origin.coords .* scale.factors), mesh.spacing .* scale.factors)
end

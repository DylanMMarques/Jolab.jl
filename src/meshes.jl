abstract type Domain{C, Dim, T} end

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

for (cart, polar) in zip((:NSX_NSY_λ, :X_Y_λ, :X_Y_t), (:NSR_NSθ_λ, :R_θ_λ))
    eval(quote
        function $cart(coords::$polar{T}) where T
            car = CartesianFromPolar()(Polar(coords[1], coords[2]))
            $cart(car.x, car.y, coords[3])
        end
        function $polar(coords::$cart{T}) where T
            car = PolarFromCartesian()(coords[1:2])
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

area(grid::CartesianGrid{C, 2}, ind::Integer) where C = grid.spacing[1] * grid.spacing[2]
volume(grid::CartesianGrid{C, 3}, ind::Integer) where C = grid.spacing[1] * grid.spacing[2] * grid.spacing[3]
centroid(grid::CartesianGrid{C, Dim}, ind::CartesianIndex{Dim}) where {C, Dim} = C(grid.origin .+ grid.spacing .* (ind.I .- 1 ./ 2))

struct CylindricalGrid{C, Dim,T} <: Domain{C, Dim, T}
    spacing::NTuple{Dim, T}
    origin::C
    lengths::NTuple{Dim, Int}
    function CylindricalGrid(lengths::NTuple{Dim, Int}, origin::C, spacing::NTuple{Dim, T}) where {C, Dim, T}
        new{C, Dim, T}(spacing, origin, lengths)
    end
end
function area(grid::CylindricalGrid{C,2}, ind::Integer) where C
    pos = centroid(grid, ind)
    grid.spacing[1] * grid.spacing[2] * pos[1]
end

function volume(grid::CylindricalGrid{C,3}, ind::Integer) where C
    pos = centroid(grid, ind)
    grid.spacing[1] * grid.spacing[2] * pos[1] * grid.spacing[3]
end

function centroid(grid::CylindricalGrid{C,3}, ind::CartesianIndex{3}) where C
    vec3 = grid.origin .+ grid.spacing .* (ind.I .- 1 ./ 2)
    C(vec3)
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

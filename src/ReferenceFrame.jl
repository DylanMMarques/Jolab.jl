export ReferenceFrame

struct ReferenceFrame{T}
    origin::Point{X_Y_Z, 3, T}
    direction::Point{X_Y_Z, 3, T}
    function ReferenceFrame(::Type{T}, origin::Point{X_Y_Z, 3}, direction::Point{X_Y_Z, 3}) where T
        new{T}(origin, direction)
    end
end

function ReferenceFrame(::Type{T}, origin, direction) where T
    ReferenceFrame(T, Point(X_Y_Z, origin), Point(X_Y_Z, direction))
end
ReferenceFrame(origin, direction) = ReferenceFrame(Float64, origin, direction)

Base.convert(::Type{ReferenceFrame{T}}, frame::ReferenceFrame) where T = ReferenceFrame{T}(frame.origin, frame.direction)

function Base.isapprox(frame1::ReferenceFrame, frame2::ReferenceFrame; kwargs...)
    isapprox(frame1.origin, frame2.origin; kwargs...) && isapprox(frame1.direction, frame2.direction; kwargs...)
end

RotXYZ(direction::FieldVector{3}) = Rotations.RotXYZ(direction.x, direction.y, direction.z)
RotXYZ(x, y, z) = Rotations.RotXYZ(x, y, z)
RotXYZ(p::Point{X_Y_Z, 3}) = Rotations.RotXYZ(p.coords[1], p.coords[2], p.coords[3])
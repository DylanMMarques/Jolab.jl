export ReferenceFrame
struct Point3D{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T 
end

struct Point2D{T} <: FieldVector{2,T}
    x::T
    y::T
end

struct ReferenceFrame{T}
    origin::Point3D{T}
    direction::Point3D{T}
end

function ReferenceFrame(::Type{T}, origin, direction) where T
    ReferenceFrame{T}(origin, direction)
end
ReferenceFrame(origin, direction) = ReferenceFrame(Float64, origin, direction)

Base.convert(::Type{ReferenceFrame{T}}, frame::ReferenceFrame) where T = ReferenceFrame{T}(frame.origin, frame.direction)

Base.isapprox(frame1::ReferenceFrame, frame2::ReferenceFrame; kwargs...) = isapprox(frame1.origin, frame2.origin; kwargs...) && isapprox(frame1.direction, frame2.direction; kwargs...)
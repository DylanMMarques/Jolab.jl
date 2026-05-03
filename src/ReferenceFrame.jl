export ReferenceFrame

struct ReferenceFrame{T}
    origin::X_Y_Z{T}
    direction::X_Y_Z{T}
    function ReferenceFrame(::Type{T}, origin::X_Y_Z, direction::X_Y_Z) where {T}
        new{T}(origin, direction)
    end
end

function ReferenceFrame(::Type{T}, origin, direction) where {T}
    ReferenceFrame(T, X_Y_Z(origin), X_Y_Z(direction))
end
ReferenceFrame(origin, direction) = ReferenceFrame(Float64, origin, direction)

Base.convert(::Type{ReferenceFrame{T}}, frame::ReferenceFrame) where {T} =
    ReferenceFrame(T, convert(X_Y_Z{T}, frame.origin), convert(X_Y_Z{T}, frame.direction))

function Base.isapprox(frame1::ReferenceFrame, frame2::ReferenceFrame; kwargs...)
    isapprox(frame1.origin, frame2.origin; kwargs...) &&
        isapprox(frame1.direction, frame2.direction; kwargs...)
end

function Base.:(+)(frame1::ReferenceFrame, frame2::ReferenceFrame)
    if !isapprox(frame1.direction, frame2.direction)
        throw(
            ArgumentError(
                "Cannot add ReferenceFrames with different directions. frame1.direction = $(frame1.direction), frame2.direction = $(frame2.direction)",
            ),
        )
    end
    ReferenceFrame(frame1.origin + frame2.origin, frame1.direction)
end

_RotXYZ(direction::FieldVector{3}) = Rotations.RotXYZ(direction.x, direction.y, direction.z)
_RotXYZ(x, y, z) = Rotations.RotXYZ(x, y, z)

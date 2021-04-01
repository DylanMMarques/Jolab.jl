mutable struct ReferenceFrame{T<:Real}
	x::T
	y::T
	z::T
	θ::T
	ϕ::T
	ReferenceFrame{T}(x, y, z, θ = 0, ϕ = 0) where T = new{T}(x, y, z, θ, ϕ)
end
ReferenceFrame(x, y, z, θ = 0, ϕ = 0) = ReferenceFrame{Float64}(x, y, z, θ, ϕ)
Base.eltype(ref::ReferenceFrame{T}) where T = T

Base.copy(ref::ReferenceFrame{T}) where T = ReferenceFrame{T}(copy(ref.x), copy(ref.y), copy(ref.z), copy(ref.θ), copy(ref.ϕ))

function Base.:(==)(ref1::ReferenceFrame, ref2::ReferenceFrame)::Bool
	return (abs(ref1.x - ref2.x) < @tol) && (abs(ref1.y - ref2.y) < @tol) && (abs(ref1.z - ref2.z) < @tol) && (abs(ref1.θ - ref2.θ) < @tol) && (abs(ref1.ϕ - ref2.ϕ) < @tol);
end

function Base.:(!=)(ref1::ReferenceFrame, ref2::ReferenceFrame)::Bool
	return !(ref1 == ref2);
end

function Base.:(+)(ref1::ReferenceFrame{A}, ref2::ReferenceFrame{B}) where {A,B}
	T = promote_type(A,B)
	checkorientation(ref1, ref2) && (return ReferenceFrame{T}(ref1.x .+ ref2.x, ref1.y .+ ref2.y, ref1.z .+ ref2.z, ref1.θ, ref1.ϕ))
	error("Sum of referenceframe is only valid for referenceframes with the same orientation (θ and ϕ)")
end

function checkorientation(ref1::ReferenceFrame, ref2::ReferenceFrame)::Bool
	return (abs(ref1.θ - ref2.θ) < @tol) && (abs(ref1.ϕ - ref2.ϕ) < @tol);
end

function checkinline(ref1::ReferenceFrame, ref2::ReferenceFrame)::Bool
	checkorientation(ref1, ref2) || return false
	if ref1.θ > @tol
		t = (ref2.z - ref1.z) / cos(ref1.θ)
		(abs(ref1.x + t * sin(ref1.θ) * cos(ref1.ϕ) - ref2.x) < @tol) || return false
		(abs(ref1.y + t * sin(ref1.θ) * sin(ref1.ϕ) - ref2.y) < @tol) || return false
	else
		(abs(ref1.x - ref2.x) < @tol) || return false
		(abs(ref1.y - ref2.y) < @tol) || return false
	end
	return true
end

function checkposition(ref1::ReferenceFrame, ref2::ReferenceFrame)::Bool
	return (abs(ref1.x - ref2.x) < @tol) && (abs(ref1.y - ref2.y) < @tol) && (abs(ref1.z - ref2.z) < @tol);
end

function checkinplane(ref1::ReferenceFrame, ref2::ReferenceFrame)
	# check if ref2 is in the plane xy of ref1

	# checkorientation(ref1, ref2) || return false
	(a1, b1, c1) = rotatecoordinatesto(0, 0, 1, ref1.θ, ref1.ϕ)
	d1 = a1 * ref1.x + b1 * ref1.y + c1 * ref1.z;

	# consider the orientation of ref1 as we just want to see if the are in plane
	d2 = a1 * ref2.x + b1 * ref2.y + c1 * ref2.z;
	(abs(d1 - d2) < @tol) || return false
	return true
end

function rotatecoordinatesto(sx, sy, sz, θ, ϕ)
	sxr = sx * cos(θ) * cos(ϕ) + sz * cos(ϕ) * sin(θ) - sy * sin(ϕ)
	syr = sy * cos(ϕ) + (sx * cos(θ) + sz * sin(θ)) * sin(ϕ)
	szr = sz * cos(θ) - sx * sin(θ)
	return (sxr, syr, szr)
end

function rotatecoordinatesfrom(sx, sy, sz, θ, ϕ)
	sxr = -sz * sin(θ) + cos(θ) * (sx * cos(ϕ) + sy * sin(ϕ))
	syr = sy * cos(ϕ) - sx * sin(ϕ)
	szr = sz * cos(θ) + sin(θ) * (sx * cos(ϕ) + sy * sin(ϕ));
	return (sxr, syr, szr);
end

@inline function rotatecoordinatesfrom(sx::T, sy::T, sz::T, θ::Real, ϕ::Real)::Tuple{T,T,T} where T
	sxr = -real(sz) * sin(θ) + cos(θ) * (sx * cos(ϕ) + sy * sin(ϕ))
	syr = sy * cos(ϕ) - sx * sin(ϕ)
	szr = real(sz) * cos(θ) + sin(θ) * (sx * cos(ϕ) + sy * sin(ϕ))
	return (sxr, syr, szr)
end

@inline function rotatecoordinatesto!(sx::AbstractArray{<:Number}, sy::AbstractArray{<:Number}, sz::AbstractArray{<:Number}, θ::Real, ϕ::Real)
	length(sx) == length(sy) == length(sz) || error("Arrays must have the same length")
	@inbounds @simd for i in eachindex(sxr)
		tmpsx = sx[i] * cos(θ) * cos(ϕ) + real(sz[i]) * cos(ϕ) * sin(θ) - sy[i] * sin(ϕ);
		tmpsy = sy[i] * cos(ϕ) + (sx[i] * cos(θ) + real(sz[i]) * sin(θ)) * sin(ϕ);
		sz[i] = real(sz[i]) * cos(θ) - sx[i] * sin(θ);
		sx[i] = tmpsx;
		sy[i] = tmpsx;
	end
end

@inline function rotatecoordinatesfrom!(sx::AbstractArray{<:Number}, sy::AbstractArray{<:Number}, sz::AbstractArray{<:Number}, θ::Real, ϕ::Real)
	length(sx) == length(sy) == length(sz) || error("Arrays must have the same length")
	@inbounds @simd for i in eachindex(sx)
		tmpsx = -real(sz[i]) * sin(θ) + cos(θ) * (sx[i] * cos(ϕ) + sy[i] * sin(ϕ));
		tmpsy = sy[i] * cos(ϕ) - sx[i] * sin(ϕ);
		sz[i] = real(sz[i]) * cos(θ) + sin(θ) * (sx[i] * cos(ϕ) + sy[i] * sin(ϕ));
		sx[i] = tmpsx;
		sy[i] = tmpsy;
	end
end

function rotatecoordinatesfromto!(sx::AbstractArray{<:Number}, sy::AbstractArray{<:Number}, sz::AbstractArray{<:Number}, θi::Real, ϕi::Real, θf::Real, ϕf::Real)
	length(sx) == length(sy) == length(sz) || error("Arrays must have the same length")
	@inbounds @simd for i in eachindex(sx)
		(sx[i], sy[i], sz[i]) = rotatecoordinatesfromto(sx[i], sy[i], sz[i], θi, ϕi, θf, ϕf)
	end
end

function rotatecoordinatesfromto(sx, sy, sz, θi, ϕi, θf, ϕf)
	tmpsx = -sz * sin(θi) + cos(θi) * (sx * cos(ϕi) + sy * sin(ϕi));
	tmpsy = sy * cos(ϕi) - sx * sin(ϕi);
	tmpsz = sz * cos(θi) + sin(θi) * (sx * cos(ϕi) + sy * sin(ϕi));

	sx2 = tmpsx * cos(θf) * cos(ϕf) + tmpsz * cos(ϕf) * sin(θf) - tmpsy * sin(ϕf);
	sy2 = tmpsy * cos(ϕf) + (tmpsx * cos(θf) + tmpsz * sin(θf)) * sin(ϕf);
	sz2 = tmpsz * cos(θf) - tmpsx * sin(θf)
	return (sx2, sy2, sz2)
end

@inline rotatecoordinatesfromto(sx, sy, sz, refold::ReferenceFrame, refnew::ReferenceFrame) = rotatecoordinatesfromto(sx, sy, sz, refold.θ, refold.ϕ + π, refnew.θ, refnew.ϕ + π);

@inline function rotatecoordinatesfromto!(sx::AbstractArray{<:Number}, sy::AbstractArray{<:Number}, sz::AbstractArray{<:Number}, refold::ReferenceFrame, refnew::ReferenceFrame)
	return rotatecoordinatesfromto!(sx, sy, sz, -refold.θ, -refold.ϕ, -refnew.θ, -refnew.ϕ);
end

distance(ref1::ReferenceFrame, ref2::ReferenceFrame) = √((ref1.x - ref2.x)^2 + (ref1.y - ref2.y)^2 + (ref1.z - ref2.z)^2)

Base.convert(::Type{ReferenceFrame{A}}, ref::ReferenceFrame{B}) where {A<:Real, B<:Real} = return ReferenceFrame{A}(A(ref.x), A(ref.y), A(ref.z), A(ref.θ), A(ref.ϕ))

function planelineintersection(refPlane::ReferenceFrame, refLine::ReferenceFrame)
	# equation of the plane ax + by + cy = d
	(a, b, c) = rotatecoordinatesto(0, 0, 1, refPlane.θ, refPlane.ϕ)
	d = a * refPlane.x + b * refPlane.y + c * refPlane.z;

	# line definition (x,y,z) = (x0,y0,z0) + t(Δx,Δy,Δz)
	(Δx, Δy, Δz) = rotatecoordinatesto(0, 0, 1, refLine.θ, refLine.ϕ)
	# intersection with the plane of refPlane
	t = (d - a * refLine.x - b * refLine.y - c * refLine.z) / (a * Δx + b * Δy + c * Δz)
	return ReferenceFrame(refLine.x + t * Δx, refLine.y + t * Δy, refLine.z + t * Δz, refLine.θ, refLine.ϕ)
end

function rotatereferenceframe!(ref_ref::ReferenceFrame, ref::ReferenceFrame, θ::Real, ϕ::Real)
	(ref.x, ref.y, ref.z) = rotatecoordinatesto(ref.x - ref_ref.x, ref.y - ref_ref.y, ref.z - ref_ref.z, θ, ϕ)
	ref.x += ref_ref.x
	ref.y += ref_ref.y
	ref.z += ref_ref.z

	(sx, sy, sz) = rotatecoordinatesto(sin(ref.θ) * cos(ref.ϕ), sin(ref.θ) * sin(ref.ϕ), cos(ref.θ), θ, ϕ)
	ref.θ = acos(sz)
	ref.ϕ = atan(sy, sx)
end

function translatereferenceframe!(ref::ReferenceFrame, x::Real, y::Real, z::Real)
	ref.x += x
	ref.y += y
	ref.z += z
end

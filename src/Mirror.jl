struct Mirror{T} <: AbstractOpticalElement{T}
    r
    media
    n2::Medium
    frame::ReferenceFrame{T}
end

function Mirror(::Type{T}, r, frame)
    Mirror{T}(r, frame)
end
Mirror(r, frame) = Mirror(Float64, r, frame)

## Need to be done
reflectivity(mirror::Mirror, λ, nsx, nsy) = mirror.r(λ, nsx, nsy)

function light_interaction(m::Mirror, pw::PlaneWaveScalar{T}) where T
    (r, t) = rtss(mirror, dir, pw.λ, pw.nsx, pw.nsy)
    rpw = PlaneWaveScalar(pw.nsx, pw.nsy, pw.e * r, pw.λ, mirror.media[2], mirror.frame)
    tpw = PlaneWaveScalar(pw.nsx, pw.nsy, pw.e * t, pw.λ, mirror.media[1], mirror.frame)
    (rpw, tpw)
end

function rtss(mirror::Mirror, dir::AbstractDirection, nsx, nsy, λ) 
    ind_media = dir <: Forward ? (1, 2) : (2, 1)
    n1 = √(mirror.media[ind_media[1]].n(λ))
    n2 = √(mirror.media[ind_media[2]].n(λ))
    tλ = √(1 - reflectivity(mirror, λ, nsx, nsy)) * n1 / n2
    rλ = √reflectivity(mirror, λ, nsx, nsy) * (dir <: Forward ? 1 : -1)
    (rλ, tλ)
end

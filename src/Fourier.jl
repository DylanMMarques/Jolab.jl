struct Fourier{T,S} <: AbstractOpticalElement{T}
    frames::Tuple{ReferenceFrame{T}, ReferenceFrame{T}}
    solver::S
    function Fourier(::Type{T}, frames::Tuple{ReferenceFrame{T}, ReferenceFrame{T}}; solver::S = nothing) where {T,S}
        new{T,S}(frames, solver)
    end
end
Fourier(frames; kwargs...) = Fourier(Float64, frames; kwargs...)

@inline function t(::Type{<:MeshedBeam}, fourier::Fourier, coord_in::NSR_NSθ_λ, coord_out::R_θ_λ, area)
    (nsr, nsθ, λ) = coord_in.coords
    (r, θ, tmp) = coord_out.coords
    im * 2π * (2π / λ)^2 * besselj0(2π / λ * nsr * r) * area
end

@inline function t(::Type{<:MeshedBeam}, fourier::Fourier, coord_in::R_θ_λ, coord_out::R_θ_λ, area)
    (nsr, nsθ, λ) = coord_out.coords
    (r, θ, tmp) = coord_in.coords
    - im / 2π * besselj0(2π / λ * nsr * r) * area
end

check_input_field(f::Fourier, field_i::MeshedBeam) = zero(UInt64)
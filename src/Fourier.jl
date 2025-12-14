struct Fourier{T,S} <: AbstractOpticalElement{T}
    frames::Tuple{ReferenceFrame{T}, ReferenceFrame{T}}
    solver::S
    function Fourier(::Type{T}, frames::Tuple{ReferenceFrame{T}, ReferenceFrame{T}}; solver::S = nothing) where {T,S}
        new{T,S}(frames, solver)
    end
end
Fourier(frames; kwargs...) = Fourier(Float64, frames; kwargs...)


struct FourierFFT end

function _light_interaction!(field_b::MeshedBeam{<:Any,Backward}, field_f::MeshedBeam{<:Any,Forward}, fourier::FourierFFT, field_i::MeshedBeam{T,D,NSX_NSY_λ}) where {T,D}
    (field_r, field_t) = reverse_if_backward(D, (field_b, field_f))
    fill!(field_r.e, 0)

    (nsx, nsy, _λ) = get_ranges(field_i.mesh)
    (x, y, _) = get_ranges(field_t.mesh)
    λ = only(_λ)

    tmp = field_i.e .* exp.((im * 2π / λ) .* (nsx .* first(x) .+ (nsy)' .* first(y)))
        
    ifftshift!(field_t.e, tmp, (1,2))

    bfft!(field_t.e, (1,2))

    field_t.e .*= (step(nsx) * step(nsy) * (2π / λ)^2)

    (field_b, field_f)
end

function _light_interaction!(field_b::MeshedBeam{<:Any,Backward}, field_f::MeshedBeam{<:Any,Forward}, fourier::FourierFFT, field_i::MeshedBeam{T,D,X_Y_λ}) where {T,D}
    (field_r, field_t) = reverse_if_backward(D, (field_b, field_f))
    fill!(field_r.e, 0)

    tmp = fft(field_i.e, (1,2))
    fftshift!(field_t.e, tmp, (1,2))

    (x, y, _λ) = get_ranges(field_i.mesh)
    (nsx, nsy, _) = get_ranges(field_t.mesh)
    λ = only(_λ)

    field_t.e .*= exp.(-(im * 2π / λ) .* (first(x) .* nsx .+ first(y) .* (nsy)')) .*
            (step(x) * step(y) / (2π)^2)

    (field_b, field_f)
end

function check_input_field(f::FourierFFT, field_i::MeshedBeam{T,D,C}) where {T,D,C}
    if !(C <: X_Y_λ || C <: NSX_NSY_λ)
        throw(ArgumentError("FourierFFT only supports X_Y_λ and NSX_NSY_λ coordinate types"))
    end
    zero(UInt64)
end

function forward_backward_field(fourier::FourierFFT, field_i::MeshedBeam{T,D,C}) where {T,D,C}
    e_t = similar(field_i.e, Complex{T})
    e_r = Zeros(T, size(field_i.e))

    mesh = if C <: AngularSpectrumCoords
        if C <: NSX_NSY_λ
           (nsx, nsy, λ) = get_ranges(field_i.mesh)
            dkx = step(nsx) / only(λ)
            dky = step(nsy) / only(λ)
            x = fftshift(fftfreq(length(nsx), 1 / dkx))
            y = fftshift(fftfreq(length(nsy), 1 / dky))

            CartesianGrid(size(field_i.mesh),
                          X_Y_λ(first(x), first(y), first(λ)),
                          (step(x), step(y), field_i.mesh.spacing[3]))
        else
            error("Not implemented for R_θ_λ and NSR_NSθ_λ coordinates")
        end
    else
        if C <: X_Y_λ
            (x, y, λ) = get_ranges(field_i.mesh)
            nsx = fftshift(fftfreq(length(x), 1 / step(x))) * only(λ)
            nsy = fftshift(fftfreq(length(y), 1 / step(y))) * only(λ)
            CartesianGrid(size(field_i.mesh),
                NSX_NSY_λ(first(nsx), first(nsy), first(λ)),
                (step(nsx), step(nsy), field_i.mesh.spacing[3]))
        else
            error("Not implemented for R_θ_λ and NSR_NSθ_λ coordinates")
        end
    end

    (dir_r, dir_t) = reverse_if_backward(D, (Backward, Forward))
    field_r = MeshedBeam{T, dir_r, C}(field_i.mesh, e_r, field_i.medium, field_i.frame)
    field_t = MeshedBeam{T, dir_t, get_transmitted_coord_type(FourierFFT, C)}(mesh, e_t, field_i.medium, field_i.frame)
    reverse_if_backward(D, (field_r, field_t))
end

function get_transmitted_coord_type(::Type{FourierFFT}, ::Type{NSX_NSY_λ}) 
    X_Y_λ
end,
function get_transmitted_coord_type(::Type{FourierFFT}, ::Type{NSR_NSθ_λ}) 
    R_θ_λ
end,
function get_transmitted_coord_type(::Type{FourierFFT}, ::Type{X_Y_λ})
    NSX_NSY_λ
end,
function get_transmitted_coord_type(::Type{FourierFFT}, ::Type{R_θ_λ})
    NSR_NSθ_λ
end


@inline function t(::Type{<:MeshedBeam}, fourier::Fourier, coord_in::NSR_NSθ_λ, coord_out::R_θ_λ, area)
    (nsr, nsθ, λ) = coord_in
    (r, θ, tmp) = coord_out
    im * 2π * (2π / λ)^2 * besselj0(2π / λ * nsr * r) * area
end

# @inline function t(::Type{<:MeshedBeam}, fourier::Fourier, coord_in::R_θ_λ, coord_out::R_θ_λ, area)
#     (r, θ, tmp) = coord_in
#     (nsr, nsθ, λ) = coord_out
#     - im / 2π * besselj0(2π / λ * nsr * r) * area
# end

@inline function t(::Type{<:MeshedBeam{T}}, fourier::Fourier, coord_in::NSX_NSY_λ, coord_out::X_Y_λ, area) where T
    (nsx, nsy, λ) = coord_in
    (x, y, λ) = coord_out
    exp(im * T(2π) / λ * (nsx * x + nsy * y)) * area
end

@inline function t(::Type{<:MeshedBeam{T}}, fourier::Fourier, coord_in::X_Y_λ, coord_out::NSX_NSY_λ, area) where T
    (x, y, λ) = coord_in
    (nsx, nsy, λ) = coord_out
    exp(im * T(2π) / λ * (nsx * x + nsy * y)) / T(4π^2) * area
end

check_input_field(f::Fourier, field_i::MeshedBeam) = zero(UInt64)

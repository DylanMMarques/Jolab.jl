struct Fibre{T, F} <: AbstractOpticalElement{T}
    modes::Float64
    frames::NTuple{2, ReferenceFrame{T}}
end

struct GaussianMode{T} <: AbstractMode{T}
    sigma_squared::T
    function GaussianMode(sigma::T) where T
        sigma_squared = sigma^2
        new{T}(sigma_squared)
    end
end

# In space
function coupling(m1::GaussianMode, beam)
    f(e, dA, x, y) = dA * e * exp(- (x^2 + y^2) / 2m1.sigma_squared)
    error("NOT TESTED")
    cons = m1.sigma_squared^2 / T(π)
    mapreduce(f, +, beam.modes.e, beam.modes.dA, beam.modes.x, beam.modes.y) * cons
end

function coupling(m1::GaussianMode, beam::ScalarAngularSpectrumBeam{T}) where T
    f(e, dA, λ, nsx, nsy) = dA * e * exp(- (2π / λ)^2 * m1.sigma_squared * 32 * (nsx^2 + nsy^2)) * (2π / λ)^2
    cons = 64 * m1.sigma_squared / π
    mapreduce(f, +, beam.modes.e, beam.modes.dA, beam.modes.wavelength, beam.modes.nsx, beam.modes.nsy) * cons
end
using Jolab

## Test mode normalization of fibre in angspe
function f(c)
    nsx = range(-1, 1, length=2000) .+ zeros(2000)'
    λ = 1550E-9
    beam = Jolab.monochromatic_angularspectrum(Float64, nsx, nsx', (nsx .* nsx)' .+ 0im, λ, Medium(1.0 + im), ReferenceFrame((0,0,0), (0,0,0)))
    beam.modes.e .= 1
    Jolab.coupling(Jolab.GaussianMode(c), beam) ≈ 2
end
@test f(10E-6)
@test f(20E-6)

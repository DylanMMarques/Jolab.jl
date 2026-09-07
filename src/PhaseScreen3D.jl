# Implementation of the phase screen method from the paper "Model-based wavefront shaping microscopy" of Abhilash Thendiyammal, Gerwin Osnabrugge, Tom Knop, and Ivo M. Vellekoop 
struct PhaseScreenSolver{F1,F2,F3,N,P,C,C2,Z} <: AbstractSolver
    field_b::F2
    field_f::F3
    field_i::F1
    n::N
    plan_fft::P
    tmp_array::C
    tmp_field_i::C
    kr_squared::C2
    z_interfaces::Z
end

function _light_interaction!(field_b::MeshedBeam{<:Any, Backward}, field_f::MeshedBeam{<:Any, Forward}, solver::PhaseScreenSolver, field_i::MeshedBeam{<:Any, D}) where D
    field_r, field_t = reverse_if_backward(D, (field_b, field_f))

    _n_slices = eachslice(solver.n, dims = 5)
    n_slices = D == Forward ? _n_slices : (@view _n_slices[end:-1:1])

    _, _, λ = get_ranges(field_i.mesh)

    (z_i, z_f) = reverse_if_backward(D, (solver.z_interfaces[1], solver.z_interfaces[end]))
    deltaz = step(solver.z_interfaces)
    scale = D == Forward ? 1im : -1im
    solver.tmp_field_i .= field_i.e .* exp.(scale .* z_i .* sqrt.(complex.((2 .* pi .* field_i.medium.n ./ λ).^2 .- solver.kr_squared))) # Translate the field from the interface to the tip of the first slice

    for slice in eachindex(n_slices)
        refractive_index_slice = n_slices[slice]
        mean_refractive_index = mean(refractive_index_slice)

        field_t.e .= exp.(0.5im .* deltaz .* sqrt.(complex.((2 .* pi .* mean_refractive_index ./ λ).^2 .- solver.kr_squared))) .* solver.tmp_field_i
        
        mul!(solver.tmp_array, solver.plan_fft, field_t.e)

        solver.tmp_array .*= exp.(im .* deltaz .* 2pi ./ λ .* (refractive_index_slice .- mean_refractive_index))
        
        ldiv!(solver.tmp_field_i, solver.plan_fft, solver.tmp_array)

        solver.tmp_field_i .*= exp.(0.5im .* deltaz .* sqrt.(complex.((2 .* pi .* mean_refractive_index ./ λ).^2 .- solver.kr_squared)))
    end
    print(z_f)
    field_t.e .= solver.tmp_field_i .* exp.(-scale .* z_f .* sqrt.(complex.((2 .* pi .* field_t.medium.n ./ λ).^2 .- solver.kr_squared))) # Translate the field back to the interface instead of the tip of the last slice

    (field_b, field_f)
end

function forward_backward_field(solver::PhaseScreenSolver, field_i::MeshedAngularSpectrum) 
    solver.field_b, solver.field_f
end


function PhaseScreenSolver(comp::RoughInterface, field_i::MeshedAngularSpectrum{T,D,C,P}, steps) where {T,D,C,P}
    n1, n2 = comp.mat[1].n, comp.mat[2].n
    
    nsx, nsy, wavelength = get_ranges(field_i.mesh)
    λ = only(wavelength)
    x = fftfreq(length(nsx), λ / step(nsx))
    y = fftfreq(length(nsy), λ / step(nsy))

    topography = comp.Δz.(x, y')

    n = similar(field_i.e, Complex{T}, (length(nsx), length(nsy), length(λ), 1, steps))

    z_s = range(minimum(topography), maximum(topography), length = steps)
    
    for (i, z_i) in enumerate(z_s)
        n[:, :, :, :, i] .= ifelse.(topography .> z_i, n2, n1)
    end

    p_fft = plan_ifft(field_i.e, (1, 2))
    tmp_array = similar(field_i.e)
    tmp_field_i = similar(field_i.e)
    kr_squared = similar(field_i.e, T)
    kr_squared .= (nsx.^2 .+ nsy'.^2) .* (2π ./ λ).^2

    e_b, e_f = reverse_if_backward(D, (Zeros(field_i.e), similar(field_i.e)))

    PhaseScreenSolver(
        MeshedBeam{T, Backward,C,P}(field_i.mesh, e_b, comp.mat[1], field_i.frame),
        MeshedBeam{T, Forward,C,P}(field_i.mesh, e_f, comp.mat[2], field_i.frame),
        field_i,
        n,
        p_fft,
        tmp_array,
        tmp_field_i,
        kr_squared,
        z_s
    )
end

function check_input_field(solver::PhaseScreenSolver, field_i::MeshedAngularSpectrum{T,D,C,P}) where {T,D,C,P}
    msg_code = zero(UInt64)
    field_i.frame ≈ solver.field_i.frame || (msg_code |= 1 << INVALID_FRAME)
    field_i.medium ≈ solver.field_i.medium || (msg_code |= 1 << INVALID_MEDIUM)
    msg_code
end

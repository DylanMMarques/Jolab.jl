abstract type AbstractWaveguideProfile{T} end
abstract type AbstractWaveguideMode{T,D} <: AbstractMode{T} end
struct Fibre{T,R, M, VM <: AbstractVector{M}}
    refractive_index_profile::R
    length::T
    frames::Tuple{ReferenceFrame{T}, ReferenceFrame{T}}
    media::Tuple{Medium{T}, Medium{T}}
    modes::Dict{T, VM}
    function Fibre(::Type{T}, profile::R, length, media, frames) where {T,R}
        M = mode_type(R)
        modes = Dict{T, StructVector{M}}()
        new{T, R, M, StructVector{M}}(profile, length, frames, media, modes)
    end
end
Fibre(profile, length, media, frames) = Fibre(Float64, profile, length, media, frames)


function check_input_field(fibre::Fibre, beam::ScalarSpatialBeam{T,D}) where {T,D}
    (n_fibre, frame_fibre) = D == Forward ? first.((fibre.media, fibre.frames)) : last.((fibre.media, fibre.frames))
    @argcheck all(broadcast((ni, λi) -> ni(λi) == n_fibre(λi), beam.modes.medium, beam.modes.wavelength)) ArgumentError("The medium of the beam and outside the fibre are not the same.")
    @argcheck all(isapprox.(Ref(frame_fibre), beam.modes.frame)) ArgumentError("The reference frame of the beam and the fibre are not the same.")
    true
end
check_input_field(fibre::Fibre, beam::Beam) = false

export Fibre, CircularStepIndexProfile

α1(na, ncore, λ, β) = √((2π / λ * ncore)^2 - β^2)
α2(na, ncore, λ, β) = √(β^2 + (2π / λ)^2 * (na^2 - ncore^2))
nclad(ncore, na) = √(ncore^2 - na^2)

function find_wavefunction_solutions(profile::AbstractWaveguideProfile{T}, m, λ) where T
    tol = T(1E-12)
    _nclad = nclad(profile.ncore.n, profile.na)
    if m > 0
        α2_min = (T(1E50) * √(T(2m / π))) ^(T(-1/m)) * 2m / T(exp(1)) / profile.radius
        βmin = √(α2_min^2 + (T(2π) * (_nclad + tol) / λ)^2);
    else
        βmin = real(2π / λ * (_nclad + tol))
    end
    βmax = real(2π / λ * (profile.ncore.n));
    
    initial_guess = 100
    initial_diff = 5000
    i = 1
    β = range(βmin, βmax, length = initial_diff * i)

    f(β) = modecondition(profile, λ, m, β)^2 # Square the function because it is easier to find the zeros

    rough_βs = find_local_minima(f, β, initial_guess)
    cur_len = length(rough_βs)
    old_length = 0
    while cur_len != old_length
        old_length = cur_len
        β = range(βmin, βmax, length = initial_diff * i)
        rough_βs = find_local_minima(f, β, initial_guess)
        cur_len = length(rough_βs)
        i += 1
    end
    map(i -> wavefunction_solutions(profile, λ, i, m), rough_βs)
end

function findmodes(profile::P, _λ) where {P<:AbstractWaveguideProfile{T}} where T
    λ = T(_λ)
    m = 0
    sizeA = 10000
    inc_size_i = 1
    β_vec = Vector{T}(undef, sizeA)
    m_vec = Vector{Int}(undef, sizeA)
    C_vec = Vector{T}(undef, sizeA)
    D_vec = Vector{T}(undef, sizeA)
    modeNumber = 0
    while true
        condition(β) = modecondition(profile, λ, m, β)
        β_solutions = find_wavefunction_solutions(profile, m, λ)
        
        isempty(β_solutions) && break
        number_solutions = length(β_solutions)
        for i in 0:(number_solutions-1)
            if modeNumber > (sizeA * inc_size_i) - 2
                map(i -> resize!(i, sizeA * inc_size_i), (β_vec, m_vec, C_vec, D_vec))
                inc_size_i += 1
            end

            βᵢ = β_solutions[i + 1]
            modeNumber += 1;
            β_vec[modeNumber] = βᵢ
            m_vec[modeNumber] = m
            (C_vec[modeNumber], D_vec[modeNumber]) = modeconstant(profile, λ, m, βᵢ)

            (m == 0) && continue
            βᵢ = β_solutions[i + 1]
            modeNumber += 1;
            β_vec[modeNumber] = βᵢ
            m_vec[modeNumber] = m
            (C_vec[modeNumber], D_vec[modeNumber]) = modeconstant(profile, λ, m, βᵢ)
        end
        m += 1
    end
    map(i -> resize!(i, modeNumber), (β_vec, m_vec, C_vec, D_vec))
    StructVector{CircularStepIndexMode{T, Bothway, T}}((Ones(modeNumber), Fill(λ, modeNumber), m_vec, β_vec, C_vec, D_vec, Fill(profile, modeNumber), Fill(ReferenceFrame((0,0,0), (0,0,0)), modeNumber)))
end

function findmodes!(fibre::Fibre, λ)
    if λ ∉ keys(fibre.modes)
        push!(fibre.modes, λ => findmodes(fibre.refractive_index_profile, λ))
    end
end

function modes(fibre, λ)
    @argcheck haskey(fibre.modes, λ) ErrorException("Mode for that wavelength not yet calculated. Use `findmodes!(fibre, λ)` to pre calculate the modes")
    fibre.modes[λ]
end

## CicurlarStepIndexMode
struct CircularStepIndexProfile{T, M<:Medium{T}} <: AbstractWaveguideProfile{T}
    radius::T
    na::T
    ncore::M
    CircularStepIndexProfile(::Type{T}, r, na, ncore::M) where {T,M} = new{T,M}(r, na, ncore)
end
CircularStepIndexProfile(r, na, ncore) = CircularStepIndexProfile(Float64, r, na, ncore)

struct CircularStepIndexMode{T,Dir<:AbstractDirection,T2<:RealOrComplex{T}} <: AbstractFieldMode{T,Dir}
    e::T2
    wavelength::T
    m::Int
    β::T
    C::T
    D::T
    profile::CircularStepIndexProfile{T}
    frame::ReferenceFrame{T}
end
function CircularStepIndexMode(::Type{T}, ::Type{Dir}, e::E, wavelength, r, m, β, C, D, frame) where {T,Dir,E}
    T2 = E <: Complex ? Complex{T} : T
    CircularStepIndexMode{T,Dir,T2}(e, wavelength, r, m, β, C, D, frame)
end
CircularStepIndexMode(::Type{Dir}, e, wavelength, r, m, β, C, D, frame) where Dir = CircularStepIndexMode(Float64, Dir, e, wavelength, r, m, β, C, D, frame)

mode_type(::Type{<:CircularStepIndexProfile{T}}) where T = CircularStepIndexMode{T, Bothway, T}

function modecondition(profile::CircularStepIndexProfile, λ, m, β)
    α_1 = α1(profile.na, profile.ncore.n, λ, β) # Doesn't work for dispersiveModes
    α_2 = α2(profile.na, profile.ncore.n, λ, β) # Doens't work for dispersiveModes
    besselj(m-1, profile.radius * α_1) / besselj(m, profile.radius * α_1) + α_2 / α_1 * besselkx(m-1, profile.radius * α_2) / besselkx(m, profile.radius * α_2)
end

function modeconstant(profile::CircularStepIndexProfile{T}, λ::Real, m::Integer, β::Number) where {T<:Real}
    ncore = profile.ncore.n
    α_1 = α1(profile.na, ncore, λ, β)
    α_2 = α2(profile.na, ncore, λ, β)
    f1(r) = besselj(m, r[1] * α_1)^2 * r[1];
    (p1, tmp) = hcubature(f1, SVector(zero(T)), SVector(profile.radius); rtol = 1E-8);

    f2(r) = besselk(m, r[1] * α_2)^2 * r[1];
    (p2, tmp) = hcubature(f2, SVector(profile.radius), SVector(2profile.radius); rtol = 1E-8);

    F = T(2π) * p1[1] + T(2π) * p2[1] * besselj(m, profile.radius * α_1)^2 / besselk(m, profile.radius * α_2)^2

    C = 1 / √(F)
    D = C * besselj(m, profile.radius * α_1) / besselk(m, profile.radius * α_2);
    (C, D)
end

function mode_field(mode::CircularStepIndexMode, coord::Polar)
    α_1 = α1(mode.profile.na, mode.profile.ncore.n, mode.wavelength, mode.β)
    α_2 = α2(mode.profile.na, mode.profile.ncore.n, mode.wavelength, mode.β)

    if coord.r < mode.profile.radius
        mode.C * besselj(mode.m, α_1 * coord.r) * exp(im * mode.m * coord.θ)
    else
        mode.D * besselk(mode.m, α_2 * coord.r) * exp(im * mode.m * coord.θ)
    end
end

function mode_field(mode::CircularStepIndexMode, x, y)
    radial_coord = PolarFromCartesian()(Point2D(x, y))
    mode_field(mode, radial_coord)
end

function mode_coupling(mode::CircularStepIndexMode, field, x, y, dA)
   overlap_integral(field, (x,y) -> mode_field(mode, x, y), x, y, dA) 
end

mode_coupling(mode::AbstractWaveguideMode, field::ScalarSpatialBeam) = mode_coupling(mode, field.e, field.x, field.y, field.dA)

function forward_backward_field(fibre::Fibre{<:Any, <:CircularStepIndexProfile, <:CircularStepIndexMode}, field::ScalarSpatialBeam{T,D,<:Any,M}) where {T,D,M}
    allequal(field.modes.wavelength) || error("not done yet")
    number_modes = length(fibre.modes[field.modes.wavelength[1]])
    modes_e = similar(field.modes.e, Complex{T}, number_modes)
    modes_wavelength = Fill(field.modes.wavelength[1], number_modes)
    modes_m = similar(field.modes.e, Int, number_modes)
    modes_β = similar(field.modes.e, T, number_modes)
    modes_C = similar(field.modes.e, T, number_modes)
    modes_D = similar(field.modes.e, T, number_modes)
    i = 1
    for λi in unique(field.modes.wavelength)
        modes = fibre.modes[λi]
        number_modes_λ = length(modes)
        modes_e[i:i+number_modes_λ-1] .= modes.e
        modes_m[i:i+number_modes_λ-1] .= modes.m
        modes_β[i:i+number_modes_λ-1] .= modes.β
        modes_C[i:i+number_modes_λ-1] .= modes.C
        modes_D[i:i+number_modes_λ-1] .= modes.D
    end
    frame = Fill((D == Forward ? first : last)(fibre.frames), number_modes)
    field_t = Beam(StructVector{CircularStepIndexMode{T,Forward,Complex{T}}}((modes_e, modes_wavelength, modes_m, modes_β, modes_C, modes_D, Fill(fibre.refractive_index_profile, number_modes), frame)))
    field_r = Beam(StructArray{ScalarPointSource{T,Backward,T,M}}((field.modes.x, field.modes.y, Zeros(size(field.modes)), field.modes.wavelength, field.modes.medium, field.modes.frame, field.modes.dA)))
    reverse_if_backward(D, (field_r, field_t))
end

# TODO check_input_field(fibre, field)

function _light_interaction!(back_beam::Beam, forw_beam::Beam, fibre::Fibre, ifield::ScalarSpatialBeam{T,D}) where {T,D}
    (r_field, t_field) = reverse_if_backward(D, (back_beam, forw_beam))

    f(modei) = mode_coupling(modei, ifield.modes.e, ifield.modes.x, ifield.modes.y, ifield.modes.dA)
    t_field.modes.e .= f.(t_field.modes)
    (back_beam, forw_beam)
end
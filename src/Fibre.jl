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
        SM = struct_type(R)
        modes = Dict{T, M}()
        new{T, R, M, SM}(profile, length, frames, media, modes)
    end
end
Fibre(profile, length, media, frames) = Fibre(Float64, profile, length, media, frames)


function check_input_field(fibre::Fibre, beam::MeshedSpatialBeam{<:Any,D}) where {D}
    return true
    (n_fibre, frame_fibre) = D == Forward ? first.((fibre.media, fibre.frames)) : last.((fibre.media, fibre.frames))
    @argcheck beam.medium ≈ n_fibre ArgumentError("The medium of the beam and outside the fibre are not the same.")
    @argcheck frame_fibre ≈ beam.frame ArgumentError("The reference frame of the beam and the fibre are not the same.")
    true
end
check_input_field(fibre::Fibre, beam::MeshedBeam) = false

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
function wavefunction_solutions() end

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
    a = (Ones(modeNumber), Fill(λ, modeNumber), m_vec, β_vec, C_vec, D_vec, Fill(profile, modeNumber), Fill(ReferenceFrame((0,0,0), (0,0,0)), modeNumber))
    StructVector{CircularStepIndexMode{T, Bothway, T, Medium{T,T}}}(a) 
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

struct CircularStepIndexMode{T,Dir<:AbstractDirection,T2<:RealOrComplex{T},M} <: AbstractFieldMode{T,Dir}
    e::T2
    wavelength::T
    m::Int
    β::T
    C::T
    D::T
    profile::CircularStepIndexProfile{T,M}
    frame::ReferenceFrame{T}
    function CircularStepIndexMode(e::T2, wavelength::T, m::Int, β::T, C::T, D::T, profile::CircularStepIndexProfile{T,M}, frame::ReferenceFrame{T}) where {T,Dir,T2<:RealOrComplex{T}, M}
        new{T,Forward,T2,M}(e, wavelength, m, β, C, D, profile, frame)
    end
end
function CircularStepIndexMode(::Type{T}, ::Type{Dir}, e::E, wavelength, m, β, C, D, profile::CircularStepIndexProfile{<:Any, M}, frame) where {T,Dir,E,M}
    T2 = E <: Complex ? Complex{T} : T
    CircularStepIndexMode{T,Dir,T2,M}(e, wavelength, m, β, C, D, profile, frame)
end
CircularStepIndexMode(::Type{Dir}, e, wavelength, m, β, C, D, profile, frame) where Dir = CircularStepIndexMode(Float64, Dir, e, wavelength, m, β, C, D, profile, frame)

mode_type(::Type{<:CircularStepIndexProfile{T}}) where T = CircularStepIndexMode{T, Bothway, T, Medium{T,T}}
struct_type(::Type{<:CircularStepIndexProfile{T}}) where T =  StructVector{CircularStepIndexMode{T, Bothway, T, Medium{T, T}}, @NamedTuple{e::Ones{T, 1, Tuple{Base.OneTo{Int}}}, wavelength::Fill{T, 1, Tuple{Base.OneTo{Int}}}, m::Vector{Int}, β::Vector{T}, C::Vector{T}, D::Vector{T}, profile::Fill{CircularStepIndexProfile{T, Medium{T, T}}, 1, Tuple{Base.OneTo{Int}}}, frame::Fill{ReferenceFrame{T}, 1, Tuple{Base.OneTo{Int}}}}, Int}



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

function mode_field(mode::CircularStepIndexMode, coord::R_θ_λ)
    α_1 = α1(mode.profile.na, mode.profile.ncore.n, mode.wavelength, mode.β)
    α_2 = α2(mode.profile.na, mode.profile.ncore.n, mode.wavelength, mode.β)

    if coord[1] < mode.profile.radius
        mode.C * besselj(mode.m, α_1 * coord[1]) * exp(im * mode.m * coord[2])
    else
        mode.D * besselk(mode.m, α_2 * coord[1]) * exp(im * mode.m * coord[2])
    end
end
mode_field(mode::CircularStepIndexMode, coord::X_Y_λ) = mode_field(mode, R_θ_λ(coord))

function mode_coupling(mode::CircularStepIndexMode, field::MeshedSpatialBeam{T,D,C}) where {T,D,C}
   overlap_integral(field.e, (coord) -> mode_field(mode, C(coord)), field.mesh)
end

function forward_backward_field(fibre::Fibre{<:Any, <:CircularStepIndexProfile, <:CircularStepIndexMode}, field::MeshedBeam{T,D,C}) where {T,D,C}
    isone(size(field.mesh)[3]) || error("not done yet")
    λ = centroid(field.mesh, 1).coords[3]

    number_modes = length(fibre.modes[λ])
    modes_e = similar(field.e, Complex{T}, number_modes)
    modes_wavelength = Fill(λ, number_modes)
    modes_m = similar(field.e, Int, number_modes)
    modes_β = similar(field.e, T, number_modes)
    modes_C = similar(field.e, T, number_modes)
    modes_D = similar(field.e, T, number_modes)
    i = 1
    modes = fibre.modes[λ]
    number_modes_λ = length(modes)
    modes_e .= modes.e
    modes_m .= modes.m
    modes_β .= modes.β
    modes_C .= modes.C
    modes_D .= modes.D
    frame = Fill((D == Forward ? first : last)(fibre.frames), number_modes)
    field_t = Beam(StructVector{CircularStepIndexMode{T,Forward,Complex{T}, Medium{T,T}}}((modes_e, modes_wavelength, modes_m, modes_β, modes_C, modes_D, Fill(fibre.refractive_index_profile, number_modes), frame)))
    field_r = MeshedBeam{T,D,C}(field.mesh, Zeros(T, size(field.mesh)), field.medium, field.frame)
    reverse_if_backward(D, (field_r, field_t))
end

# TODO check_input_field(fibre, field)

function _light_interaction!(back_beam::MeshedSpatialBeam, forw_beam::Beam, fibre::Fibre, ifield::MeshedSpatialBeam{T,Forward,C}) where {T,C}
    f(modei) = mode_coupling(modei, ifield)
    for i in eachindex(forw_beam.modes)
        forw_beam.modes.e[i] = f(forw_beam.modes[i])
    end
    (back_beam, forw_beam)
end

function _light_interaction!(back_beam::Beam, forw_beam::MeshedSpatialBeam, fibre::Fibre, ifield::MeshedSpatialBeam{T,Backward,C}) where {T,C}
    _light_interaction!(forw_beam, back_beam, fibre, ifield)
end

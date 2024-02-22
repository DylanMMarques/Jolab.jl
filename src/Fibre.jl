abstract type AbstractWaveguideProfile{T} end
struct Fibre{T,R,M}
    refractive_index_profile::R
    length::T
    frames::Tuple{ReferenceFrame{T}, ReferenceFrame{T}}
    media::Tuple{Medium{T}, Medium{T}}
    modes::Vector{M}
    function Fibre(::Type{T}, profile::R, length, media, frames) where {T,R}
        M = mode_storage_type(R)
        modes = Vector{M}(undef, 0)
        new{T, R, M}(profile, length, frames, media, modes)
    end
end
Fibre(profile, length, media, frames) = Fibre(Float64, profile, length, media, frames)

export Fibre, CircularStepIndexProfile

struct CircularStepIndexProfile{T, M<:Medium{T}} <: AbstractWaveguideProfile{T}
    r::T
    na::T
    ncore::M
    CircularStepIndexProfile(::Type{T}, r, na, ncore::M) where {T,M} = new{T,M}(r, na, ncore)
end
CircularStepIndexProfile(r, na, ncore) = CircularStepIndexProfile(Float64, r, na, ncore)

mode_storage_type(::Type{<:CircularStepIndexProfile{T}}) where T = StructVector{CircularStepIndexMode{T}, @NamedTuple{λ::Fill{T, 1, Tuple{Base.OneTo{Int}}}, m::Vector{Int}, β::Vector{T}, C::Vector{T}, D::Vector{T}}, Int} 

struct CircularStepIndexMode{T} <: AbstractMode{T}
    λ::T
    m::Int
    β::T
    C::T
    D::T
    CircularStepIndexModes(::Type{T}, r, ncore, na, λ, m, β, C, D) where T = new{T}(r, ncore, na, λ, m, β, C, D)
end

α1(na, ncore, λ, β) = √((2π / λ * ncore)^2 - β^2)
α2(na, ncore, λ, β) = √(β^2 + (2π / λ)^2 * (na^2 - ncore^2))
nclad(ncore, na) = √(ncore^2 - na^2)

function modecondition(profile::CircularStepIndexProfile, λ, m, β)
    α_1 = α1(profile.na, profile.ncore.n, λ, β) # Doesn't work for dispersiveModes
    α_2 = α2(profile.na, profile.ncore.n, λ, β) # Doens't work for dispersiveModes
    besselj(m-1, profile.r * α_1) / besselj(m, profile.r * α_1) + α_2 / α_1 * besselkx(m-1, profile.r * α_2) / besselkx(m, profile.r * α_2)
end

function find_wavefunction_solutions(profile::CircularStepIndexProfile{T}, m, λ) where T
    tol = T(1E-12)
    _nclad = nclad(profile.ncore.n, profile.na)
    if m > 0
        α2_min = (T(1E50) * √(T(2m / π))) ^(T(-1/m)) * 2m / T(exp(1)) / profile.r
        βmin = √(α2_min^2 + (T(2π) * (_nclad + tol) / λ)^2);
    else
        βmin = real(2π / λ * (_nclad + tol))
    end
    βmax = real(2π / λ * (profile.ncore.n));
    
    initial_guess = 100
    initial_diff = 500
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
    map(i -> wavefunction_solutions(profile, i, m), rough_βs)
end

function findmodes(profile::P, _λ) where {P<:AbstractWaveguideProfile{T}} where T
    λ = T(_λ)
    m = 0
    sizeA = 10000
    inc_size_i = 1
    m_A = Vector{Int}(undef, sizeA)
    β_A = Vector{T}(undef, sizeA)
    C_A = Vector{T}(undef, sizeA)
    D_A = Vector{T}(undef, sizeA)
    modeNumber = 0
    while true
        condition(β) = modecondition(profile, λ, m, β)

        β_solutions = find_wavefunction_solutions(profile, m, λ)

        isempty(β_solutions) && break
        number_solutions = length(β_solutions)
        for i in 0:(number_solutions-1)
            if modeNumber > (sizeA * inc_size_i) - 2
                map(i -> resize!(i, sizeA * inc_size_i), (m_A, β_A, C_A, D_A))
                inc_size_i += 1
            end

            βᵢ = β_solutions[i + 1]
            modeNumber += 1;
            m_A[modeNumber] = m
            β_A[modeNumber] = βᵢ
            (C_A[modeNumber], D_A[modeNumber]) = modeconstant(profile, λ, m, βᵢ)

            (m == 0) && continue
            βᵢ = β_solutions[i + 1]
            modeNumber += 1;
            m_A[modeNumber] = -m
            β_A[modeNumber] = βᵢ
            (C_A[modeNumber], D_A[modeNumber]) = modeconstant(profile, λ, -m, βᵢ)
        end
        m += 1
    end
    resize!(m_A, modeNumber)
    resize!(β_A, modeNumber)
    resize!(C_A, modeNumber)
    resize!(D_A, modeNumber)
    StructArray{CircularStepIndexMode{T}}((Fill(λ, modeNumber), m_A, β_A, C_A, D_A))
end

function Jolab.findmodes!(fibre::Fibre, λ)
    push!(fibre.modes, findmodes(fibre.refractive_index_profile, λ))
end

function modeconstant(profile::CircularStepIndexProfile{T}, λ::Real, m::Integer, β::Number) where {T<:Real}
    ncore = profile.ncore.n
    α_1 = α1(profile.na, ncore, λ, β)
    α_2 = α2(profile.na, ncore, λ, β)
    f1(r) = besselj(m, r[1] * α_1)^2 * r[1];
    (p1, tmp) = hcubature(f1, SVector(zero(T)), SVector(profile.r); rtol = 1E-8);

    f2(r) = besselk(m, r[1] * α_2)^2 * r[1];
    (p2, tmp) = hcubature(f2, SVector(profile.r), SVector(2profile.r); rtol = 1E-8);

    F = T(2π) * p1[1] + T(2π) * p2[1] * besselj(m, profile.r * α_1)^2 / besselk(m, profile.r * α_2)^2

    C = 1 / √(F)
    D = C * besselj(m, profile.r * α_1) / besselk(m, profile.r * α_2);
    (C, D)
end
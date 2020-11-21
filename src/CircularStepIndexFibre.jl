mutable struct CircularStepIndexModes{T} <: AbstractModes{T}
    r::T
    ncore::T
    na::T
    λ::T
    m::Vector{Int64}
    β::Vector{Complex{T}}
    C::Vector{Complex{T}}
    D::Vector{Complex{T}}
    CircularStepIndexModes{T}(r, ncore, na, λ, m, β, C, D) where T = new{T}(r, ncore, na, λ, m, β, C, D)
end
CircularStepIndexModes(r, ncore, na, λ, m, β, C, D) = CircularStepIndexModes{Float64}(r, ncore, na, λ, m, β, C, D)

mutable struct CircularStepIndexFibre{T} <: AbstractOpticalComponent{T}
    r::T
    ncore::JolabFunction1D{T,Complex{T}}
    na::T
    n₁::JolabFunction1D{T,Complex{T}}
    ref₁::ReferenceFrame{T}
    n₂::JolabFunction1D{T,Complex{T}}
    ref₂::ReferenceFrame{T}
    modes::Vector{CircularStepIndexModes{T}}
    function CircularStepIndexFibre{T}(r, na, ncore, n₁, ref₁, n₂, length) where T
        modes = Vector{CircularStepIndexModes{T}}(undef, 0)
        ref₂ = ReferenceFrame(ref₁.x + length * sin(ref₁.θ) * cos(ref₁.ϕ), ref₁.y + length * sin(ref₁.θ) * sin(ref₁.ϕ), ref₁.z + length * cos(ref₁.θ), ref₁.θ, ref₁.ϕ)
        new{T}(r, ncore, na, n₁, ref₁, n₂, ref₂, modes)
    end
end
CircularStepIndexFibre(r, na, ncore, n₁, ref₁, n₂, length) = CircularStepIndexFibre{Float64}(r, na, ncore, n₁, ref₁, n₂, length)

function rotatestructure!(fib::CircularStepIndexFibre, ref_ref::ReferenceFrame, θ::Real, ϕ::Real)
    rotatereferenceframe!(ref_ref, fib.ref₁, θ, ϕ)
    rotatereferenceframe!(ref_ref, fib.ref₂, θ, ϕ)
end

α1(na, ncore, λ, β) = √((2π / λ * ncore)^2 - β^2)
α2(na, ncore, λ, β) = √(β^2 + (2π / λ)^2 * (na^2 - ncore^2))
nclad(ncore, na) = √(ncore^2 - na^2)

function circularstepindex_modecondition(r::Real, na::Real, ncore::Number, λ::Real, m::Integer, β::Number)::Real
    α_1 = α1(na, ncore, λ, β) # Doesn't work for dispersiveModes
    α_2 = α2(na, ncore, λ, β) # Doens't work for dispersiveModes
    return real(besseljx(m-1, r * α_1) / besseljx(m, r * α_1) + α_2 / α_1 * besselkx(m-1, r * α_2) / besselkx(m, r * α_2))
end

function findmodes!(fibre::CircularStepIndexFibre{T}, λ::Real) where T<:Real
    dispersiveModes = false # can't calculate disperive modes for now
    m = 0;
    sizeA = 10000;
    inc_sizeA = 10000;
    β_A = Vector{Complex{T}}(undef, sizeA)
    m_A = Vector{Int64}(undef, sizeA)
    C_A = Vector{Complex{T}}(undef, sizeA)
    D_A = Vector{Complex{T}}(undef, sizeA)
    modeNumber = 0
    ncore = real(fibre.ncore(λ));
    tmp_nclad = nclad(ncore, fibre.na)
    iWithoutModes = 0;
    while true
        condition(β) = circularstepindex_modecondition(fibre.r, fibre.na, ncore, λ, m, β)

        if m > 0
            α2_min = (1E50 * √(2m / π)) ^(-1/m) * 2m / exp(1) / fibre.r
            βmin = √(α2_min^2 + (2π * (tmp_nclad + @tol) / λ)^2);
        else
            βmin = real(2π / λ * (tmp_nclad + @tol))
        end
        βmax = real(2π / λ * (ncore));

        roots = find_zeros(condition, βmin, βmax, k = 100)
        #(m > 0) && (rootsminus = find_zeros(conditionminus, βmin, βmax))

        (length(roots) == 0) && break
        numberRoots = length(roots)
        for i in 0:(numberRoots-1)
            if modeNumber > sizeA - 2
                resize!(β_A, sizeA + inc_sizeA);
                resize!(m_A, sizeA + inc_sizeA);
                resize!(C_A, sizeA + inc_sizeA);
                resize!(D_A, sizeA + inc_sizeA);
            end

            βᵢ = roots[numberRoots - i]
            # This removes the zeros due to the assymptotas
            if abs(condition(βᵢ)) < 1
                modeNumber += 1;
                m_A[modeNumber] = m
                β_A[modeNumber] = βᵢ
                (C_A[modeNumber], D_A[modeNumber]) = circularstepindex_modeconstant(fibre, λ, m, βᵢ)
            end

            (m == 0) && continue
            βᵢ = roots[numberRoots - i]
            if abs(condition(βᵢ)) < 1
                modeNumber += 1;
                m_A[modeNumber] = -m
                β_A[modeNumber] = βᵢ
                (C_A[modeNumber], D_A[modeNumber]) = circularstepindex_modeconstant(fibre, λ, -m, βᵢ)
            end
        end
        #@show m modeNumber
        m += 1
    end
    resize!(m_A, modeNumber);
    resize!(β_A, modeNumber);
    resize!(C_A, modeNumber);
    resize!(D_A, modeNumber);

    push!(fibre.modes, CircularStepIndexModes{T}(fibre.r, ncore, fibre.na, λ, m_A, β_A, C_A, D_A))
end

function checkapplicability(fibre::CircularStepIndexFibre, field::FieldSpace{T,1}) where T
    checkinplane(ref1(fibre), field.ref) || error("Field space and fibre reference frame must be inplane")
    checkorientation(ref1(fibre), field.ref) || error("Field space and fibre reference frame must have the same orientation")
    isapprox(fibre.n₁(field.λ), field.n, atol = @tol) || error("Field and Fibre must have the same refractive index")
    return true
end

function checkapplicability(fibre::CircularStepIndexFibre, field::FieldSpace{T,-1}) where {T}
    checkinplane(ref2(fibre), field.ref) || error("Field space and fibre reference frame must be inplane")
    checkorientation(ref2(fibre), field.ref) || error("Field space and fibre reference frame must have the same orientation")
    isapprox(fibre.n₂(field.λ), field.n, atol = @tol) || error("Field and Fibre must have the same refractive index")
    return true
end

function getfields_lr(fibre::CircularStepIndexFibre{T}, fieldi::FieldSpace{T,1,Y}) where {T,Y}
    fieldl = FieldSpace{T,-1,Y}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
    fieldr = FieldSpace{T,1,Y}(deepcopy(fieldi.x_X) .- mean(fieldi.x_X), deepcopy(fieldi.y_Y) .- mean(fieldi.y_Y), 0deepcopy(fieldi.e_SXY), fieldi.λ, fibre.n₂(fieldi.λ), ref2(fibre))
    return (fieldl, fieldr)
end

function getfields_lr(fibre::CircularStepIndexFibre{T}, fieldi::FieldSpace{T,-1,Y}) where {T,Y}
    fieldl = FieldSpace{T,-1,Y}(deepcopy(fieldi.x_X) .- mean(fieldi.x_X), deepcopy(fieldi.y_Y) .- mean(fieldi.y_Y), 0deepcopy(fieldi.e_SXY), fieldi.λ, fibre.n₁(fieldi.λ), ref1(fibre))
    fieldr = FieldSpace{T,1,Y}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
    return (fieldl, fieldr)
end

function get_scatteringmatrixtype(fibre::CircularStepIndexFibre{T}, field::FieldSpace{T,D,Y}) where {T,D,Y}
    (fieldl, fieldr) = getfields_lr(fibre, field)
    (sizeX, sizeY) = (length(field.x_X), length(field.y_Y))
    r12 = nothing
    r21 = nothing
    t12 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeX * sizeY)
    t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeX * sizeY)
    return ScatteringMatrix{T,FieldSpace{T,-1,Y}, FieldSpace{T,1,Y}, Nothing, Matrix{Complex{T}}}(r12,t12,r21,t21,fieldl,fieldr)
end

function circularstepindex_modeconstant(fibre::CircularStepIndexFibre{T}, λ::Real, m::Integer, β::Number) where {T<:Real}
    lim = 2fibre.r;
    ncore = real(fibre.ncore(λ));
    α_1 = α1(fibre.na, ncore, λ, β);
    α_2 = α2(fibre.na, ncore, λ, β);
    f1(r) = besselj(m, r[1] * α_1)^2 * r[1];
    (p1, tmp) = hcubature(f1, SVector(zero(T)), SVector(fibre.r); rtol = 1E-8);

    f2(r) = besselk(m, r[1] * α_2)^2 * r[1];
    (p2, tmp) = hcubature(f2, SVector(fibre.r), SVector(2fibre.r); rtol = 1E-8);

    F = 2π * p1[1] + 2π * p2[1] * besselj(m, fibre.r * α_1)^2 / besselk(m, fibre.r * α_2)^2

    C = 1 / √(real(F))
    D = C * besselj(m, fibre.r * α_1) / besselk(m, fibre.r * α_2);
    return (C, D)
end

@inline ref1(fibre::CircularStepIndexFibre) = fibre.ref₁
@inline ref2(fibre::CircularStepIndexFibre) = fibre.ref₂

function Base.getindex(modes::CircularStepIndexModes{T}, i) where T
    return CircularStepIndexModes{T}(modes.r, modes.ncore, modes.na, modes.λ, [modes.m[i]], [modes.β[i]], [modes.C[i]], [modes.D[i]])
end

Base.lastindex(modes::CircularStepIndexModes{T}) where T = numberofmodes(modes)
numberofmodes(modes::CircularStepIndexModes{T}) where T = length(modes.m)

function getmodes(fibre::CircularStepIndexFibre, λ::Real)
    @inbounds for i in eachindex(fibre.modes)
        isapprox(fibre.modes[i].λ, λ, atol = @tol) && return fibre.modes[i]
    end
    findmodes!(fibre, λ)
    return fibre.modes[end]
end

function coefficient_general(fibre::CircularStepIndexFibre{T}, field::FieldSpace{X,D,Y}) where {T,X,Y,D}
    checkapplicability(fibre, field)
    modes = getmodes(fibre, field.λ)

    fieldi = changereferenceframe(field, dir(field) > 0 ? ref1(fibre) : ref2(fibre))
    scat = get_scatteringmatrixtype(fibre, field)

    sizeM = numberofmodes(modes)
    (sizeX, sizeY) = (length(fieldi.x_X), length(fieldi.y_Y))

    t12 = Matrix{Complex{T}}(undef, sizeM, sizeX * sizeY)
    t23 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeM)
    t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeM)
    t32 = Matrix{Complex{T}}(undef, sizeM, sizeX * sizeY)

    fibrelength = (fibre.ref₂.z - fibre.ref₁.z) / cos(fibre.ref₁.θ)

    @inbounds Threads.@threads for iM in 1:sizeM
        i = 1
        propagationmode_term = exp(im * modes.β[iM] * fibrelength)
        for iY in eachindex(fieldi.y_Y)
            for iX in eachindex(fieldi.x_X)
                t12[iM,i] = t_coupling(modes, fieldi, iX, iY, iM)
                t12[iM,i] *= propagationmode_term
                t23[i,iM] = t_fieldspace(modes, dir(field) > 0 ? scat.fieldr : scat.fieldl, iX, iY, iM)

                t32[iM,i] = t_coupling(modes, dir(field) > 0 ? scat.fieldr : scat.fieldl, iX, iY, iM)
                t32[iM,i] *= propagationmode_term
                t21[i,iM] = t_fieldspace(modes, fieldi, iX, iY, iM)
                i += 1
            end
        end
    end
    if dir(field) > 0
        mul!(scat.t₁₂, t23, t12)
        mul!(scat.t₂₁, t21, t32)
    else
        mul!(scat.t₂₁, t23, t12)
        mul!(scat.t₁₂, t21, t32)
    end

    return scat
end

function t_fieldspace(modes::CircularStepIndexModes{T}, field::FieldSpace{T}, iX, iY, iM) where T
    α_1 = α1(modes.na, modes.ncore, modes.λ, modes.β[iM])
    α_2 = α2(modes.na, modes.ncore, modes.λ, modes.β[iM])

    r = √(field.x_X[iX]^2 + field.y_Y[iY]^2)
    ϕ = atan(field.y_Y[iY], field.x_X[iX])
    if r < modes.r
        return modes.C[iM] * besselj(modes.m[iM], α_1 * r) * exp(im * modes.m[iM] * ϕ)
    else
        return modes.D[iM] * besselk(modes.m[iM], α_2 * r) * exp(im * modes.m[iM] * ϕ)
    end
end

function t_coupling(modes::CircularStepIndexModes{T}, field::FieldSpace{T}, iX, iY, iM) where T
	(xmin, xmax) = integralExtremes(field.x_X, iX)
	(ymin, ymax) = integralExtremes(field.y_Y, iY)

    α_1 = α1(modes.na, modes.ncore, modes.λ, modes.β[iM])
    α_2 = α2(modes.na, modes.ncore, modes.λ, modes.β[iM])

    r = √(field.x_X[iX]^2 + field.y_Y[iY]^2)
    if r < modes.r
        fabs = modes.C[iM] * besselj(modes.m[iM], α_1 * r)
    else
        fabs = modes.D[iM] * besselk(modes.m[iM], α_2 * r)
    end
    function fangle(x, y)::T
        ϕ = atan(y, x)
        return - modes.m[iM] * ϕ
    end
    return fabs * exp(im * fangle(field.x_X[iX], field.y_Y[iY])) * (xmax - xmin) * (ymax - ymin)
    return fabs * integrate_exp_xy_x_y(fangle, xmin, xmax, ymin, ymax)
end

function lightinteraction(fibre::CircularStepIndexFibre, field::FieldSpace{T}) where T
    checkapplicability(fibre, field)
    fieldi_newref = changereferenceframe(field, dir(field) > 0 ? ref1(fibre) : ref2(fibre))
    fibre_length = (fibre.ref₂.z - fibre.ref₁.z) / cos(fibre.ref₁.θ)

	modes = getmodes(fibre, field.λ)
    sizeM = numberofmodes(modes)

	(fieldl, fieldr) = getfields_lr(fibre, field)
    fieldout = dir(field) > 0 ? fieldr : fieldl
    fieldl.e_SXY .= zero(Complex{T})
    fieldr.e_SXY .= zero(Complex{T})

    @inbounds Threads.@threads for iM in 1:sizeM
        modeamplitude = zero(Complex{T})
        for iY in eachindex(fieldi_newref.y_Y)
            for iX in eachindex(fieldi_newref.x_X)
                modeamplitude += fieldi_newref.e_SXY[1,iX,iY] * t_coupling(modes, fieldi_newref, iX, iY, iM)
            end
        end
        modeamplitude *= exp(im * modes.β[iM] * fibre_length)
        for iY in eachindex(fieldout.y_Y)
            for iX in eachindex(fieldout.x_X)
                fieldout.e_SXY[1,iX,iY] += modeamplitude * t_fieldspace(modes, fieldout, iX, iY, iM)
            end
        end
    end
    return (fieldl, fieldr)
end

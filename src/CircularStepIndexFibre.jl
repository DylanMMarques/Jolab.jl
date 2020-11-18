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

function getfields_lr(fibre::CircularStepIndexFibre{T}, field::FieldSpace{T,1}) where {T}
    fieldl = FieldSpace{T,-1,Y}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, field.n, field.ref)
    fieldr = FieldSpace{T,1,Y}(copy(fieldi.x_X) .- xmean, copy(fieldi.y_Y) .- ymean, fieldi.e_SXY, fieldi.λ, fibre.n₂(fieldi.λ), ref2(fibre))
    return (fieldl, fieldr)
end

function getfields_lr(fibre::CircularStepIndexFibre{T}, field::FieldSpace{T,-1}) where {T}
    fieldl = FieldSpace{T,-1,Y}(copy(fieldi.x_X) .- xmean, copy(fieldi.y_Y) .- ymean, fieldi.e_SXY, fieldi.λ, field.n, field.ref)
    fieldr = FieldSpace{T,1,Y}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, fibre.n₂(fieldi.λ), ref2(fibre))
    return (fieldl, fieldr)
end

function get_scatteringmatrixtype(fibre::CircularStepIndexFibre{T}, field::FieldSpace{T,D,Y}) where {T,D,Y}
    (fieldl, fieldr) = getfields_lr(fibre, field)
    (sizeX, sizeY) = (length(field.x_X), length(field.y_Y))
    r12 = nothing
    r21 = nothing
    t12 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeX * sizeY)
    t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeX * sizeY)
    if dir(field) > 0
        fieldl = FieldSpace{T,-1,Y}(copy(fieldi.x_X), copy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, field.n, field.ref)
        fieldr = FieldSpace{T,1,Y}(copy(fieldi.x_X) .- xmean, copy(fieldi.y_Y) .- ymean, deepcopy(fieldi.e_SXY), fieldi.λ, fibre.n₂(fieldi.λ), ref2(fibre))
    else
        fieldl = FieldSpace{T,-1,Y}(copy(fieldi.x_X) .- xmean, copy(fieldi.y_Y) .- ymean, fieldi.e_SXY, fieldi.λ, field.n, field.ref)
        fieldr = FieldSpace{T,1,Y}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, fibre.n₂(fieldi.λ), ref2(fibre))
    end
    return ScatteringMatrix{T,FieldSpace{T,-1,Y}, FieldSpace{T,1,Y}, Nothing, Matrix{Complex{T}}}(r12,t12,r21,t21,fieldl,fiedlr)
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

function modefield(fieldmodes::CircularStepIndexModes{T}, index, x, y)::Complex{T} where T
    α_1 = α1(fieldmodes.modes.na, fibre.ncore(fieldmodes.modes.λ), fieldmodes.modes.λ, fieldmodes.modes.β[index])
    α_2 = α2(fieldmodes.modes.na, fibre.ncore(fieldmodes.modes.λ), fieldmodes.modes.λ, fieldmodes.modes.β[index])
    r_var = √(x^2 + y^2)
    ϕ = atan(y, x)
    if (r_var < fieldmodes.r)
        return fieldmodes.C[index] * besselj(fieldmodes.m[index], α_1 * r_var) * exp(im * fieldmodes.m[index] * ϕ)
    else
        return  fieldmodes.D[index] * besselk(fieldmodes.m[index], α_2 * r_var) * exp(im * fieldmodes.m[index] * ϕ)
    end
end

function fieldmodestospace!(space::FieldSpace{T}, fieldmodes::FieldModes) where T
    space.e_SXY .= zero(Complex{T})
    @inbounds Threads.@threads for iM in eachindex(fieldmodes.modes.β)
        @simd for iY in eachindex(y_Y)
            for iX in eachindex(x_X)
                space.e_SXY[1,iX,iY] += modefield(fibre.modes, iM, space.x_X[iX], space.y_Y[iY])
            end
        end
    end
end

function calculatecoupling!(fieldmodes::FieldModes{T}, space::FieldSpace) where T
    modeshape_XY = Array{Complex{T}}(undef, length(space.x_X), length(space.y_Y))

    @inbounds for iM in eachindex(fieldmodes.modes.β)
        Threads.@threads for iY in eachindex(y_Y)
            for iX in eachindex(x_X)
                modeshape_XY[iX,iY] = conj(modefield(fibre.modes, iM, space.x_X[iX], space.y_Y[iY]))
                modeshape_XY[iX,iY] *= space.E_SXY
            end
        end
        fieldmodes.modesamplitude[iM] = ∫∫(modeshape_XY, field.x_X, field.y_Y)
    end
end
function circularstepindex_calculatecoupling!(emodes_A::AbstractArray{Complex{T}}, modes::CircularStepIndexModes{T}, e_SXY::AbstractArray{<:Number,3}, x_X::AbstractVector{<:Real}, y_Y::AbstractVector{<:Real}) where T

    size(e_SXY, 1) == 1 || error("Wrong sizes")
    size(e_SXY, 2) == length(x_X) || error("Wrong sizes")
    size(e_SXY, 3) == length(y_Y) || error("Wrong sizes")

    modeshape_SXY = Array{Complex{T}}(undef, 1, length(x_X), length(y_Y))
    modeshape_XY = reshape(modeshape_SXY, length(x_X), length(y_Y))

    @inbounds @simd for iMode in eachindex(modes.m)
        circularstepindex_modefield!(modeshape_SXY, modes.r, modes.ncore, modes.na, modes.λ, modes.m[iMode], modes.β[iMode], modes.C[iMode], modes.D[iMode], 1, x_X, y_Y, 0)
        conj!(modeshape_SXY)
        @simd for i in eachindex(e_SXY)
            modeshape_SXY[i] *= e_SXY[i]
        end
        emodes_A[iMode] = ∫∫(modeshape_XY, x_X, y_Y)
    end
end

@inline ref1(fibre::CircularStepIndexFibre) = fibre.ref₁
@inline ref2(fibre::CircularStepIndexFibre) = fibre.ref₂

@inline function circularstepindex_calculatecoupling(modes::CircularStepIndexModes{T}, e_SXY::AbstractArray{<:Number,3}, x_X::AbstractVector{<:Real}, y_Y::AbstractVector{<:Real}) where T
    emodes_A = Vector{Complex{T}}(undef, length(modes.m))
    circularstepindex_calculatecoupling!(emodes_A, modes, e_SXY, x_X, y_Y)
    return emodes_A
end

function calculatefieldspace(fieldmodes::FieldModes{T}, x_X::AbstractVector{<:Real}, y_Y::AbstractArray{<:Real}) where T<:Real
    modesNumber = length(fieldmodes.modes.m);
    e_SXY = zeros(Complex{T}, length(x_X) * length(y_Y));
    tmp_e_SXY = Array{Complex{T}, 3}(undef, 1, length(x_X), length(y_Y));
    @inbounds @simd for i in 1:modesNumber
       circularstepindex_modefield!(tmp_e_SXY, fieldmodes.modes.r, fieldmodes.modes.ncore, fieldmodes.modes.na, fieldmodes.modes.λ, fieldmodes.modes.m[i], fieldmodes.modes.β[i], fieldmodes.modes.C[i], fieldmodes.modes.D[i], fieldmodes.modesamplitude[i], x_X, y_Y, 0)
       @simd for i in eachindex(e_SXY)
           e_SXY[i] += tmp_e_SXY[i]
       end
   end
    e_SXY = reshape(e_SXY, 1, length(x_X), length(y_Y))
end

function lightinteraction_interface(fibre::CircularStepIndexFibre{T}, fieldspace::FieldSpace{T,D}) where {T,D}
    reffib = dir(fieldspace) > 0 ? fibre.ref₁ : fibre.ref₂

    changereferenceframe!(fieldspace, reffib);

    modeArg = -1
    @inbounds for i in eachindex(fibre.modes)
        (abs(fibre.modes[i].λ - fieldspace.λ) < @tol) && (modeArg = i; break)
    end
    if modeArg == -1
        findmodes!(fibre, fieldspace.λ)
        modeArg = length(fibre.modes)
    end
    modes = getmodes(fibre, fieldspace.λ)

    modesamplitude = circularstepindex_calculatecoupling(modes, fieldspace.e_SXY, fieldspace.x_X, fieldspace.y_Y);
    modesfield = FieldModes{T,D}(modesamplitude, modes, fieldspace.ref);
end

function lightinteraction_interface(fibre::CircularStepIndexFibre{A}, fieldmodes::FieldModes{T,D}, x_X::AbstractVector{<:Real}, y_Y::AbstractVector{<:Real}) where {T<:Real,A, D}
    dir(fieldmodes) > 0 ? ref = fibre.ref₂ : ref = fibre.ref₁

    changereferenceframe!(fieldmodes, ref);
    e_SXY = calculatefieldspace(fieldmodes, x_X, y_Y);

    n = dir(fieldmodes) > 0 ? fibre.n₂(fieldmodes.modes.λ) : fibre.n₁(fieldmodes.modes.λ)
    return FieldSpace{T,D}(x_X, y_Y, e_SXY, fieldmodes.modes.λ, n, ref)
end

function Base.getindex(modes::CircularStepIndexModes{T}, i) where T
    return CircularStepIndexModes{T}(modes.r, modes.ncore, modes.na, modes.λ, [modes.m[i]], [modes.β[i]], [modes.C[i]], [modes.D[i]])
end

Base.lastindex(modes::CircularStepIndexModes{T}) where T = numberofmodes(modes)
numberofmodes(modes::CircularStepIndexModes{T}) where T = length(modes.m)

function propagationmatrix(fieldl::L, fieldr::L) where {L <: FieldModes{T,D,M} where {T,D,M <:CircularStepIndexModes}}
	t = im * (fieldr.ref.z - fieldl.ref.z) / cos(refold.θ)
	modeweigth_M .*= exp.(tmp_cte .* β_M)
	propMatrix = Diagonal(Vector{Complex{T}}(undef, length(fieldl.modesamplitude)))
	@inbounds @simd for i in eachindex(fieldl.modesamplitude)
		propMatrix.diag[i] = exp(tmp_cte * fieldl.modes.β_M[i])
	end
end

function getmodes(fibre::CircularStepIndexFibre, λ::Real)
    @inbounds for i in eachindex(fibre.modes)
        isapprox(fibre.modes[i].λ, λ, atol = @tol) && return fibre.modes[i]
    end
    findmodes!(fibre, λ)
    return fibre.modes[end]
end

function coefficient_general(fibre::CircularStepIndexFibre{T}, field::FieldSpace{X,D,Y}) where {T,X,Y,D}
    fieldi = changereferenceframe(field, dir(field) > 0 ? ref1(fibre) : ref2(fibre))
    modes = getmodes(fibre, fieldi.λ)

    sizeM = numberofmodes(modes)
    (sizeX, sizeY) = (length(fieldi.x_X), length(fieldi.y_Y))
    dir(fieldi) > 0 || errorToDo()

    t12 = Matrix{Complex{T}}(undef, sizeM, sizeX * sizeY)
    t23 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeM)

    fibrelength = (fibre.ref₂.z - fibre.ref₁.z) / cos(fibre.ref₁.θ)

    Δx = ΔIntegrationTrap(fieldi.x_X)
    Δy = ΔIntegrationTrap(fieldi.y_Y)
    xmean = sum(fieldi.x_X) / sizeX
    ymean = sum(fieldi.y_Y) / sizeY
    @inbounds Threads.@threads for iM in 1:sizeM
        α_1 = α1(modes.na, modes.ncore, modes.λ, modes.β[iM])
        α_2 = α2(modes.na, modes.ncore, modes.λ, modes.β[iM])
        phaseTerm = exp(im * modes.β[iM] * fibrelength)
        i = 1
        for iY in eachindex(fieldi.y_Y)
            for iX in eachindex(fieldi.x_X)
                r_var = √(fieldi.x_X[iX]^2 + fieldi.y_Y[iY]^2)
                ϕ = atan(fieldi.y_Y[iY], fieldi.x_X[iX])
                if (r_var < modes.r)
                    t12[iM,i] = conj(Δx[iX] * Δy[iY] * modes.C[iM] * besselj(modes.m[iM], α_1 * r_var) * exp(im * modes.m[iM] * ϕ)) * phaseTerm
                else
                    t12[iM,i] = conj(Δx[iX] * Δy[iY] * modes.D[iM] * besselk(modes.m[iM], α_2 * r_var) * exp(im * modes.m[iM] * ϕ)) * phaseTerm
                end
                r_var = √((fieldi.x_X[iX] - xmean)^2 + (fieldi.y_Y[iY] - ymean)^2)
                ϕ = atan(fieldi.y_Y[iY] - ymean, fieldi.x_X[iX] - xmean)
                if (r_var < modes.r)
                    t23[i,iM] = modes.C[iM] * besselj(modes.m[iM], α_1 * r_var) * exp(im * modes.m[iM] * ϕ)
                else
                    t23[i,iM] = modes.D[iM] * besselk(modes.m[iM], α_2 * r_var) * exp(im * modes.m[iM] * ϕ)
                end
                i += 1
            end
        end
    end
    t13 = t23 * t12
    if dir(field) > 0
        fieldl = FieldSpace{T,-1,Y}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, field.n, field.ref)
        fieldr = FieldSpace{T,1,Y}(copy(fieldi.x_X) .- xmean, copy(fieldi.y_Y) .- ymean, fieldi.e_SXY, fieldi.λ, fibre.n₂(fieldi.λ), ref2(fibre))
    else
        fieldl = FieldSpace{T,-1,Y}(copy(fieldi.x_X) .- xmean, copy(fieldi.y_Y) .- ymean, fieldi.e_SXY, fieldi.λ, fibre.n₁(fieldi.λ), ref1(fibre))
        fieldr = FieldSpace{T,1,Y}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, field.ref, field.ref)
    end
    r = UniformScaling(zero(Complex{T}))

    return ScatteringMatrix{T, FieldSpace{T,-1,Y}, FieldSpace{T,1,Y}, typeof(r), Matrix{Complex{T}}}(r, t13, r, t13, fieldl, fieldr)
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
    ϕ = atan(field.y_Y[iY], field.x_X[iX])
    if r < modes.r
        fabs = modes.C[iM] * besselj(modes.m[iM], α_1 * r)
    else
        fabs = modes.D[iM] * besselk(modes.m[iM], α_2 * r)
    end
    function fangle(modes, x, y)::T
        ϕ = atan(y, x)
        return - modes.m[iM] * ϕ
    end
    return fabs * integrate_exp_xy_xy_y(fangle, xmin, xmax, ymin, ymax)
end

function lightinteraction(fibre::CircularStepIndexFibre, space::FieldSpace{T}) where T
    checkapplicability(fibre, space)
    fieldi_newref = changereferenceframe(space, dir(space) > 0 ? ref1(fibre) : ref2(fibre))
    fibre_length = (fibre.ref₂.z - fibre.ref₁.z) / cos(fibre.ref₁.θ)

	modes = getmodes(fibre, field.λ)
    sizeM = numberofmodes(modes)

	(fieldl, fieldr) = getfields_lr(fibre, field)
    fieldl.e_SXY .= zero(Complex{T})
    fieldr.e_SXY .= zero(Complex{T})

    for iM in 1:sizeM
        for iY in eachindex(fieldi.y_Y)
            for iX in eachindex(fieldi.x_X)
                modesamplitude += t_coupling(modes, fieldi_newref, iX, iY, iM)
            end
        end
        modesamplitude *= exp(im * modes.β[iM] * fibre_length)
        for iY in eachindex(fieldi.y_Y)
            for iX in eachindex(fieldi.x_X)
                if dir(space) > 0
                    fieldr.e_SXY[1,iX,iY] += modesamplitude * t_fieldspace(modes, fieldi_newref, iX, iY, iM)
                else
                    fieldl.e_SXY[1,iX,iY] += modesamplitude * t_fieldspace(modes, fieldi_newref, iX, iY, iM)
                end
            end
        end
    end
    return (fieldl, fieldr)
end

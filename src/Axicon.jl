export Axicon 
using KernelAbstractions, Atomix

struct Axicon{T, M1<:Medium{T}, M2<:Medium{T}, S} <: AbstractOpticalElement{T}
    α::T
    axicon_medium::M1
    medium::M1
    frame::ReferenceFrame{T}
    solver::S
    function Axicon(::Type{T}, α::T, axicon_medium::M1, medium::M2, frame::ReferenceFrame{T}; solver::S = nothing) where {T,M1,M2,S}
        new{T,M1,M2,S}(α, axicon_medium, medium, frame, solver)
    end
end
Axicon(α, axicon_medium, medium, frame; kwargs...) = Axicon(Float64, α, axicon_medium, medium, frame; kwargs...)

# Meshed beam is input type because needs checking for radial symmetry
@inline function t(::Type{<:MeshedBeam}, axicon::Axicon, β, coord_in::NSR_NSθ_λ, coord_out::R_θ_λ, area)
    (nsr, nsθ, λ) = coord_in
    (r, θ, tmp) = coord_out
    return im * (2π / λ)^2 * exp(-im * 2π / λ * β * r) * besselj0(2π / λ * nsr * r) * area
end

@inline function t(::Type{<:MeshedBeam}, axicon::Axicon, β, coord_in::R_θ_λ, coord_out::NSR_NSθ_λ, area)
    (nsr, nsθ, λ) = coord_out
    (r, θ, tmp) = coord_in
    return - im / 4π^2 * exp(-im * 2π / λ * β * r) * besselj0(2π / λ * nsr * r) * area
end

function _light_interaction!(field_b::FB, field_f::FF, axicon::O, field_i::F) where {FB<:MeshedBeam{TB, Backward, CB}, FF<:MeshedBeam{TF,Forward,CF}, F<:MeshedBeam{T,D,C}, O<:Union{Fourier, Axicon}} where {T,D,C, TF, CF, TB, CB}
    (field_r, field_t) = reverse_if_backward(D, (field_b, field_f))
    field_r.e .= 0
    
    if O <: Axicon
        β = (axicon.axicon_medium.n − axicon.medium.n) * axicon.α
    end

    kernel! = my_kernel_3!(CPU())
    kernel!(field_t, field_i, axicon, β, ndrange = (length(field_t.e),))
    synchronize(CPU())

    (field_b, field_f)
end

@kernel function my_kernel_3!(field_out, @Const(field_i::F), @Const(axicon), @Const(β)) where F
    ind_out = @index(Global)
    out_c = centroid(field_out.mesh, ind_out)
    @inbounds field_out.e[ind_out] = 0.0
    @inbounds for ind_in in eachindex(field_i.e)
        field_out.e[ind_out] += t(F, axicon, β, centroid(field_i.mesh, ind_in), out_c, volume(field_i.mesh, ind_in) / field_i.mesh.spacing[3]) * field_i.e[ind_in]
    end
end
# Radial version
function forward_backward_field(axicon::Union{Axicon, Fourier}, field_i::F) where F<:MeshedBeam{T,D,C} where {T,D,C<:Union{X_Y_λ, NSX_NSY_λ}}
    (symbol, C_T) = F <: MeshedAngularSpectrum ? ((:x, :y), X_Y_λ) : ((:nsx, :nsy), NSX_NSY_λ)
    
    (haskey(axicon.solver, symbol[1]) && haskey(axicon.solver, symbol[2])) || error("Solver must have cartesian coordinate (field x and y or nsx and nsy)")
    x = getfield(axicon.solver, symbol[1])
    y = getfield(axicon.solver, symbol[2])
    e_t = similar(field_i.e, Complex{T}, (length(x), length(y), size(field_i.e)[3]))
    e_r = Zeros(T, size(field_i.e))

    lengths = (length(x), length(y), size(field_i.mesh)[3]) 
    spacing = (step(x), step(y), field_i.mesh.spacing[3])
    origin = C_T(first(x), first(y), field_i.mesh.origin[3])
    
    mesh = CartesianGrid(lengths, origin, spacing)
    
    (dir_r, dir_t) = reverse_if_backward(D, (Backward, Forward))
    field_r = MeshedBeam{T, dir_r, C}(deepcopy(field_i.mesh), e_r, field_i.medium, field_i.frame)
    field_t = MeshedBeam{T, dir_t, C_T}(mesh, e_t, field_i.medium, field_i.frame)
    reverse_if_backward(D, (field_r, field_t))
end

# Radial version
function forward_backward_field(axicon::Union{Axicon, Fourier}, field_i::F) where F<:MeshedBeam{T,D,C,P} where {T,D,C<:Union{R_θ_λ, NSR_NSθ_λ},P}
    (symbol, C_T) = F <: MeshedAngularSpectrum ? (:r, R_θ_λ) : (:nsr, NSR_NSθ_λ)
    
    haskey(axicon.solver, symbol) || error("Solver must have radial coordinate (field r or nsr)")
    r = getfield(axicon.solver, symbol)
    e_t = similar(field_i.e, Complex{T}, (length(r), size(field_i.e)[2], size(field_i.e)[3]))
    e_r = Zeros(T, size(field_i.e))

    lengths = (length(r), size(field_i.mesh)[2], size(field_i.mesh)[3]) 
    spacing = (step(r), field_i.mesh.spacing[2], field_i.mesh.spacing[3])
    origin = C_T(first(r), field_i.mesh.origin[2], field_i.mesh.origin[3])
    
    mesh = CylindricalGrid(lengths, origin, spacing)
    
    (dir_r, dir_t) = reverse_if_backward(D, (Backward, Forward))
    field_r = MeshedBeam{T, dir_r, C,P}(deepcopy(field_i.mesh), e_r, field_i.medium, field_i.frame)
    field_t = MeshedBeam{T, dir_t, C_T,P}(mesh, e_t, field_i.medium, field_i.frame)
    reverse_if_backward(D, (field_r, field_t))
end

function check_input_field(axicon::Axicon, field_i::MeshedBeam)
    code = zero(UInt64)
    isapprox(axicon.frame, field_i.frame) || (code |=  1 << INVALID_FRAME)
    isapprox(axicon.medium, field_i.medium) || (code |= 1 << INVALID_MEDIUM)
    code
end

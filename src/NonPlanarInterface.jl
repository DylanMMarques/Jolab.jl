struct RoughInterface{T, T2<:RealOrComplex{T}, F} <: AbstractOpticalElement{T}
    Δz::F
    mat::Tuple{Medium{T2}, Medium{T2}}
    frame::ReferenceFrame{T}
    function RoughInterface{T}(media::N, Δz::F, frame::ReferenceFrame) where {T, T2<:RealOrComplex{T}, F, N<:Tuple{Medium{T2}, Medium{T2}}}
        new{T,T2,F}(Δz, media, frame)
    end
end

struct RoughInterfaceRRSolver{F1, F2, F3, Y, X, Z, P} <: AbstractSolver
    field_b::F2
    field_f::F3
    field_i::F1
    r12::Y
    t12::Y
    ir12::X
    sr12::X
    it12::X
    st12::X
    r21::Y
    t21::Y
    ir21::X
    sr21::X
    it21::X
    st21::X
    Δz::Z
    plan_fft::P
    tmp_array::X
    tmp_array_2::X
    function RoughInterfaceRRSolver(field_b::F2, field_f::F3, field_i::F1, r12::Y, t12::Y, ri12::X, ti12::X, rs12::X, ts12::X, r21::Y, t21::Y, ri21::X, ti21::X, rs21::X, ts21::X, z::Z, plan_fft::P, tmp_array, tmp_array_2) where {F1, F2, F3, Y, X, Z, P}
        new{F1, F2, F3, Y, X, Z, P}(field_b, field_f, field_i, r12, t12, ri12, rs12, ti12, ts12, r21, t21, ri21, rs21, ti21, ts21, z, plan_fft, tmp_array, tmp_array_2)
    end

end

function rr_1(n1, nsz1, n2, nsz2, λ)
    k0 = 2π / λ
    ru = reflectioncoefficient_interfaces(n1, nsz1, n2, nsz2)
    tu = transmissioncoefficient_interfaces(n1, nsz1, n2, nsz2)
    r_i = -im * k0 / 2 * (n2^2 - n1^2) * (1 + ru)
    r_s = (1 + ru) / nsz1
    t_i = -im * k0 / 2 * (n2^2 - n1^2) * (1 + ru)
    t_s = (1 - ru) / nsz2
    return (ru, tu, r_i, r_s, t_i, t_s)
end

function rr_coeffiecients(n1, nsx1, nsy1, n2, λ)
    nsr_sq = nsx1^2 + nsy1^2
    nsz1 = √(complex(n1^2 - nsr_sq))
    nsz2 = √(complex(n2^2 - nsr_sq))
    return rr_1(n1, nsz1, n2, nsz2, λ)
end

function RoughInterfaceRRSolver(comp::RoughInterface, field_i::MeshedAngularSpectrum{T,D,P}) where {T,D,P}
    field_b, field_f = forward_backward_field(comp, field_i)
    
    (nsx, nsy, _λ) = get_ranges(field_i.mesh)
    λ = only(_λ)
    dkx = step(nsx) / λ
    dky = step(nsy) / λ
    x = fftshift(fftfreq(length(nsx), 1 / dkx))
    y = fftshift(fftfreq(length(nsy), 1 / dky))
    
    tmp_array = similar(field_i.e, Complex{T})
    tmp_array_2 = similar(field_i.e, Complex{T})

    p_fft = plan_fft(field_i.e, (1,2))
    z_numeric = comp.Δz.(x, y')

    r12, t12, ir12, sr12, it12, st12 = ntuple(i -> similar(field_i.e, Complex{T}), 6)
    r21, t21, ir21, sr21, it21, st21 = ntuple(i -> similar(field_i.e, Complex{T}), 6)

    struc_array = StructArray{NTuple{6, Complex{T}}}((r12, t12, ir12, sr12, it12, st12))
    struc_array .= rr_coeffiecients.(first(comp.mat).n, nsx, nsy', last(comp.mat).n, λ)

    RoughInterfaceRRSolver(field_b, field_f, field_i, r12, t12, ir12, it12, sr12, st12, r21, t21, ir21, it21, sr21, st21, z_numeric, p_fft, tmp_array, tmp_array_2)
end

function light_interaction(comp::RoughInterface, field_i)
    msg_code = check_input_field(comp, field_i)
    msg_code == 0 || throw_error_msg(msg_code)
    fourier = RoughInterfaceRRSolver(comp, field_i)
    _light_interaction!(fourier.field_b, fourier.field_f, fourier, field_i)
end


function _light_interaction!(field_b::MeshedBeam{<:Any,Backward}, field_f::MeshedBeam{<:Any,Forward}, solver::RoughInterfaceRRSolver, field_i::MeshedBeam{T,D,NSX_NSY_λ,P}) where {T,D,P<:PolarizationScalar}
    (field_r, field_t) = reverse_if_backward(D, (field_b, field_f))
    (ru, tu, ris, tis, rss, tss) = if D == Forward 
        (solver.r12, solver.t12, solver.ir12, solver.it12, solver.sr12, solver.st12)
    else    
        (solver.r21, solver.t21, solver.ir21, solver.it21, solver.sr21, solver.st21)
    end

    solver.tmp_array .= field_i.e .* ris
    mul!(solver.tmp_array_2, solver.plan_fft, solver.tmp_array)
    solver.tmp_array_2 .*= solver.Δz
    ldiv!(field_r.e, solver.plan_fft, solver.tmp_array_2)
    field_r.e .*= rss
    field_r.e .+= field_i.e .* ru
    
    solver.tmp_array .= field_i.e .* tis
    mul!(solver.tmp_array_2, solver.plan_fft, solver.tmp_array)
    solver.tmp_array_2 .*= solver.Δz
    ldiv!(field_t.e, solver.plan_fft, solver.tmp_array_2)
    field_t.e .*= tss
    field_t.e .+= field_i.e .* tu

    (field_b, field_f)
end

function RoughInterface(::Type{T}, media::Tuple, Δz, frame) where T
    RoughInterface{T}(media, Δz, frame)
end,
function RoughInterface(::Type{T}, media::AbstractVector, Δz, frame) where T
    length(media) == 2 || error("media should have two elements")
    RoughInterface(Float64, (media[1], media[2]), Δz, frame)
end, 
function RoughInterface(media::Tuple, Δz, frame)
    RoughInterface(Float64, media, Δz, frame)
end,
function RoughInterface(media::AbstractVector, Δz, frame)
    RoughInterface(Float64, (media[1], media[2]), Δz, frame)
end

function check_input_field(comp::RoughInterface, field_i::MeshedAngularSpectrum{T,D}) where {D,T}
    msg_code = zero(UInt64)
    field_i.frame ≈ comp.frame || (msg_code |= 1 << INVALID_FRAME)
    field_i.medium ≈ (D == Forward ? first : last)(comp.mat) || (msg_code |= 1 << INVALID_MEDIUM)
    msg_code
end

function forward_backward_field(comp::RoughInterface, field_i::MeshedAngularSpectrum{T,D,C,P}) where {T,D, C, P}
    frame = comp.frame
    medium_b = first(comp.mat)
    medium_f = last(comp.mat)
    e_b = similar(field_i.e, Complex{T})
    e_f = similar(field_i.e, Complex{T})

    field_b = MeshedBeam{T, Backward, C, P}(field_i.mesh, e_b, medium_b, frame)
    field_t = MeshedBeam{T, Forward, C, P}(field_i.mesh, e_f, medium_f, frame)
    (field_b, field_t)
end




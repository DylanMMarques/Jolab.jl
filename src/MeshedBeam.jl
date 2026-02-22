
struct MeshedBeam{T, D, C, P<:AbstractPolarization, V, E<:AbstractArray{<:RealOrComplex{T},3}, M<:Medium{T,<:RealOrComplex{T}}, T2<:RealOrComplex{T}} <: AbstractField{T, D}
    mesh::V
    e::E
    medium::M
    frame::ReferenceFrame{T}
    function MeshedBeam{T,D,C,P}(mesh::V, e::E, medium::Medium{<:Any, M}, frame::ReferenceFrame) where {T,V,E <: AbstractArray{T2},D,C,P,M} where T2
        @argcheck size(e)[1:2] == size(mesh)[1:2] DimensionMismatch
        vectorial_size = number_components(P)

        @argcheck size(e, 3) in (1, vectorial_size) DimensionMismatch
        M1 = M <: Complex ? Complex{T} : T
        T3 = T2 <: Complex ? Complex{T} : T
        new{T,D,C,P,V,E,Medium{T, M1},T3}(mesh, e, medium, frame)
    end
end

function polarization_system(beam::MeshedBeam{T,D,C,P}) where {T,D,C,P}
    P
end

function MonochromaticAngularSpectrum(::Type{T}, ::Type{D}, nsx::AbstractRange, nsy::AbstractRange, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {T, D}
    mesh = CartesianGrid((length(nsx), length(nsy), 1),
        NSX_NSY_λ(first(nsx), first(nsy), λ - eps_factor * eps(T) / 2), 
        (step(nsx), step(nsy), eps_factor * eps(T)))
    e ./= sqrt(eps_factor * eps(T))
    MeshedBeam{T, D, NSX_NSY_λ, PolarizationScalar}(mesh, reshape(e, size(mesh)[1:2]..., 1), medium, frame)
end
MonochromaticAngularSpectrum(::Type{D}, nsx::AbstractRange, nsy::AbstractRange, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {D} = MonochromaticAngularSpectrum(Float64, D, nsx, nsy, e, λ, medium, frame)

function gaussianbeam_electricfield_angspe(::Type{T}, nsr_squared, ω, λ, n) where T
    norm = T(ω) / (4 * T(√2) * T(π)^T(3/2))
    exp(-T(ω)^2 * T(π)^2 * T(nsr_squared) / 4 / T(λ)^2) * norm
end,
function gaussianbeam_electricfield_angspe(::Type{T}, nsx, nsy, ω, λ, n) where T
gaussianbeam_electricfield_angspe(T, (nsx^2 + nsy^2), ω, λ, n)
end

function MonochromaticAngularSpectrum_gaussian(::Type{T}, ::Type{D}, nsx::AbstractRange, nsy::AbstractRange, ω, λ, medium, frame) where {T,D}
    e = complex.(gaussianbeam_electricfield_angspe.(T, nsx, nsy', ω, λ, medium.n))
    MonochromaticAngularSpectrum(T, D, nsx, nsy, e, λ, medium, frame)
end,
function MonochromaticAngularSpectrum_gaussian(D, nsx, nsy, ω, λ, medium, frame) 
    MonochromaticAngularSpectrum_gaussian(Float64, D, nsx, nsy, ω, λ, medium, frame) 
end

function MonochromaticAngularSpectrumRadialSymmetric(::Type{T}, ::Type{D}, nsr::AbstractRange, e::AbstractArray, λ, medium, frame) where {T,D}
    mesh = CylindricalGrid((length(nsr), 1, 1),
        NSR_NSθ_λ(first(nsr), T(0), λ - eps_factor * eps(T) / 2), 
        (step(nsr), 2π, eps_factor * eps(T)))
    e ./= sqrt(eps_factor * eps(T))
    MeshedBeam{T, D, NSR_NSθ_λ, PolarizationScalar}(mesh, reshape(e, size(mesh)[1:2]..., 1), medium, frame)
end

function MonochromaticAngularSpectrumRadialSymmetric_gaussian(::Type{T}, ::Type{D}, nsr::AbstractRange, ω, λ, medium, frame) where {T,D}
    e = complex.(gaussianbeam_electricfield_angspe.(T, nsr.^2, ω, λ, medium.n))
    MonochromaticAngularSpectrumRadialSymmetric(T, D, nsr, e, λ, medium, frame)
end,
function MonochromaticAngularSpectrumRadialSymmetric_gaussian(D, nsr, ω, λ, medium, frame) 
    MonochromaticAngularSpectrumRadialSymmetric_gaussian(Float64, D, nsr, ω, λ, medium, frame) 
end

function MonochromaticSpatialBeam(::Type{T}, ::Type{D}, x::AbstractVector, y::AbstractVector, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {T, D}
    mesh = CartesianGrid((length(x), length(y), 1),
        X_Y_λ(first(x), first(y), λ - eps_factor * eps(T) / 2),  # The -0.5 is to center the point on the face
        (step(x), step(y), eps_factor * eps(T)))
    e ./= sqrt(eps_factor * eps(T))
    MeshedBeam{T, D, X_Y_λ,PolarizationScalar}(mesh, reshape(e, size(mesh)[1:2]..., 1), medium, frame)
end
MonochromaticSpatialBeam(::Type{D}, x::AbstractRange, y::AbstractRange, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {D} = MonochromaticSpatialBeam(Float64, D, x, y, e, λ, medium, frame)

function MonochromaticSpatialBeamVectorial(::Type{T}, ::Type{D}, x::AbstractVector, y::AbstractVector, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {T, D}
    mesh = CartesianGrid((length(x), length(y), 1),
        X_Y_λ(first(x), first(y), λ - eps_factor * eps(T) / 2),  # The -0.5 is to center the point on the face
        (step(x), step(y), eps_factor * eps(T)))
    e .= e ./ sqrt(eps_factor * eps(T))
    MeshedBeam{T, D, X_Y_λ,PolarizationXYZ}(mesh, e, medium, frame)
end,
function MonochromaticSpatialBeamVectorial(::Type{D}, x::AbstractRange, y::AbstractRange, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {D}
    MonochromaticSpatialBeamVectorial(Float64, D, x, y, e, λ, medium, frame)
end




function gaussianbeam_electricfield_space(::Type{T}, r_squared, ω, λ, n) where T
    norm = 2 * T(√(2 / π)) / T(ω)
    T(exp(T(r_squared) * (-4 / T(ω)^2)) * norm)
end,
function gaussianbeam_electricfield_space(::Type{T}, x, y, ω, λ, n) where T
    gaussianbeam_electricfield_space(T, (x^2 + y^2), ω, λ, n)
end

function MonochromaticSpatialBeam_gaussian(::Type{T}, ::Type{D}, x::AbstractVector, y::AbstractVector, ω, λ, medium, frame) where {T,D}
    e = complex.(gaussianbeam_electricfield_space.(T, x, y', ω, λ, medium.n))
    MonochromaticSpatialBeam(T, D, x, y, e, λ, medium, frame)
end,
function MonochromaticSpatialBeam_gaussian(D, x, y, ω, λ, medium, frame) 
    MonochromaticSpatialBeam_gaussian(Float64, D, x, y, ω, λ, medium, frame) 
end
const eps_factor = 1000
function MonochromaticSpatialBeamRadialSymmetric(::Type{T}, ::Type{D}, r::AbstractRange, e::AbstractArray, λ, medium, frame) where {T,D}
    mesh = CylindricalGrid((length(r), 1, 1),
        R_θ_λ(first(r), T(0), λ - eps_factor * eps(T) / 2), 
        (step(r), 2π, eps_factor * eps(T)))
    e ./= sqrt(eps_factor * eps(T))
    MeshedBeam{T, D, R_θ_λ,PolarizationScalar}(mesh, reshape(e, size(mesh)[1:2]..., 1), medium, frame)
end,
function MonochromaticSpatialBeamRadialSymmetric(::Type{D}, r, e, λ, medium, frame) where D
    MonochromaticSpatialBeamRadialSymmetric(Float64, D, r, e, λ, medium, frame) 
end

function MonochromaticSpatialBeamRadialSymmetric_gaussian(::Type{T}, ::Type{D}, r::AbstractRange, ω, λ, medium, frame) where {T,D}
    e = complex.(gaussianbeam_electricfield_space.(T, r.^2, ω, λ, medium.n))
    MonochromaticSpatialBeamRadialSymmetric(T, D, r, e, λ, medium, frame)
end,
function MonochromaticSpatialBeamRadialSymmetric_gaussian(::Type{D}, r, ω, λ, medium, frame) where D
    MonochromaticSpatialBeamRadialSymmetric_gaussian(Float64, D, r, ω, λ, medium, frame) 
end

function intensity(beam::MeshedBeam{T,D,C}) where {T,D,C}
    norm(ind) = C <: AngularSpectrumCoords ? T(16π^4) / (centroid(beam.mesh, ind)[3])^2 : 1
    f(e, ind) = abs2(e) * volume(beam.mesh, ind) * norm(ind)
    mapreduce(f, +, vec(values_nonzeros(beam.e)), eachindex_nonzeros(beam.e))
end

const AngularSpectrumCoords = Union{NSX_NSY_λ, NSR_NSθ_λ}
const MeshedAngularSpectrum{T,D,C<:AngularSpectrumCoords} = MeshedBeam{T,D,C}

const SpatialCoords = Union{X_Y_λ, R_θ_λ}
const MeshedSpatialBeam{T,D,C<:SpatialCoords} = MeshedBeam{T,D,C}

function MeshedPlaneWaveScalar(::Type{T}, ::Type{D}, nsx, nsy, e::T2, λ, medium, frame) where {T,D,T2}
    mesh = CartesianGrid((1, 1, 1),
        NSX_NSY_λ(nsx - eps_factor * eps(T) / 2, nsy - eps_factor * eps(T) / 2, λ - eps_factor * eps(T) / 2), 
        eps_factor .* (eps(T), eps(T), eps(T))
        )
    TE = T2 <: Complex ? Complex{T} : T
    MeshedBeam{T, D, NSX_NSY_λ, PolarizationScalar}(mesh, (@SArray [TE(e);;;]), medium, frame)
end
MeshedPlaneWaveScalar(::Type{D}, nsx, nsy, e, λ, medium, frame) where D = MeshedPlaneWaveScalar(Float64, D, nsx, nsy, e, λ, medium, frame)


function Base.isapprox(beam1::MeshedBeam{T1,D,C}, beam2::MeshedBeam{T2,D,C}; kwargs...) where {T1, T2, D, C}
    isapprox(beam1.mesh, beam2.mesh; kwargs...) &&
    isapprox(beam1.frame, beam2.frame; kwargs...) &&
    isapprox(beam1.medium, beam2.medium; kwargs...) &&
    isapprox(beam1.e, beam2.e; kwargs...)
end

function translate_referenceframe(beam::MeshedAngularSpectrum{T,D,C}, new_origin::X_Y_Z) where {T,D,C}
    dir = beam.frame.direction # should not matter
    frames = reverse_if_backward(D, (beam.frame, ReferenceFrame(T, new_origin, dir)))
    prop = Propagation(T, frames, beam.medium)
    light_interaction(prop, beam)
end

function nsz_nocomplex(n, nsx, nsy) 
    tmp = (n^2 - nsx^2 - nsy^2)
    signbit(real(tmp)) && throw(ArgumentError("The reference frame of evasnescent waves cannot be rotated."))
    return √tmp
end


# function rotate_referenceframe(pw::PlaneWaveScalar{T,D}, new_angles::Point) where {T,D}
#     is_complex_medium(pw.medium) && throw(ArgumentError("Thereference frame of a angular spectrum defined in a medium with a complex refractive index is not defined."))
    
#     rot_matrix = RotXYZ(pw.frame.direction.x, pw.frame.direction.y, pw.frame.direction.z)

#     nsz_val = nsz_nocomplex(pw.medium.n, pw.nsx, pw.nsy) 
#     (new_nsx, new_nsy, new_nsz) = inv(_RotXYZ(new_angles.x, new_angles.y, new_angles.z)) * (rot_matrix * Point{T}(pw.nsx, pw.nsy, nsz_val))
    
#     PlaneWaveScalar(T, D, new_nsx, new_nsy, pw.e, pw.wavelength, pw.medium, ReferenceFrame(pw.frame.origin, new_angles))
# end,
# rotate_referenceframe(pw::Union{AbstractFieldMode{T}, Beam{T}}, new_angles) where T = rotate_referenceframe(pw, convert(Point{T}, new_angles))
translate_referenceframe(pw::Union{AbstractFieldMode{T}, MeshedBeam{T}}, new_origin) where T = translate_referenceframe(pw, X_Y_Z(new_origin))

## Light light_interaction

function light_interaction!(field_b, field_f, comp, beam)
    @argcheck check_input_field(comp, beam) ArgumentError
    @argcheck check_output_fields(field_b, field_f, comp, beam) ArgumentError
    _light_interaction!(field_b, field_f, comp, beam)
end

function light_interaction(comp, beam)
    msg_code = check_input_field(comp, beam)
    msg_code == 0 || throw_error_msg(msg_code)
    (field_b, field_f) = forward_backward_field(comp, beam)
    _light_interaction!(field_b, field_f, comp, beam)
end


function Base.view(beam::MeshedBeam{T,D,C,P}, x, y, z) where {T,D,C,P}
    mesh = @view beam.mesh[_int_to_unitrange(x), _int_to_unitrange(y), _int_to_unitrange(z)]
    e = @view beam.e[_int_to_unitrange(x), _int_to_unitrange(y), _int_to_unitrange(z)]
    return MeshedBeam{T,D,C,P}(mesh, e, beam.medium, beam.frame)
end

_int_to_unitrange(i::Int) = UnitRange(i, i)
_int_to_unitrange(i) = i


function change_polarization_basis(field::MeshedBeam{T,D,C,PolarizationSP}, ::Type{PolarizationXYZ}) where {T,D,C<:AngularSpectrumCoords}
    field_i_s = view(field.e, :, :, 1)
    field_i_p = view(field.e, :, :, 2)

    field_new = similar(field.e, Complex{T}, (size(field.e, 1), size(field.e, 2), 3))
    field_new_x = view(field_new, :, :, 1)
    field_new_y = view(field_new, :, :, 2)
    field_new_z = view(field_new, :, :, 3)

    map(eachindex(field.mesh)) do ind
        (nsr, ϕ, λ) = NSR_NSθ_λ(centroid(field.mesh, ind))
        θ = acos(√(complex(1 - (nsr / field.medium.n)^2)))

        n_s = unit_vector(R_θ_ϕ(zero(T), θ, ϕ), Val(:ϕ))
        n_p = unit_vector(R_θ_ϕ(zero(T), θ, ϕ), Val(:θ))

        comp_s = field_i_s[ind]
        comp_p = field_i_p[ind]

        xyz_field = comp_s * n_s + comp_p * n_p

        field_new_x[ind] = xyz_field[1]
        field_new_y[ind] = xyz_field[2]
        field_new_z[ind] = xyz_field[3]
    end

    return MeshedBeam{T,D,C,PolarizationXYZ}(field.mesh, field_new, field.medium, field.frame)
end

function _unchecked_add!(field1::MeshedBeam{T,D,C,P}, field2::MeshedBeam{T,D,C,P}) where {T,D,C,P}
    field1.e .+= field2.e
    field1
end

function fill_zeros!(field::MeshedBeam{T,D,C,P}) where {T,D,C,P}
    fill!(field.e, zero(eltype(field.e)))
    field
end


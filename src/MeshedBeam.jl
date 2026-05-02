struct MeshedBeam{T, D, C, P<:AbstractPolarization, V, E<:AbstractArray{<:RealOrComplex{T},4}, M<:Medium{T,<:RealOrComplex{T}}, T2<:RealOrComplex{T}} <: AbstractField{T, D}
    mesh::V
    e::E
    medium::M
    frame::ReferenceFrame{T}
    function MeshedBeam{T,D,C,P}(mesh::V, e::E, medium::Medium{<:Any, M}, frame::ReferenceFrame) where {T,V,E <: AbstractArray{T2},D,C,P,M} where T2
        @argcheck size(e)[1:3] == size(mesh)[1:3] DimensionMismatch
        @argcheck size(e, 4) == number_components(P) DimensionMismatch

        M1 = M <: Complex ? Complex{T} : T
        T3 = T2 <: Complex ? Complex{T} : T
        new{T,D,C,P,V,E,Medium{T, M1},T3}(mesh, e, medium, frame)
    end
end

const AngularSpectrumCoords = Union{NSX_NSY_λ, NSR_NSθ_λ}
const MeshedAngularSpectrum{T,D,C<:AngularSpectrumCoords} = MeshedBeam{T,D,C}

const SpatialCoords = Union{X_Y_λ, R_θ_λ}
const MeshedSpatialBeam{T,D,C<:SpatialCoords} = MeshedBeam{T,D,C}

struct Gaussian{T}
    ω::T
end
Base.convert(::Type{Gaussian{T2}}, g::Gaussian{T1}) where {T1, T2} = Gaussian(T2(g))
_electric_field_f(gauss::Gaussian, coord::R_θ_λ) = exp(-(coord[1] * 2 / gauss.ω)^2) * 2 * √(2 / π) / gauss.ω
_electric_field_f(gauss::Gaussian, coord::X_Y_λ) = _electric_field_f(gauss, R_θ_λ(coord))

_electric_field_f(gauss::Gaussian, coord::NSR_NSθ_λ) = exp(-gauss.ω^2 * π^2 * coord[1]^2 / 4 / coord[3]^2) * gauss.ω / (4 * √2 * π^(3/2))
_electric_field_f(gauss::Gaussian, coord::NSX_NSY_λ) = _electric_field_f(gauss, NSR_NSθ_λ(coord))
function electric_field!(e::AbstractArray, gauss::Gaussian{T}, mesh::Domain) where {T}
    e .= map(i -> _electric_field_f(gauss, centroid(mesh, i)), eachindex(mesh))
end



function polarization_system(beam::MeshedBeam{T,D,C,P}) where {T,D,C,P}
    P
end

function _SpatialBeam(::Type{T}, ::Type{D}, grid::Domain{C,3}, e::AbstractArray, medium::Medium, frame::ReferenceFrame) where {C<:SpatialCoords, D<:AbstractDirection, T<:Real}
    e_reshaped = reshape(e, size(grid)[1:3]..., 1)
    MeshedBeam{T, D, C, PolarizationScalar}(grid, e_reshaped, medium, frame)
end,
function _SpatialBeam(::Type{T}, ::Type{D}, grid::Domain{C,3}, g::Gaussian, medium::Medium, frame::ReferenceFrame) where {C<:SpatialCoords, D<:AbstractDirection, T<:Real}
    e = zeros(Complex{T}, size(grid)[1:3]..., 1)
    electric_field!(vec(view(e, :,:,:,1)), g, grid)
    MeshedBeam{T, D, C, PolarizationScalar}(grid, e, medium, frame)
end
function SpatialBeam(args...)
    _SpatialBeam(Float64, args...)
end,
function SpatialBeam(::Type{T}, args...) where T<:Real
    _SpatialBeam(T, args...)
end

# Constructor that accepts a Domain (CartesianGrid or CylindricalGrid) directly for AngularSpectrum
function _AngularSpectrum(::Type{T}, ::Type{D}, grid::Domain{C,3}, e::AbstractArray, medium::Medium, frame::ReferenceFrame) where {C<:AngularSpectrumCoords, D<:AbstractDirection, T<:Real}
    e_reshaped = reshape(e, size(grid)[1:3]..., 1)
    MeshedBeam{T, D, C, PolarizationScalar}(grid, e_reshaped, medium, frame)
end,
function _AngularSpectrum(::Type{T}, ::Type{D}, grid::Domain{C,3}, g::Gaussian, medium::Medium, frame::ReferenceFrame) where {C<:AngularSpectrumCoords, D<:AbstractDirection, T<:Real}
    e = zeros(Complex{T}, size(grid)[1:3]..., 1)
    electric_field!(vec(view(e, :,:,:,1)), g, grid)
    MeshedBeam{T, D, C, PolarizationScalar}(grid, e, medium, frame)
end
function AngularSpectrum(args...)
    _AngularSpectrum(Float64, args...)
end,
function AngularSpectrum(::Type{T}, args...) where T<:Real
    _AngularSpectrum(T, args...)
end


function intensity(beam::MeshedBeam{T,D,C}) where {T,D,C}
    norm(ind) = C <: AngularSpectrumCoords ? T(16π^4) / (centroid(beam.mesh, ind)[3])^2 : 1
    f(e, ind) = abs2(e) * volume(beam.mesh, ind) * norm(ind)
    mapreduce(f, +, vec(values_nonzeros(beam.e)), eachindex_nonzeros(beam.e))
end


function Base.isapprox(beam1::MeshedBeam{T1,D,C1}, beam2::MeshedBeam{T2,D,C2}; kwargs...) where {T1, T2, D, C1, C2}
    !same_coordinate_type(C1, C2) &&
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


function Base.view(beam::MeshedBeam{T,D,C,P}, x, y, z, p) where {T,D,C,P}
    mesh = @view beam.mesh[_int_to_unitrange(x), _int_to_unitrange(y), _int_to_unitrange(z)]
    e = @view beam.e[_int_to_unitrange(x), _int_to_unitrange(y), _int_to_unitrange(z), _int_to_unitrange(p)]
    return MeshedBeam{T,D,C,P}(mesh, e, beam.medium, beam.frame)
end

_int_to_unitrange(i::Int) = UnitRange(i, i)
_int_to_unitrange(i) = i


function change_polarization_basis(field::MeshedBeam{T,D,C,PolarizationSP}, ::Type{PolarizationXYZ}) where {T,D,C<:AngularSpectrumCoords}
    field_i_s = view(field.e, :, :, :, 1)
    field_i_p = view(field.e, :, :, :, 2)

    field_new = similar(field.e, Complex{T}, (size(field.e)[1:3]..., 3))
    field_new_x = view(field_new, :, :, :, 1)
    field_new_y = view(field_new, :, :, :, 2)
    field_new_z = view(field_new, :, :, :, 3)

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
    fill!(field.e, 0)
    field
end


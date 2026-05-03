struct GaussianMode{T,D,T2<:RealOrComplex{T}} <: AbstractFieldMode{T,D}
    e::T2
    sigma::T
    wavelength::T
    frame::ReferenceFrame{T}
end
function GaussianMode(::Type{T}, ::Type{D}, e::T2, sigma, λ, frame) where {T,D,T2}
    T3 = T2 <: Complex ? Complex{T} : T
    GaussianMode{T,D,T3}(e, sigma, λ, frame) 
end

struct SingleModeProfile{T} <: AbstractWaveguideProfile{T}
    mode_field_diameter::T
    SingleModeProfile(::Type{T}, mfd) where T = new{T}(mfd)
end

const SingleModeFibre{T,M,VM} = Fibre{T, SingleModeProfile{T},M,VM}
export SingleModeFibre, SingleModeProfile

mode_type(::Type{SingleModeProfile{T}}) where T = GaussianMode{T, Bothway, Complex{T}}
struct_type(::Type{SingleModeProfile{T}}) where T =  SVector{1, GaussianMode{T, Bothway, Complex{T}}}

function SingleModeFibre(::Type{T}, mfd, length, media, frames) where T
    profile = SingleModeProfile(T, mfd)
    Fibre(T, profile, length, media, frames)
end
SingleModeFibre(mfd, length, media, frames) = SingleModeFibre(Float64, mfd, length, media, frames)

@inline function mode_field(mode::GaussianMode{T}, coord::SpatialCoords) where T 
    gauss = Gaussian(T(mode.sigma))
    _electric_field_f(gauss, coord)
end
@inline function mode_field(mode::GaussianMode{T}, coord::AngularSpectrumCoords) where T
    gauss = Gaussian(T(mode.sigma))
    (2π / coord[3])^2 * 4π^2 * _electric_field_f(gauss, coord)
end

function findmodes(profile::SingleModeProfile{T}, _λ) where T
    return @SVector [GaussianMode(T, Bothway, complex(1), profile.mode_field_diameter, _λ, ReferenceFrame((0,0,0), (0,0,0)))]
end

function check_input_field(fib::Fibre{<:Any, <:SingleModeProfile}, field::MeshedBeam{<:Any, D, <:Any, P}) where {D, P<:AbstractPolarization}
    (frame, medium) = D == Forward ? (fib.frames[1], fib.media[1]) : (fib.frames[2], fib.media[2])
    code = zero(UInt64)
    (P == PolarizationScalar) || (code |= 1 << INVALID_BEAM_TYPE)
    frame ≈ field.frame || (code |= 1 << INVALID_FRAME)
    medium ≈ field.medium || (code |= 1 << INVALID_MEDIUM)
    code
end

function forward_backward_field(fibre::SingleModeFibre, field::MeshedBeam{T, D, C, P}) where {D,T,C,P}
    wavelength = only(get_ranges(field.mesh)[3])
    _modes = modes(fibre, wavelength)
    modes_e = @MVector zeros(Complex{T}, 1)
    modes_wavelength = @SVector [wavelength]
    modes_sigma = @SVector [_modes[1].sigma]

    frames = @SVector [field.frame]
    mode_vec = StructVector{GaussianMode{T, D, Complex{T}}}((modes_e, modes_sigma, modes_wavelength, frames))
    field_t = Beam(mode_vec)
    field_r = MeshedBeam{T,!D,C,P}(field.mesh, Zeros(T, size(field.mesh)..., number_components(P)), field.medium, field.frame)
    reverse_if_backward(D, (field_r, field_t))
end


function get_frame_medium(fibre::SingleModeFibre, ::Type{D}) where D
    if D == Forward
        first.((fibre.frames, fibre.media))
    else
        last.((fibre.frames, fibre.media))
    end
end
function MonochromaticAngularSpectrum(::Type{T}, ::Type{D}, fibre::SingleModeFibre, nsx, nsy, λ) where {T,D}
    mode = modes(fibre, λ)[1]
    frame, medium = get_frame_medium(fibre, D)
    MonochromaticAngularSpectrum_gaussian(T, D, nsx, nsy, mode.sigma, λ, medium, frame)
end
function MonochromaticAngularSpectrumRadialSymmetric(::Type{T}, ::Type{D}, fibre::SingleModeFibre, nsr, λ) where {T,D}
    mode = modes(fibre, λ)[1]
    frame, medium = get_frame_medium(fibre, D)
    MonochromaticAngularSpectrumRadialSymmetric_gaussian(T, D, nsr, mode.sigma, λ, medium, frame)
end
function MonochromaticSpatialBeam(::Type{T}, ::Type{D}, fibre::SingleModeFibre, x, y, λ) where {T,D}
    mode = modes(fibre, λ)[1]
    frame, medium = get_frame_medium(fibre, D)
    MonochromaticSpatialBeam_gaussian(T, D, x, y, mode.sigma, λ, medium, frame)
end
function MonochromaticSpatialBeamRadialSymmetric(::Type{T}, ::Type{D}, fibre::SingleModeFibre, r, λ) where {T,D}
    mode = modes(fibre, λ)[1]
    frame, medium = get_frame_medium(fibre, D)
    MonochromaticSpatialBeamRadialSymmetric_gaussian(T, D, r, mode.sigma, λ, medium, frame)
end

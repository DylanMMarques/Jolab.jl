
struct GaussianMode{T,D,T2<:RealOrComplex{T}} <: AbstractFieldMode{T,D}
    e::T2
    sigma::T
    wavelength::T
    frame::ReferenceFrame{T}
    function GaussianMode(::Type{T}, ::Type{D}, e::T2, sigma, λ, frame) where {T,D,T2}
        T3 = T2 <: Complex ? Complex{T} : T
        new{T,D,T3}(e, sigma, λ, frame) 
    end
end

struct SingleModeProfile{T} <: AbstractWaveguideProfile{T}
    mode_field_diameter::T
    SingleModeProfile(::Type{T}, mfd) where T = new{T}(mfd)
end

const SingleModeFibre{T,M,VM} = Fibre{T, SingleModeProfile{T},M,VM}
export SingleModeFibre, SingleModeProfile

mode_type(::Type{SingleModeProfile{T}}) where T = GaussianMode{T, Bothway, T}
struct_type(::Type{SingleModeProfile{T}}) where T =  SVector{1, GaussianMode{T, Bothway, T}}

function SingleModeFibre(::Type{T}, mfd, length, media, frames) where T
    profile = SingleModeProfile(T, mfd)
    Fibre(T, profile, length, media, frames)
end
SingleModeFibre(mfd, length, media, frames) = SingleModeFibre(Float64, mfd, length, media, frames)

@inline mode_field(mode::GaussianMode{T}, coord::R_θ_λ) where T = gaussianbeam_electricfield_space(T, coord.coords[1]^2, mode.sigma, mode.wavelength, 1)
@inline mode_field(mode::GaussianMode, coord::X_Y_λ) = mode_field(mode, R_θ_λ(coord))

@inline mode_field(mode::GaussianMode{T}, coord::NSR_NSθ_λ) where T = (2π / coord.coords[3])^2 * 4π^2 * gaussianbeam_electricfield_angspe(T, coord.coords[1]^2, mode.sigma, mode.wavelength, 1)
@inline mode_field(mode::GaussianMode, coord::NSX_NSY_λ) = mode_field(mode, NSR_NSθ_λ(coord))

function findmodes(profile::SingleModeProfile{T}, _λ) where T
    return @SVector [GaussianMode(T, Bothway, 1, profile.mode_field_diameter, _λ, ReferenceFrame((0,0,0), (0,0,0)))]
end

function check_input_field(fib::Fibre{<:Any, <:SingleModeProfile}, field::MeshedAngularSpectrum{<:Any, D}) where {D}
    (frame, medium) = D == Forward ? (fib.frames[1], fib.media[1]) : (fib.frames[2], fib.media[2])
    code = zero(UInt64)
    frame ≈ field.frame || (code |= 1 << INVALID_FRAME)
    medium ≈ field.medium || (code |= 1 << INVALID_MEDIUM)Fibre{<:Any, <:SingleModeProfile}
    code
end

function forward_backward_field(fibre::SingleModeFibre, field::MeshedBeam{T, D, C}) where {D,T,C}
    wavelength = get_ranges(field.mesh)[3]
    _modes = map(i -> modes(fibre, i), wavelength)
    number_wavelengths = size(field.mesh, 3)
    modes_e = @MVector zeros(Complex{T}, number_wavelengths)
    modes_wavelength = get_ranges(field.mesh)[3]
    modes_sigma = map(i -> i[1].sigma, _modes)# _modes.mode_field_diameter

    frames = Fill(field.frame, number_wavelengths)

    field_t = Beam(StructVector{GaussianMode{T, D, Complex{T}}}((modes_e, modes_sigma, modes_wavelength, frames)))
    field_r = MeshedBeam{T,!D,C}(field.mesh, Zeros(T, size(field.mesh)), field.medium, field.frame)
    reverse_if_backward(D, (field_r, field_t))
end


function get_frame_medium(fibre::SingleModeFibre, ::Type{D}) where D
    map(i -> D == Forward ? first(i) : last(i), (fibre.media, fibre.frames))
end
function MonochromaticAngularSpectrum(::Type{T}, ::Type{D}, fibre::SingleModeFibre, nsx, nsy, λ) where {T,D}
    mode = modes(fibre, λ)[1]
    MonochromaticAngularSpectrum_gaussian(T, D, nsx, nsy, mode.sigma, λ, get_frame_medium(fibre, D)...)
end
function MonochromaticAngularSpectrumRadialSymmetric(::Type{T}, ::Type{D}, fibre::SingleModeFibre, nsr, λ) where {T,D}
    mode = modes(fibre, λ)[1]
    MonochromaticAngularSpectrumRadialSymmetric_gaussian(T, D, nsr, mode.sigma, λ, get_frame_medium(fibre, D)...)
end
function MonochromaticSpatialBeam(::Type{T}, ::Type{D}, fibre::SingleModeFibre, x, y, λ) where {T,D}
    mode = modes(fibre, λ)[1]
    MonochromaticSpatialBeam_gaussian(T, D, x, y, mode.sigma, λ, get_frame_medium(fibre, D)...)
end
function MonochromaticSpatialBeamRadialSymmetric(::Type{T}, ::Type{D}, fibre::SingleModeFibre, r, λ) where {T,D}
    mode = modes(fibre, λ)[1]
    MonochromaticSpatialBeamRadialSymmetric_gaussian(T, D, r, mode.sigma, λ, get_frame_medium(fibre, D)...)
end
struct Beam{T, D, M<:StructArray{<:AbstractFieldMode{T,D}}} <: AbstractField{T,D}
    modes::M
    function Beam(modes::M) where {M<:StructArray{<:AbstractFieldMode{T,D}}} where {T,D}
        new{T, D, M}(modes)
    end
end

const ScalarAngularSpectrumBeam{T,D,E,M,N,V,C} = Beam{T, D, StructArray{PlaneWaveScalar{T,D,E,M},N,V,C}} 
const VectorialAngularSpectrumBeam{T,D,E,M,N,V,C} = Beam{T, D, StructArray{PlaneWaveVectorial{T,D,E,M},N,V,C}} 
const AngularSpectrumBeam{T,D,E,M,N,V,C} = Union{ScalarAngularSpectrumBeam{T,D,E,M,N,V,C}, VectorialAngularSpectrumBeam{T,D,E,M,N,V,C}}

const ScalarSpatialBeam{T,D,E,M,N,V,C} = Beam{T, D, StructArray{ScalarPointSource{T,D,E,M},N,V,C}} 
const VectorialSpatialBeam{T,D,E,M,N,V,C} = Beam{T, D, StructArray{VectorialPointSource{T,D,E,M},N,V,C}} 
const SpatialBeam{T,D,E,M,N,V,C} = Union{ScalarSpatialBeam{T,D,E,M,N,V,C}, VectorialSpatialBeam{T,D,E,M,N,V,C}}


function monochromatic_angularspectrum(::Type{T}, ::Type{D}, nsx, nsy, e::AbstractArray{E}, wavelength, medium::M, frame::ReferenceFrame{T}) where {T, M<:Medium, E,D}
    E2 = E <: Complex ? Complex{T} : T

    number_modes = length(e)
    @argcheck all(i -> length(i) == number_modes, (nsx, nsy)) DimensionMismatch

    modes = StructArray{PlaneWaveScalar{T,D,E2,M}}((vec(nsx), vec(nsy), vec(e), Fill(wavelength, number_modes), Fill(medium, number_modes), Fill(frame, number_modes), Fill((nsx[2] - nsx[1]) * (nsy[1,2] - nsy[1]), number_modes)))
    Beam(modes)
end
monochromatic_angularspectrum(::Type{D}, nsx, nsy, e, wavelength, medium, frame) where D = monochromatic_angularspectrum(Float64, D, nsx, nsy, e, wavelength, medium, frame)

function monochromatic_fieldspace(::Type{T}, ::Type{D}, x, y, e::AbstractArray{E}, wavelength, medium::M, frame::ReferenceFrame{T}) where {T, M<:Medium, E,D}
    E2 = E <: Complex ? Complex{T} : T

    number_modes = size(e)
    @argcheck all(i -> size(i) == number_modes, (x, y)) DimensionMismatch
    modes = StructArray{ScalarPointSource{T,D,E2,M}}((x, y, e, Fill(wavelength, number_modes), Fill(medium, number_modes), Fill(frame, number_modes), Fill((x[2] - x[1]) * (y[1,2] - y[1]), number_modes)))
    Beam(modes)
end
monochromatic_fieldspace(::Type{D}, x, y, e, wavelength, medium, frame) where D = monochromatic_fieldspace(Float64, D, x, y, e, wavelength, medium, frame)

function intensity(beam::ScalarAngularSpectrumBeam) 
    f(e, dA, mat) = abs2(e) * (dA * mat.n)
    mapreduce(f, +, beam.modes.e, beam.modes.dA, beam.modes.medium)
end


## Rotation and translate_referenceframe

function translate_referenceframe(angspe::Beam{T}, new_origin::Point3D) where {T}
    Δpos_f(origin, ref2) = origin - ref2.origin 
    Δpos = broadcast(i -> Δpos_f(new_origin, i), angspe.modes.frame)
    
    rot_matrix_f(ref) = RotXYZ(ref.direction.x, ref.direction.y, ref.direction.z) 
    rot_matrix = rot_matrix_f.(angspe.modes.frame)
    rΔpos = inv.(rot_matrix) .* Δpos

    nsz_nocomplex(m::Medium, nsx, nsy) = complex(m.n^2 - nsx^2 - nsy^2)^(1/2)

    phase_term(λ, rΔpos, nsx, nsy, medium) = begin 
        nsz = nsz_nocomplex(medium, nsx, nsy)
        exp(im * 2T(π) / λ * dot(rΔpos, (nsx, nsy, nsz)))
    end
    new_e = angspe.modes.e .* phase_term.(angspe.modes.wavelength, rΔpos, angspe.modes.nsx, angspe.modes.nsy, angspe.modes.medium)
    new_nsx = copy(angspe.modes.nsx)
    new_nsy = copy(angspe.modes.nsy)
    ref_new = broadcast(i -> ReferenceFrame(new_origin, i.direction), angspe.modes.frame)
    new_modes = StructArray{eltype(angspe.modes)}((new_nsx, new_nsy, new_e, angspe.modes.wavelength, angspe.modes.medium, ref_new, angspe.modes.dA))
    Beam(new_modes)
end

function rotate_referenceframe(angspe::Beam, refnew::ReferenceFrame)
    error("TO DO")
end


## Plane wave
function translate_referenceframe(pw::PlaneWaveScalar{T,D}, new_origin::Point3D) where {T,D}
    Δpos = new_origin - pw.frame.origin
    
    rot_matrix = RotXYZ(pw.frame.direction.x, pw.frame.direction.y, pw.frame.direction.z)
    rΔpos = inv(rot_matrix) * Δpos

    nsz_val = (pw.medium.n^2 - pw.nsx^2 - pw.nsy^2)^(T(1/2)) 
    e_new = pw.e * exp(im * 2T(π) / pw.wBeamavelength * dot(rΔpos, (pw.nsx, pw.nsy, nsz_val)))

    PlaneWaveScalar(T, D, pw.nsx, pw.nsy, e_new, pw.wavelength, pw.medium, ReferenceFrame(new_origin, pw.frame.direction)) 
end

function rotate_referenceframe(pw::PlaneWaveScalar{T,D}, new_angles::Point3D) where {T,D}
    is_complex_medium(pw.medium) && throw(ArgumentError("Thereference frame of a angular spectrum defined in a medium with a complex refractive index is not defined."))
    
    rot_matrix = RotXYZ(pw.frame.direction.x, pw.frame.direction.y, pw.frame.direction.z)

    nsz_val = nsz_nocomplex(pw.medium.n, pw.nsx, pw.nsy) 
    (new_nsx, new_nsy, new_nsz) = inv(RotXYZ(new_angles.x, new_angles.y, new_angles.z)) * (rot_matrix * Point3D{T}(pw.nsx, pw.nsy, nsz_val))

    PlaneWaveScalar(T, D, new_nsx, new_nsy, pw.e, pw.wavelength, pw.medium, ReferenceFrame(pw.frame.origin, new_angles))
end

rotate_referenceframe(pw::Union{AbstractFieldMode{T}, Beam{T}}, new_angles) where T = rotate_referenceframe(pw, convert(Point3D{T}, new_angles))
translate_referenceframe(pw::Union{AbstractFieldMode{T}, Beam{T}}, new_origin) where T = translate_referenceframe(pw, convert(Point3D{T}, new_origin))
()
number_modes(beam::Beam) = length(beam.modes)


## Light light_interaction

function light_interaction!(field_b, field_f, comp, beam)
    @argcheck check_input_field(comp, beam) ArgumentError
    @argcheck check_output_fields(field_b, field_f, comp, beam) ArgumentError
    _light_interaction!(field_b, field_f, comp, beam)
end

function light_interaction(comp, beam)
    @argcheck check_input_field(comp, beam) ArgumentError
    (field_b, field_f) = forward_backward_field(comp, beam)
    _light_interaction!(field_b, field_f, comp, beam)
end

function light_interaction(comp, beam::PlaneWaveScalar)
    @argcheck check_input_field(comp, beam) ArgumentError
    _light_interaction(comp, beam)
end

function Base.isapprox(a::Beam{T1,D1,<:StructArray{M1}}, b::Beam{T2,D2,<:StructArray{M2}}; kwargs...) where {T1, T2, D1, D2, M1, M2}
    D1 == D2 || return false
    same_mode_type(M1, M2) || return false
    fields = fieldnames(M1)
    for field in fields
        all(isapprox.(component(a.modes, field), component(b.modes, field); kwargs...)) || return false
    end
    return true
end
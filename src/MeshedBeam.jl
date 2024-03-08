coord_list = (:X_Y_λ,
    :X_Y_t,
    :X_NSY_λ,
    :X_NSY_t,
    :NSX_Y_λ,
    :NSX_Y_t,
    :NSX_NSY_λ,
    :NSX_NSY_t,
    :R_θ_λ,
    :R_θ_t,
    :NSR_NSθ_λ,
    :NSR_NSθ_t
)
for coordi in coord_list
    eval(quote
        struct $coordi{T}
            coords::Vec{3, T}
        end
        Base.getindex(c::$coordi, i) = c.coords[i]
        $coordi(coords::Point3) = $coordi(coords.coords)
    end)
end

for (cart, polar) in zip((:NSX_NSY_λ, :NSX_NSY_t, :X_Y_λ, :X_Y_t), (:NSR_NSθ_λ, :NSR_NSθ_t, :R_θ_λ, :R_θ_t))
    eval(quote
        function $cart(coords::$(polar){T}) where T
            car = CartesianFromPolar()(Polar(coords[1], coords[2]))
            $cart{T}((car.x, car.y, coords[3]))
        end
        function $polar(coords::$(cart){T}) where T
            car = PolarFromCartesian()(coords[1:2])
            $polar{T}((car.r, car.θ, coords[3]))
        end
    end)
end


struct MeshedBeam{T, D, C, V, E<:AbstractArray{<:RealOrComplex{T},3}, M<:Medium{T,<:RealOrComplex{T}}, T2<:RealOrComplex{T}} <: AbstractField{T, D}
    mesh::V
    e::E
    medium::M
    frame::ReferenceFrame{T}
    function MeshedBeam{T,D,C}(mesh::V, e::E, medium::Medium{<:Any, M}, frame::ReferenceFrame) where {T,V,E <: AbstractArray{T2},D,C,M} where T2
        @argcheck size(e) == size(mesh) DimensionMismatch
        M1 = M <: Complex ? Complex{T} : T
        T3 = T2 <: Complex ? Complex{T} : T
        new{T,D,C,V,E,Medium{T, M1},T3}(mesh, e, medium, frame)
    end
end

function MonochromaticAngularSpectrum(::Type{T}, ::Type{D}, nsx::AbstractRange, nsy::AbstractRange, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {T, D}
    mesh = CartesianGrid((length(nsx), length(nsy), 1),
        Point3(first(nsx), first(nsy), λ - eps(T) / 2), 
        (step(nsx), step(nsy), eps(T)))
    MeshedBeam{T, D, NSX_NSY_λ}(mesh, reshape(e, size(e)..., 1), medium, frame)
end
MonochromaticAngularSpectrum(::Type{D}, nsx::AbstractRange, nsy::AbstractRange, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {D} = MonochromaticAngularSpectrum(Float64, D, nsx, nsy, e, λ, medium, frame)

function MonochromaticSpatialBeam(::Type{T}, ::Type{D}, x::AbstractVector, y::AbstractVector, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {T, D}
    @argcheck size(e) == (length(x), length(y)) DimensionMismatch
    mesh = CartesianGrid((length(x), length(y), 1),
        Point3(first(x), first(y), λ - eps(T) / 2),  # The -0.5 is to center the point on the face
        (step(x), step(y), eps(T)))
    MeshedBeam{T, D, X_Y_λ}(mesh, reshape(e, size(e)..., 1), medium, frame)
end
MonochromaticSpatialBeam(::Type{D}, x::AbstractRange, y::AbstractRange, e::AbstractArray, λ, medium::Medium, frame::ReferenceFrame) where {D} = MonochromaticSpatialBeam(Float64, D, x, y, e, λ, medium, frame)

function intensity(beam::MeshedBeam)
    f(e, ind) = abs2(e) * (volume(beam.mesh, ind) * beam.medium.n)
    mapreduce(f, +, vec(beam.e), eachindex(beam.e))
end

const AngularSpectrumCoords = Union{NSX_NSY_λ, NSR_NSθ_λ}
const MeshedAngularSpectrum{T,D,C<:AngularSpectrumCoords} = MeshedBeam{T,D,C}

const SpatialCoords = Union{X_Y_λ, R_θ_λ, X_Y_t, R_θ_t}
const MeshedSpatialBeam{T,D,C<:SpatialCoords} = MeshedBeam{T,D,C}

function MeshedPlaneWaveScalar(::Type{T}, ::Type{D}, nsx, nsy, e::T2, λ, medium, frame) where {T,D,T2}
    mesh = CartesianGrid((1, 1, 1),
        Point{3,T}(nsx - eps(T) / 2, nsy - eps(T) / 2, λ - eps(T) / 2), 
        (eps(T), eps(T), eps(T))
        )
    TE = T2 <: Complex ? Complex{T} : T
    MeshedBeam{T, D, NSX_NSY_λ}(mesh, (@SArray [TE(e);;;]), medium, frame)
end
MeshedPlaneWaveScalar(::Type{D}, nsx, nsy, e, λ, medium, frame) where D = MeshedPlaneWaveScalar(Float64, D, nsx, nsy, e, λ, medium, frame)

function Base.isapprox(beam1::MeshedBeam{T1,D,C}, beam2::MeshedBeam{T2,D,C}; kwargs...) where {T1, T2, D, C}
    isapprox(beam1.mesh, beam2.mesh; kwargs...) || return false
    isapprox(beam1.frame, beam2.frame; kwargs...) || return false
    isapprox(beam1.medium, beam2.medium; kwargs...) || return false
    return true
end

function translate_referenceframe(beam::MeshedAngularSpectrum{T,D,C}, new_origin::Point3D) where {T,D,C}
    Δpos = new_origin - beam.frame.origin

    rot = RotXYZ(beam.frame.direction.x, beam.frame.direction.y, beam.frame.direction.z) 
    rΔpos = inv(rot) * Δpos

    function phase_term(medium, rΔpos, coord::NSX_NSY_λ)
        (nsx, nsy, λ) = coord.coords
        nsz = (medium.n^2 - nsx^2 - nsy^2)^(1/2) # for Enzyme
        exp(im * 2T(π) / λ * dot(rΔpos, (nsx, nsy, nsz)))
    end
    phase_term(medium, rΔpos, coord::NSR_NSθ_λ) = phase_term(medium, rΔpos, convert(NSX_NSY_λ, coord))
    phase_term(medium, rΔpos, index) = phase_term(medium, rΔpos, C(centroid(beam.mesh, index).coords))

    new_e = reshape(map(i -> beam.e[i] * phase_term(beam.medium, rΔpos, i), eachindex(beam.e)), size(beam.e))
    MeshedBeam{T,D,C}(deepcopy(beam.mesh), new_e, deepcopy(beam.medium), ReferenceFrame(new_origin, beam.frame.direction))
end

function nsz_nocomplex(n, nsx, nsy) 
    tmp = (n^2 - nsx^2 - nsy^2)
    signbit(real(tmp)) && throw(ArgumentError("The reference frame of evasnescent waves cannot be rotated."))
    return √tmp
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
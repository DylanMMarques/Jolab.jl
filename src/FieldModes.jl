mutable struct FieldModes{T<:Real,D,M<:AbstractModes{T}} <: AbstractFieldMonochromatic{T,D}
    modesamplitude::Vector{Complex{T}}
    modes::M
	ref::ReferenceFrame{T}
    function FieldModes{T}(modesamplitude, modes::M, dir, ref) where {T, M}
        length(modesamplitude) == length(modes.m) || error("size must be the same")
        new{T,dir,M}(modesamplitude, modes, dir, ref)
    end
end

FieldModes(modesamplitude, modes, dir, ref) where T = FieldModes{Float64}(modesamplitude, modes, dir, ref) where T
Base.eltype(modes::AbstractModes{T}) where T = T

function changereferenceframe!(modeweigth_M::AbstractVector{<:Number}, β_M::AbstractVector{<:Number}, dir::Integer, refold::ReferenceFrame, refnew::ReferenceFrame)
	checkorientation(refold,refnew) || error("Rotation on a FieldMode cannot be done")
	checkinline(refold, refnew) && translatereferenceframe!(modeweigth_M, β_M, dir, refold, refnew)
end

function translatereferenceframe!(modeweigth_M::AbstractVector{<:Number}, β_M::AbstractVector{<:Number}, dir::Integer, refold::ReferenceFrame, refnew::ReferenceFrame)
	t = (refnew.z - refold.z) / cos(refold.θ)
	tmp_cte = im * t * dir;
	modeweigth_M .*= exp.(tmp_cte .* β_M)
end

function changereferenceframe!(fieldmodes::FieldModes, ref::ReferenceFrame)
	changereferenceframe!(fieldmodes.modesamplitude, fieldmodes.modes.β, fieldmodes.dir, fieldmodes.ref, ref)
	fieldmodes.ref = ref;
end

function samedefinitions(fieldl::L, fieldr::L) where L <: FieldModes
	(checkorientation(fieldl.ref, fieldr.ref) && checkinline(fieldl.ref, fieldr.ref)) || return false
	fieldl.modes == fieldr.modes || return false
	return true
end

abstract type AbstractModes{T<:Real} end

mutable struct FieldModes{T<:Real,M<:AbstractModes{T}} # <: AbstractMonochromaticField
    modesamplitude::Vector{Complex{T}}
    modes::M
    dir::Int8
	ref::ReferenceFrame{T}
    function FieldModes{T}(modesamplitude, modes::M, dir, ref) where {T, M}
        length(modesamplitude) == length(modes.m) || error("size must be the same")
        new{T,M}(modesamplitude, modes, dir, ref)
    end
end
FieldModes(modesamplitude, modes, dir, ref) where T = FieldModes{Float64}(modesamplitude, modes, dir, ref) where T
Base.eltype(modes::AbstractModes{T}) where T = T

function changereferential!(modeweigth_M::AbstractVector{<:Number}, β_M::AbstractVector{<:Number}, dir::Integer, refold::ReferenceFrame, refnew::ReferenceFrame)
	!checkorientation(refold,refnew) && error("Rotation on a FieldMode cannot be done")
	checkinline(refold, refnew) && translatereferential!(modeweigth_M, β_M, dir, refold, refnew)
end

function translatereferential!(modeweigth_M::AbstractVector{<:Number}, β_M::AbstractVector{<:Number}, dir::Integer, refold::ReferenceFrame, refnew::ReferenceFrame)
	t = (refnew.z - refold.z) / cos(refold.θ)
	tmp_cte = im * t * dir;
	modeweigth_M .*= exp.(tmp_cte .* β_M)
end

function changereferential!(fieldmodes::FieldModes, ref::ReferenceFrame)
	changereferential!(fieldmodes.modesamplitude, fieldmodes.modes.β, fieldmodes.dir, fieldmodes.ref, ref)
	fieldmodes.ref = ref;
end

changereferential(fieldmodes::FieldModes, ref::ReferenceFrame) = changereferential!(deepcopy(fieldmodes), ref);

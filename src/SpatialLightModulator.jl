mutable struct SpatialLightModulator{T} <: AbstractOpticalComponent{T}
    x_X::AbstractVector{T}
    y_Y::AbstractVector{T}
    r_XY::AbstractArray{Complex{T},2}
    t_XY::AbstractArray{Complex{T},2}
    ref::ReferenceFrame{T}
    function SpatialLightModulator(x_X, y_Y, r_XY, t_XY, ref)
        T = promote_type(eltype(x_X), eltype(y_Y), typeof(real(t_XY[1])), typeof(real(r_XY[1])), typeof(ref.x));
        length(x_X) == size(t_XY, 1) || error("Wrong sizes")
        length(y_Y) == size(t_XY, 2) || error("Wrong sizes")
        new{T}(x_X, y_Y, r_XY, t_XY, ref)
    end
end

SpatialLightModulator_slm(x_X, y_Y, t_XY, ref) = return SpatialLightModulator(x_X, y_Y, zeros(eltype(t_XY), size(t_XY)), t_XY, ref)
SpatialLightModulator_mask(x_X, y_Y, t_XY, ref) = return SpatialLightModulator(x_X, y_Y, ones(eltype(t_XY), size(t_XY)) .- t_XY, t_XY, ref)
SpatialLightModulator_reflectivemask(x_X, y_Y, r_XY, ref) = return SpatialLightModulator(x_X, y_Y, r_XY, ones(eltype(r_XY), size(r_XY)) .- r_XY, ref)

function coefficient_general(slm::SpatialLightModulator{T}, field::FieldSpace{T}) where {T<:Real}
	checkorientation(field.ref, slm.ref) || errorToDo()
    fieldi = changereferenceframe(field, slm.ref)

    slm_itp_t = extrapolate(interpolate((slm.x_X, slm.y_Y), slm.t_XY, Gridded(Constant())), one(Complex{T}))
    slm_itp_r = extrapolate(interpolate((slm.x_X, slm.y_Y), slm.r_XY, Gridded(Constant())), zero(Complex{T}))

	(sizeX, sizeY) = (length(fieldi.x_X), length(fieldi.y_Y))
	r = Diagonal(Vector{Complex{T}}(undef, sizeX * sizeY))
	t = Diagonal(Vector{Complex{T}}(undef, sizeX * sizeY))
	i = 1
	@inbounds for iY in eachindex(fieldi.y_Y)
		for iX in eachindex(fieldi.x_X)
			t.diag[i] = slm_itp_t(fieldi.x_X[iX], fieldi.y_Y[iY])
			r.diag[i] = slm_itp_r(fieldi.x_X[iX], fieldi.y_Y[iY])
			i += 1
		end
	end
	fieldl = FieldSpace{T}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, -1, fieldi.ref)
	fieldr = FieldSpace{T}(copy(fieldi.x_X), copy(fieldi.y_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, 1, fieldi.ref)
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), Diagonal{Complex{T},Vector{Complex{T}}}, Diagonal{Complex{T},Vector{Complex{T}}}}(r, t, r, t, fieldl, fieldr)
end

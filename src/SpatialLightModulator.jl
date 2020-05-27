mutable struct SpatialLightModulator{T} <: AbstractOpticalComponent{T}
    x_X::AbstractVector{T}
    y_Y::AbstractVector{T}
    t_XY::AbstractArray{Complex{T},2}
    ref::ReferenceFrame{T}
    function SpatialLightModulator(x_X, y_Y, t_XY, ref)
        T = promote_type(eltype(x_X), eltype(y_Y), typeof(real(t_XY[1])), typeof(ref.x));
        length(x_X) == size(t_XY, 1) || error("Wrong sizes")
        length(y_Y) == size(t_XY, 2) || error("Wrong sizes")
        new{T}(x_X, y_Y, t_XY, ref)
    end
end

function lightinteraction(slm::SpatialLightModulator{T}, fieldspace::FieldSpace) where T
    tmp_fieldspace = changereferential(fieldspace, slm.ref)
    slm_itp = extrapolate(interpolate((slm.x_X, slm.y_Y), slm.t_XY, Gridded(Constant())), one(T))

    sizeS = size(tmp_fieldspace.e_SXY,1)
    tmp_e_SXY = reshape(tmp_fieldspace.e_SXY, sizeS, :)
    iXY = 1
    @inbounds for yi in tmp_fieldspace.y_Y
        for xi in tmp_fieldspace.x_X
            t = slm_itp(xi,yi)
            for iS in 1:sizeS
                tmp_e_SXY[iS,iXY] *= t
            end
            iXY +=1
        end
    end
    return tmp_fieldspace
end

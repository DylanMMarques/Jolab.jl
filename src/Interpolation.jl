struct Interpolation{T, X <: AbstractVector{T}} <: AbstractOpticalComponent{T}
    x::X
    y::X
end

function coefficient_specific(itp::Interpolation, field::Union{FieldAngularSpectrumScalar, FieldSpaceScalar})
    (fieldl, fieldr) = getfields_lr(itp, field)
    coefficient_specific!(fieldl, fieldr, itp, field)
    return (fieldl, fieldr)
end

function coefficient_specific!(fieldl::FieldSpaceScalar, fieldr::FieldSpaceScalar, itp::Interpolation, space::FieldSpaceScalar{T}) where T
    checkapplicability(fieldl, fieldr, itp, space)
    field = dir(space) > 0 ? fieldr : fieldl

    e_itp = LinearInterpolation((space.x_X, space.y_Y), reshape(space.e_SXY, length(space.x_X), length(space.y_Y)), extrapolation_bc = zero(Complex{T}))

	cart = CartesianIndices(field)
	# @inbounds @simd 
    for i in iterator_index(field)
		field.e_SXY[i] = e_itp(field.x_X[cart[i][2]], field.y_Y[cart[i][3]])
	end
    return (fieldl, fieldr)
end

function coefficient_specific!(fieldl::FieldAngularSpectrumScalar, fieldr::FieldAngularSpectrumScalar, itp::Interpolation, angspe::FieldAngularSpectrumScalar{T}) where T
    checkapplicability(fieldl, fieldr, itp, angspe)
    field = dir(angspe) > 0 ? fieldr : fieldl

    e_itp = LinearInterpolation((angspe.nsx_X, angspe.nsy_Y), reshape(angspe.e_SXY, length(angspe.nsx_X), length(angspe.nsy_Y)), extrapolation_bc = zero(Complex{T}))

	cart = CartesianIndices(field)
	# @inbounds @simd 
    for i in iterator_index(field)
		field.e_SXY[i] = e_itp(field.nsx_X[cart[i][2]], field.nsy_Y[cart[i][3]])
	end
    return (fieldl, fieldr)
end

function getfields_lr(itp::Interpolation{T1, X1}, angspe::FieldAngularSpectrumScalar{T,D,X,Y}) where {T,D,X,Y, T1,X1}
    e_SXY = zeros(Complex{T}, length(itp.x) * length(itp.y))
    if dir(angspe) > 0 
	    fieldl = FieldAngularSpectrumScalar{T,-1,X,Y}(deepcopy(angspe.nsx_X), deepcopy(angspe.nsy_Y), 0deepcopy(angspe.e_SXY), angspe.λ, angspe.n, angspe.ref)
	    fieldr = FieldAngularSpectrumScalar{T,1,X1,Vector{Complex{T}}}(deepcopy(itp.x), deepcopy(itp.y), e_SXY, angspe.λ, angspe.n, angspe.ref)
    else
	    fieldl = FieldAngularSpectrumScalar{T,-1,X1,Vector{Complex{T}}}(deepcopy(itp.x), deepcopy(itp.y), e_SXY, angspe.λ, angspe.n, angspe.ref)
	    fieldr = FieldAngularSpectrumScalar{T,1,X,Y}(deepcopy(angspe.nsx_X), deepcopy(angspe.nsy_Y), 0deepcopy(angspe.e_SXY), angspe.λ, angspe.n, angspe.ref)
    end
    return (fieldl, fieldr)
end

function getfields_lr(itp::Interpolation{T1, X1}, space::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y, T1,X1}
    e_SXY = zeros(Complex{T}, length(itp.x) * length(itp.y))
    if dir(space) > 0 
	    fieldl = FieldSpaceScalar{T,-1,X,Y}(deepcopy(space.x_X), deepcopy(space.y_Y), 0deepcopy(space.e_SXY), space.λ, space.n, space.ref)
	    fieldr = FieldSpaceScalar{T,1,X1,Vector{Complex{T}}}(deepcopy(itp.x), deepcopy(itp.y), e_SXY, space.λ, space.n, space.ref)
    else
	    fieldl = FieldSpaceScalar{T,-1,X1,Vector{Complex{T}}}(deepcopy(itp.x), deepcopy(itp.y), e_SXY, space.λ, space.n, space.ref)
	    fieldr = FieldSpaceScalar{T,1,X,Y}(deepcopy(space.x_X), deepcopy(space.y_Y), 0deepcopy(space.e_SXY), space.λ, space.n, space.ref)
    end
    return (fieldl, fieldr)
end

function checkapplicability(fieldl, fieldr, itp::Interpolation, fieldi::Union{FieldSpaceScalar, FieldAngularSpectrumScalar})
    field = dir(fieldi) > 0 ? fieldr : fieldl
    fieldi.ref == field.ref || error("The interpolated field must have the same reference frame as the original field")
    isapprox(fieldi.n, field.n, atol = @tol) || error("The interpolated field must be on the same refractive index that the original field")
    isapprox(fieldi.λ, field.λ, atol = @tol) || error("The interpolated field must have the same wavelength that the original field")
end

function interpolate(field, x_X, y_Y)
    itp = Interpolation(x_X, y_Y)
    (fieldl, fieldr) = coefficient_specific(itp, field)
    return dir(field) > 0 ? fieldr : fieldl
end
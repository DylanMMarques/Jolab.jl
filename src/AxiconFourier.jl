struct AxiconFourier{T,X<:AbstractVector{T}} <: AbstractOpticalComponent{T}
	x_X::X
	y_Y::X
	nsx_X::X
	nsy_Y::X
    β::T
    ref::ReferenceFrame{T}
	function AxiconFourier{T}(x_X, y_Y, nsx_X::X1, nsy_Y::X2, β, ref) where {T,X1,X2}
		X = promote_type(X1,X2)
		return new{T,X}(x_X, y_Y, nsx_X, nsy_Y, β, ref)
	end
end

AxiconFourier(x_X, y_Y, nsx_X, nsy_Y, β, ref) = AxiconFourier{Float64}(x_X, y_Y, nsx_X, nsy_Y, β, ref)

getsizes(axicon::AxiconFourier, field::AbstractFieldSpace) = return (length(axicon.nsx_X), length(axicon.nsy_Y))
getsizes(axicon::AxiconFourier, field::AbstractFieldAngularSpectrum) = return (length(axicon.x_X), length(axicon.y_Y))

r(x, y) = sqrt(x^2 + y^2)

function t(axicon::AxiconFourier{T}, space::FieldSpaceScalarRadialSymmetric, iX, iY, iA, iB) where T
	return - im / 2π * exp(-im * 2π / space.λ * axicon.β * space.r_R[iX]) * besselj0(2π / space.λ * axicon.nsx_X[iA] * space.r_R[iX]) * space.r_R[iX] * Δvector(space.r_R, iX)
end

function t(axicon::AxiconFourier{T}, angspe::FieldAngularSpectrumScalarRadialSymmetric, iX, iY, iA, iB) where T
	return im * 2π * (2π / angspe.λ)^2 * exp(-im * 2π / angspe.λ * axicon.β * axicon.x_X[iA]) * besselj0(2π / angspe.λ * angspe.nsr_R[iX] * axicon.x_X[iA]) * angspe.nsr_R[iX] * Δvector(angspe.nsr_R, iX)
end

@inline function t(axicon::AxiconFourier{T}, space::AbstractFieldSpace{T}, iX, iY, iA, iB)::Complex{T} where T
	(xmin, xmax) = integralExtremes(space.x_X, iX)
	(ymin, ymax) = integralExtremes(space.y_Y, iY)

	k = 2π / space.λ

	f(x,y) = -k * (axicon.β * √(x^2 + y^2) + axicon.nsx_X[iA] * x + axicon.nsy_Y[iB] * y)
	t12 = 1 / 4π^2 * integrate_exp_xy_x_y(f, xmin, xmax, ymin, ymax)
	return t12
end

@inline function tinv(axicon::AxiconFourier{T}, space::AbstractFieldSpace, iX, iY, iA, iB)::Complex{T} where T
	(nsxmin, nsxmax) = integralExtremes(axicon.nsx_X, iA)
	(nsymin, nsymax) = integralExtremes(axicon.nsy_Y, iB)
	k = 2π / space.λ

	return k^2 * exp(-im * k * axicon.β * √(space.x_X[iX]^2 + space.y_Y[iY]^2)) * integrate_exp_x_y(k * axicon.nsx_X[iA], k * axicon.nsy_Y[iB], zero(T), nsxmin, nsxmax, nsymin, nsymax)
end

@inline function t(axicon::AxiconFourier{T}, angspe::AbstractFieldAngularSpectrum, iX, iY, iA, iB)::Complex{T} where T
	(nsxmin, nsxmax) = integralExtremes(angspe.nsx_X, iX)
	(nsymin, nsymax) = integralExtremes(angspe.nsy_Y, iY)

	k = 2π / angspe.λ

	return k^2 * exp(-im * k * axicon.β * √(axicon.x_X[iA]^2 + axicon.y_Y[iB]^2)) * integrate_exp_x_y(k * axicon.x_X[iA], k * axicon.y_Y[iB], zero(T), nsxmin, nsxmax, nsymin, nsymax)
end

@inline function tinv(axicon::AxiconFourier{T}, angspe::AbstractFieldAngularSpectrum, iX, iY, iA, iB)::Complex{T} where T
	(xmin, xmax) = integralExtremes(axicon.x_X, iA)
	(ymin, ymax) = integralExtremes(axicon.y_Y, iB)
	k = 2π / angspe.λ

	f(x,y) = -k * (axicon.β * √(x^2 + y^2) + angspe.nsx_X[iX] * x + angspe.nsy_Y[iY] * y)
	t12 = 1 / 4π^2 * integrate_exp_xy_x_y(f, xmin, xmax, ymin, ymax)
	return t12
end

function coefficient_general(axicon::AxiconFourier, field::AbstractFieldMonochromatic{T}) where T
	checkapplicability(axicon, field)

	scat = get_scatteringmatrixtype(axicon, field)
	field_s = dir(field) > 0 ? scat.fieldr : scat.fieldl

	cart_i = CartesianIndices(field)
	cart_s = CartesianIndices(field_s)

	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let scat = scat, axicon = axicon, field = field
		@inbounds Threads.@threads for i_i in iterator_index(field)
			@simd for i_s in iterator_index(field_s)
				if dir(field) > 0
					scat.t₁₂[i_s, i_i] = t(axicon, field, cart_i[i_i][2], cart_i[i_i][3], cart_s[i_s][2], cart_s[i_s][3])
					scat.t₂₁[i_i, i_s] = tinv(axicon, field, cart_i[i_i][2], cart_i[i_i][3], cart_s[i_s][2], cart_s[i_s][3])
				else
					scat.t₁₂[i_s, i_i] = tinv(axicon, field, cart_i[i_i][2], cart_i[i_i][3], cart_s[i_s][2], cart_s[i_s][3])
					scat.t₂₁[i_i, i_s] = t(axicon, field, cart_i[i_i][2], cart_i[i_i][3], cart_s[i_s][2], cart_s[i_s][3])
				end
			end
		end
	end
	correctscatteringmatrix_referenceframes!(scat, axicon, field)
	return scat
end

function lightinteraction(axicon::AxiconFourier{T,X}, field::AbstractFieldMonochromatic{T}) where {T<:Real,X}
	fieldi_newref = changereferenceframe(field, axicon.ref)
	checkapplicability(axicon, fieldi_newref)


	(fieldl, fieldr) = getfields_lr(axicon, field)
	fieldl.e_SXY .= zero(Complex{T})
	fieldr.e_SXY .= zero(Complex{T})
	field_s = dir(field) > 0 ? fieldr : fieldl

	cart_i = CartesianIndices(field)
	cart_s = CartesianIndices(field_s)

	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let axicon = axicon, fieldi_newref = fieldi_newref, fieldl = fieldl, fieldr = fieldr
		@inbounds Threads.@threads 	for i_s in iterator_index(field_s)
			for i_i in iterator_index(field)
				field_s.e_SXY[i_s] += fieldi_newref.e_SXY[i_i] * t(axicon, fieldi_newref, cart_i[i_i][2], cart_i[i_i][3], cart_s[i_s][2], cart_s[i_s][3])
			end
		end
	end
	return (fieldl, fieldr)
end

function checkapplicability(axicon::AxiconFourier, space::AbstractFieldSpace)
	checkinplane(axicon.ref, space.ref) || error("The axicon and the field must be defined on the same plane")
	checkorientation(axicon.ref, space.ref) || error("The axicon and the field must be defined on reference frames with the same orientation")
	return true
end

function checkapplicability(axicon::AxiconFourier, angspe::AbstractFieldAngularSpectrum)
	return true
end

function get_scatteringmatrixtype(axicon::AxiconFourier{T,Y}, fieldi::FieldSpaceScalarRadialSymmetric{T,D,X,B}) where {T,D,X,Y,B}
	sizeI = length(fieldi.e_SXY)
	sizeS = length(axicon.nsx_X)

	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeS, sizeI)
	t21 = Matrix{Complex{T}}(undef, sizeI, sizeS)
	m = Vector{Complex{T}}(undef, sizeS)
	if dir(fieldi) > 0
		fieldl = FieldSpaceScalarRadialSymmetric{T,-1,X,B}(deepcopy(fieldi.r_R), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrumScalarRadialSymmetric{T,1,Y,B}(deepcopy(axicon.nsx_X), m, fieldi.λ, fieldi.n, axicon.ref)
		return ScatteringMatrix{T,FieldSpaceScalarRadialSymmetric{T,-1,X,B}, FieldAngularSpectrumScalarRadialSymmetric{T,1,Y,B}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldAngularSpectrumScalarRadialSymmetric{T,-1,Y,B}(deepcopy(axicon.nsx_X), m, fieldi.λ, fieldi.n, axicon.ref)
		fieldr = FieldSpaceScalarRadialSymmetric{T,1,X,B}(deepcopy(fieldi.r_R), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldAngularSpectrumScalarRadialSymmetric{T,-1,Y,B}, FieldSpaceScalarRadialSymmetric{T,1,X,B}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	end
end

function get_scatteringmatrixtype(axicon::AxiconFourier{T,Y}, fieldi::FieldAngularSpectrumScalarRadialSymmetric{T,D,X,B}) where {T,D,X,Y,B}
	sizeI = length(fieldi.e_SXY)
	sizeS = length(axicon.r_R)

	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeS, sizeI)
	t21 = Matrix{Complex{T}}(undef, sizeI, sizeS)

	m = Vector{Complex{T}}(undef, sizeS)
	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrumScalarRadialSymmetric{T,-1,Y,B}(deepcopy(fieldi.nsr_R), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpaceScalarRadialSymmetric{T,1,X,B}(deepcopy(axicon.x_X), m, fieldi.λ, fieldi.n, axicon.ref)
		return ScatteringMatrix{T, FieldAngularSpectrumScalarRadialSymmetric{T,-1,Y,B}, FieldSpaceScalarRadialSymmetric{T,1,X,B}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldSpaceScalarRadialSymmetric{T,-1,X,B}(deepcopy(axicon.x_X), m, fieldi.λ, fieldi.n, axicon.ref)
		fieldr = FieldAngularSpectrumScalarRadialSymmetric{T,1,Y,B}(deepcopy(fieldi.nsr_R), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldSpaceScalarRadialSymmetric{T,-1,X,B}, FieldAngularSpectrumScalarRadialSymmetric{T,1,Y,B}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	end
end

function get_scatteringmatrixtype(axicon::AxiconFourier{T,Y}, fieldi::FieldSpaceScalar{T,D,X,B}) where {T,D,X,Y,B}
	sizeI = length(fieldi.e_SXY)
	sizeS = length(axicon.nsx_X) * length(axicon.nsy_Y)

	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeS, sizeI)
	t21 = Matrix{Complex{T}}(undef, sizeI, sizeS)
	m = Vector{Complex{T}}(undef, sizeS)
	if dir(fieldi) > 0
		fieldl = FieldSpaceScalar{T,-1,X,B}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrumScalar{T,1,Y,B}(deepcopy(axicon.nsx_X), deepcopy(axicon.nsy_Y), m, fieldi.λ, fieldi.n, axicon.ref)
		return ScatteringMatrix{T,FieldSpaceScalar{T,-1,X,B}, FieldAngularSpectrumScalar{T,1,Y,B},Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldAngularSpectrumScalar{T,-1,Y,B}(deepcopy(axicon.nsx_X), deepcopy(axicon.nsy_Y), m, fieldi.λ, fieldi.n, axicon.ref)
		fieldr = FieldSpaceScalar{T,1,X,B}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldAngularSpectrumScalar{T,-1,Y,B}, FieldSpaceScalar{T,1,X,B}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	end
end

function get_scatteringmatrixtype(axicon::AxiconFourier{T,Y}, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,D,X,Y,B}
	sizeI = length(fieldi.e_SXY)
	sizeS = length(axicon.x_X) * length(axicon.y_Y)

	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeS, sizeI)
	t21 = Matrix{Complex{T}}(undef, sizeI, sizeS)

	m = Vector{Complex{T}}(undef, sizeS)
	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrumScalar{T,-1,Y,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpaceScalar{T,1,X,B}(deepcopy(axicon.x_X), deepcopy(axicon.y_Y), m, fieldi.λ, fieldi.n, axicon.ref)
		return ScatteringMatrix{T, FieldAngularSpectrumScalar{T,-1,Y,B}, FieldSpaceScalar{T,1,X,B}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldSpaceScalar{T,-1,X,B}(deepcopy(axicon.x_X), deepcopy(axicon.y_Y), m, fieldi.λ, fieldi.n, axicon.ref)
		fieldr = FieldAngularSpectrumScalar{T,1,Y,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldSpaceScalar{T,-1,X,B}, FieldAngularSpectrumScalar{T,1,Y,B}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	end
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldSpaceScalar{T,1,X,B}) where {T,X,A,B}
	m = Vector{Complex{T}}(undef,length(axicon.nsx_X) * length(axicon.nsy_Y))
	fieldl = FieldSpaceScalar{T,-1,X,B}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldAngularSpectrumScalar{T,1,A,B}(deepcopy(axicon.nsx_X), deepcopy(axicon.nsy_Y), m, fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldSpaceScalar{T,-1,X,B}) where {T,X,A,B}
	m = Vector{Complex{T}}(undef,length(axicon.nsx_X) * length(axicon.nsy_Y))
	fieldl = FieldAngularSpectrumScalar{T,-1,A,B}(deepcopy(axicon.nsx_X), deepcopy(axicon.nsy_Y), m, fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldSpaceScalar{T,1,X,B}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldAngularSpectrumScalar{T,1,X,B}) where {T,X,A,B}
	m = Vector{Complex{T}}(undef,length(axicon.x_X) * length(axicon.y_Y))
	fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldSpaceScalar{T,1,A,B}(deepcopy(axicon.x_X), deepcopy(axicon.y_Y), m, fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldAngularSpectrumScalar{T,-1,X,B}) where {T,X,A,B}
	m = Vector{Complex{T}}(undef,length(axicon.x_X) * length(axicon.y_Y))
	fieldl = FieldSpaceScalar{T,-1,A,B}(deepcopy(axicon.x_X), deepcopy(axicon.y_Y), m, fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldAngularSpectrumScalar{T,1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldSpaceScalarRadialSymmetric{T,1,X,B}) where {T,X,A,B}
	m = Vector{Complex{T}}(undef,length(axicon.nsx_X))
	fieldl = FieldSpaceScalarRadialSymmetric{T,-1,X,B}(deepcopy(fieldi.r_R), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldAngularSpectrumScalarRadialSymmetric{T,1,A,B}(deepcopy(axicon.nsx_X), m, fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldSpaceScalarRadialSymmetric{T,-1,X,B}) where {T,X,A,B}
	m = Vector{Complex{T}}(undef,length(axicon.nsx_X))
	fieldl = FieldAngularSpectrumScalarRadialSymmetric{T,-1,A,B}(deepcopy(axicon.nsx_X), m, fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldSpaceScalarRadialSymmetric{T,1,X,B}(deepcopy(fieldi.r_R), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldAngularSpectrumScalarRadialSymmetric{T,1,X,B}) where {T,X,A,B}
	m = Vector{Complex{T}}(undef,length(axicon.x_X))
	fieldl = FieldAngularSpectrumScalarRadialSymmetric{T,-1,X,B}(deepcopy(fieldi.nsr_R), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldSpaceScalarRadialSymmetric{T,1,A,B}(deepcopy(axicon.x_X), m, fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldAngularSpectrumScalarRadialSymmetric{T,-1,X,B}) where {T,X,A,B}
	m = Vector{Complex{T}}(undef,length(axicon.x_X))
	fieldl = FieldSpaceScalarRadialSymmetric{T,-1,A,B}(deepcopy(axicon.x_X), m, fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldAngularSpectrumScalarRadialSymmetric{T,1,X,B}(deepcopy(fieldi.nsr_R), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

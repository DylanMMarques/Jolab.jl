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

getsizes(axicon::AxiconFourier, field::FieldSpace) = return (length(axicon.nsx_X), length(axicon.nsy_Y))
getsizes(axicon::AxiconFourier, field::FieldAngularSpectrum) = return (length(axicon.x_X), length(axicon.y_Y))

function r(x::T, y::T)::T where T
	return sqrt(x^2 + y^2)
end

@inline function t(axicon::AxiconFourier{T}, space::FieldSpace{T}, iX, iY, iA, iB)::Complex{T} where T
	(xmin, xmax) = integralExtremes(space.x_X, iX)
	(ymin, ymax) = integralExtremes(space.y_Y, iY)

	k = 2π / space.λ

	bool = space.x_X[iX]^2 + space.y_Y[iY]^2 > 3E-6 && xmin * xmax > 0 && ymin * ymax > 0
	# @show space.x_X[iX]^2 + space.y_Y[iY]^2 > 3E-6
	# @show xmin * xmax > 0 && ymin * ymax > 0
	if bool
		@inline f(x,y) = -k * (axicon.β * √(x^2 + y^2) + axicon.nsx_X[iA] * x + axicon.nsy_Y[iB] * y)
		t12 = 1 / 4π^2 * integrate_exp_xy_x_y(f, xmin, xmax, ymin, ymax)
	else
		@inline phaseTerm(r) = exp(-im * k * (axicon.β * √(r[1]^2 + r[2]^2) + axicon.nsx_X[iA] * r[1] + axicon.nsy_Y[iB] * r[2]))
	 	t12 = 1 / 4π^2 * hcubature(phaseTerm, SVector(xmin, ymin), SVector(xmax, ymax), rtol = 1E-2, maxevals = 2000)[1]
	end
	return t12
end

@inline function tinv(axicon::AxiconFourier{T}, space::FieldSpace, iX, iY, iA, iB)::Complex{T} where T
	(nsxmin, nsxmax) = integralExtremes(axicon.nsx_X, iA)
	(nsymin, nsymax) = integralExtremes(axicon.nsy_Y, iB)
	k = 2π / space.λ

	return k^2 * exp(-im * k * axicon.β * √(space.x_X[iX]^2 + space.y_Y[iY]^2)) * integrate_exp_x_y(k * axicon.nsx_X[iA], k * axicon.nsy_Y[iB], zero(T), nsxmin, nsxmax, nsymin, nsymax)
end

@inline function t(axicon::AxiconFourier{T}, angspe::FieldAngularSpectrum, iX, iY, iA, iB)::Complex{T} where T
	(nsxmin, nsxmax) = integralExtremes(angspe.nsx_X, iX)
	(nsymin, nsymax) = integralExtremes(angspe.nsy_Y, iY)

	k = 2π / angspe.λ

	return k^2 * exp(-im * k * axicon.β * √(axicon.x_X[iA]^2 + axicon.y_Y[iB]^2)) * integrate_exp_x_y(k * axicon.x_X[iA], k * axicon.y_Y[iB], zero(T), nsxmin, nsxmax, nsymin, nsymax)
end

@inline function tinv(axicon::AxiconFourier{T}, angspe::FieldAngularSpectrum, iX, iY, iA, iB)::Complex{T} where T
	(xmin, xmax) = integralExtremes(axicon.x_X, iA)
	(ymin, ymax) = integralExtremes(axicon.y_Y, iB)
	k = 2π / angspe.λ

	bool = axicon.x_X[iA]^2 + axicon.y_Y[iB]^2 > 3E-6 && xmin * xmax > 0 && ymin * ymax > 0
	# Forward case
	if bool
		@inline f(x,y) = -k * (axicon.β * √(x^2 + y^2) + angspe.nsx_X[iX] * x + angspe.nsy_Y[iY] * y)
		t12 = 1 / 4π^2 * integrate_exp_xy_x_y(f, xmin, xmax, ymin, ymax)
	else
		@inline phaseTerm(r) = exp(-im * k * (axicon.β * √(r[1]^2 + r[2]^2) + angspe.nsx_X[iX] * r[1] + axicon.nsy_Y[iY] * r[2]))
	 	t12 = 1 / 4π^2 * hcubature(phaseTerm, SVector(xmin, ymin), SVector(xmax, ymax), rtol = 1E-2, maxevals=100)[1]
	end
	return t12
end

function coefficient_general(axicon::AxiconFourier, field::AbstractFieldMonochromatic{T}) where T
	checkapplicability(axicon, field)

	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(axicon, field)

	coordXY = LinearIndices((sizeX, sizeY))
	coordAB = LinearIndices((sizeA, sizeB))

	scat = get_scatteringmatrixtype(axicon, field)

	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let scat = scat, axicon = axicon, field = field
		@inbounds Threads.@threads for iX1 in 1:sizeX
	 	 	@simd for iY1 in 1:sizeY
				i1 = coordXY[iX1, iY1]
				for iA2 in 1:sizeA
					for iB2 in 1:sizeB
						i2 = coordAB[iA2, iB2]
						if dir(field) > 0
							scat.t₁₂[i2, i1] = t(axicon, field, iX1, iY1, iA2, iB2)
							scat.t₂₁[i1, i2] = tinv(axicon, field, iX1, iY1, iA2, iB2)
						else
							scat.t₁₂[i2, i1] = tinv(axicon, field, iX1, iY1, iA2, iB2)
							scat.t₂₁[i1, i2] = t(axicon, field, iX1, iY1, iA2, iB2)
						end
					end
				end
			end
		end
	end
	correctscatteringmatrix_referenceframes!(scat, axicon, field)
	return scat
end

function lightinteraction(axicon::AxiconFourier{T,X}, field::AbstractFieldMonochromatic{T}) where {T<:Real,X}
	checkapplicability(axicon, field)
	fieldi_newref = changereferenceframe(field, axicon.ref)

	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(axicon, field)
	coordXY = LinearIndices((sizeX, sizeY))
	coordAB = LinearIndices((sizeA, sizeB))

	(fieldl, fieldr) = getfields_lr(axicon, field)
	fieldl.e_SXY .= zero(Complex{T})
	fieldr.e_SXY .= zero(Complex{T})

	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let axicon = axicon, fieldi_newref = fieldi_newref, fieldl = fieldl, fieldr = fieldr
		@inbounds Threads.@threads for iA2 in 1:sizeA
			for iB2 in 1:sizeB
				i2 = coordAB[iA2, iB2]
				for iX1 in 1:sizeX
	 	 			for iY1 in 1:sizeY
						i1 = coordXY[iX1, iY1]
						if dir(fieldi_newref) > 0
							fieldr.e_SXY[i2] += fieldi_newref.e_SXY[i1] * t(axicon, field, iX1, iY1, iA2, iB2)
						else
							fieldl.e_SXY[i2] += fieldi_newref.e_SXY[i1] * t(axicon, field, iX1, iY1, iA2, iB2)
						end
					end
				end
			end
		end
	end
	return (fieldl, fieldr)
end

function checkapplicability(axicon::AxiconFourier, space::FieldSpace)
	checkinplane(axicon.ref, space.ref) || error("The axicon and the field must be defined on the same plane")
	checkorientation(axicon.ref, space.ref) || error("The axicon and the field must be defined on reference frames with the same orientation")
	return true
end

function checkapplicability(axicon::AxiconFourier, angspe::FieldAngularSpectrum)
	return true
end

function get_scatteringmatrixtype(axicon::AxiconFourier{T,Y}, fieldi::FieldSpace{T,D,X}) where {T,D,X,Y}
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	sizeA, sizeB = length(axicon.nsx_X), length(axicon.nsy_Y)

	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)
	m = Array{Complex{T},3}(undef,1,length(axicon.nsx_X), length(axicon.nsy_Y))
	if dir(fieldi) > 0
		fieldl = FieldSpace{T,-1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrum{T,1,Y}(deepcopy(axicon.nsx_X), deepcopy(axicon.nsy_Y), m, fieldi.λ, fieldi.n, axicon.ref)
		return ScatteringMatrix{T,FieldSpace{T,-1,X}, FieldAngularSpectrum{T,1,Y},Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldAngularSpectrum{T,-1,Y}(deepcopy(axicon.nsx_X), deepcopy(axicon.nsy_Y), m, fieldi.λ, fieldi.n, axicon.ref)
		fieldr = FieldSpace{T,1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldAngularSpectrum{T,-1,Y}, FieldSpace{T,1,X}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	end
end

function get_scatteringmatrixtype(axicon::AxiconFourier{T,Y}, fieldi::FieldAngularSpectrum{T,D,X}) where {T,D,X,Y}
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	sizeA, sizeB = length(axicon.x_X), length(axicon.y_Y)

	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)

	m = Array{Complex{T},3}(undef,1,length(axicon.x_X), length(axicon.y_Y))
	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrum{T,-1,Y}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpace{T,1,X}(deepcopy(axicon.x_X), deepcopy(axicon.y_Y), m, fieldi.λ, fieldi.n, axicon.ref)
		return ScatteringMatrix{T, FieldAngularSpectrum{T,-1,Y}, FieldSpace{T,1,X}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldSpace{T,-1,X}(deepcopy(axicon.x_X), deepcopy(axicon.y_Y), m, fieldi.λ, fieldi.n, axicon.ref)
		fieldr = FieldAngularSpectrum{T,1,Y}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldSpace{T,-1,X}, FieldAngularSpectrum{T,1,Y}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	end
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldSpace{T,1,X}) where {T,X,A}
	m = Array{Complex{T},3}(undef,1,length(axicon.nsx_X), length(axicon.nsy_Y))
	fieldl = FieldSpace{T,-1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldAngularSpectrum{T,1,A}(deepcopy(axicon.nsx_X), deepcopy(axicon.nsy_Y), m, fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldSpace{T,-1,X}) where {T,X,A}
	m = Array{Complex{T},3}(undef,1,length(axicon.nsx_X), length(axicon.nsy_Y))
	fieldl = FieldAngularSpectrum{T,-1,A}(deepcopy(axicon.nsx_X), deepcopy(axicon.nsy_Y), m, fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldSpace{T,1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldAngularSpectrum{T,1,X}) where {T,X,A}
	m = Array{Complex{T},3}(undef, 1, length(axicon.x_X), length(axicon.y_Y))
	fieldl = FieldAngularSpectrum{T,-1,X}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldSpace{T,1,A}(deepcopy(axicon.x_X), deepcopy(axicon.y_Y), m, fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

function getfields_lr(axicon::AxiconFourier{T,A}, fieldi::FieldAngularSpectrum{T,-1,X}) where {T,X,A}
	m = Array{Complex{T},3}(undef, 1, length(axicon.x_X), length(axicon.y_Y))
	fieldl = FieldSpace{T,-1,A}(deepcopy(axicon.x_X), deepcopy(axicon.y_Y), m, fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldAngularSpectrum{T,1,X}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

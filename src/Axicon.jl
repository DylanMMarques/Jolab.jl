struct Axicon{T} <: AbstractOpticalComponent{T}
    β::T
    ref::ReferenceFrame{T}
	function Axicon{T}(β, ref) where T
		return new{T}(β, ref)
	end
end

Axicon(β, ref) = Axicon{Float64}(β, ref)

@inline function t(axicon::Axicon{T}, space::FieldSpace, iX, iY) where T
	k = 2π / space.λ
	return exp(-im * k * axicon.β * √(space.x_X[iX]^2 + space.y_Y[iY]^2))
end

function checkapplicability(axicon::Axicon, space::FieldSpace)
	checkinplane(axicon.ref, space.ref) || error("The axicon and the field must be defined on the same plane")
	checkorientation(axicon.ref, space.ref) || error("The axicon and the field must be defined on reference frames with the same orientation")
end

function coefficient_general(axicon::Axicon{T}, space::FieldSpace{T}, nsx_A = range(-2axicon.β, 2axicon.β, length = size(space.e_SXY)[2]), nsy_B = range(-2axicon.β, 2axicon.β, length = size(space.e_SXY)[3])) where T
	checkapplicability(axicon, space)

	scat = get_scatteringmatrixtype(axicon, space)

	(sizeX, sizeY) = size(space.e_SXY)[2:3]
	coordXY = LinearIndices((sizeX, sizeY))

	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let scat = scat, axicon = axicon, space = space
		@inbounds Threads.@threads for iX1 in 1:sizeX
	 	 	@simd for iY1 in 1:sizeY
				i1 = coordXY[iX1, iY1]
				scat.t₁₂.diag[i1] = t(axicon, space, iX1, iY1)
			end
		end
	end
	return scat
end

function lightinteraction(axicon::Axicon, space::FieldSpace{T}) where T
	checkapplicability(axicon, space)

	(sizeX, sizeY) = size(space.e_SXY)[2:3]
	coordXY = LinearIndices((sizeX, sizeY))

	(fieldl, fieldr) = getfields_lr(axicon, space)

	e_SXY = dir(space) > 0 ? fieldr.e_SXY : fieldl.e_SXY
	dir(space) > 0 ? fieldl.e_SXY .= zero(Complex{T}) : fieldr.e_SXY .= zero(Complex{T})
	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let axicon = axicon, space = space
		@inbounds Threads.@threads for iX1 in 1:sizeX
	 	 	@simd for iY1 in 1:sizeY
				i1 = coordXY[iX1, iY1]
				e_SXY[i1] = space.e_SXY[i1] * t(axicon, space, iX1, iY1)
			end
		end
	end
	return (fieldl, fieldr)
end

function get_scatteringmatrixtype(axicon::Axicon{T}, fieldi::FieldSpace{T,D,X}) where {T,D,X}
	sizeXY = length(fieldi.x_X) * length(fieldi.y_Y)
	r12 = nothing
	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r21 = nothing
	t21 = t12
	sizeX = length(fieldi.x_X)
	sizeY = length(fieldi.y_Y)

	(fieldl, fieldr) = getfields_lr(axicon, fieldi)

	return ScatteringMatrix{T,FieldSpace{T,-1,X},FieldSpace{T,1,X}, Nothing, Diagonal{Complex{T},Vector{Complex{T}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

function getfields_lr(axicon::Axicon, fieldi::FieldSpace{T,D,X}) where {T,D,X}
	fieldl = FieldSpace{T,-1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldSpace{T,1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

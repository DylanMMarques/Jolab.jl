struct Axicon{T} <: AbstractOpticalComponent{T}
    β::T
    ref::ReferenceFrame{T}
	function Axicon{T}(β, ref) where T
		return new{T}(β, ref)
	end
end

Axicon(β, ref) = Axicon{Float64}(β, ref)

@inline function t(axicon::Axicon{T}, space::AbstractFieldSpace, iX, iY) where T
	k = 2π / space.λ
	return exp(-im * k * axicon.β * √(space.x_X[iX]^2 + space.y_Y[iY]^2))
end

function checkapplicability(axicon::Axicon, space::AbstractFieldSpace)
	checkinplane(axicon.ref, space.ref) || error("The axicon and the field must be defined on the same plane")
	checkorientation(axicon.ref, space.ref) || error("The axicon and the field must be defined on reference frames with the same orientation")
	axicon.ref == space.ref || error("Need to be done the transformation")
end

function coefficient_general(axicon::Axicon{T}, space::AbstractFieldSpace{T}) where T
	checkapplicability(axicon, space)

	scat = get_scatteringmatrixtype(axicon, space)

	cart = CartesianIndices(space)
	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let scat = scat, axicon = axicon, space = space
		@inbounds Threads.@threads for i in iterator_index(scat.fieldl)
			scat.t₁₂.diag[i] = t(axicon, space, cart[i][2], cart[i][3])
		end
	end
	return scat
end

function lightinteraction(axicon::Axicon, space::AbstractFieldSpace{T}) where T
	checkapplicability(axicon, space)

	(fieldl, fieldr) = getfields_lr(axicon, space)

	cart = CartesianIndices(space)
	fieldt = dir(space) > 0 ? fieldr : fieldl
	dir(space) > 0 ? fieldl.e_SXY .= zero(Complex{T}) : fieldr.e_SXY .= zero(Complex{T})
	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let axicon = axicon, space = space
		@inbounds Threads.@threads for i in iterator_index(fieldt)
			fieldt.e_SXY[i] = space.e_SXY[i] * t(axicon, space, cart[i][2], cart[i][3])
		end
	end
	return (fieldl, fieldr)
end

function get_scatteringmatrixtype(axicon::Axicon{T}, fieldi::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y}
	sizeXY = length(fieldi.e_SXY)
	r12 = nothing
	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r21 = nothing
	t21 = t12

	(fieldl, fieldr) = getfields_lr(axicon, fieldi)

	return ScatteringMatrix{T,FieldSpaceScalar{T,-1,X,Y},FieldSpaceScalar{T,1,X,Y}, Nothing, Diagonal{Complex{T},Vector{Complex{T}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

function getfields_lr(axicon::Axicon, fieldi::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y}
	fieldl = FieldSpaceScalar{T,-1,X,Y}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	fieldr = FieldSpaceScalar{T,1,X,Y}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, axicon.ref)
	return (fieldl, fieldr)
end

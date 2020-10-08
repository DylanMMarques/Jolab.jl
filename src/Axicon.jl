struct Axicon{T} <: AbstractOpticalComponent{T}
    β::T
    α::T
	dir::Int16
    ref::ReferenceFrame{T}
end

n(axicon::Axicon) = sin(axicon.β + axicon.α) / sin(axicon.α)

function coefficient_general(axicon::Axicon{T}, space::FieldSpace{T}) where T
	checkorientation(space.ref, axicon.ref) || errorToDo()
	checkinplane(space.ref, axicon.ref) || errorToDo()
	space.dir * axicon.dir < 0 && errorToDo()

	sizeXY = length(space.x_X) * length(space.y_Y)
	(sizeX, sizeY) = (length(space.x_X), length(space.y_Y))
	r = UniformScaling(zero(Complex{T}))

	k = 2π / space.λ
	t12 = Matrix{Complex{T}}(undef, sizeXY, sizeXY)
	t21 = Matrix{Complex{T}}(undef, sizeXY, sizeXY)
	cart = LinearIndices((sizeX, sizeY))

	nsx_X = range(-axicon.β * 2, axicon.β * 2, length = sizeX)
	nsy_Y = range(-axicon.β * 2, axicon.β * 2, length = sizeY)

	Threads.@threads for ixN in eachindex(space.x_X)
		(xmin, xmax) = integralExtremes(space.x_X, ixN)
		for iyN in eachindex(space.y_Y)
			(ymin, ymax) = integralExtremes(space.y_Y, iyN)
			iN = cart[ixN, iyN]
			bool =  space.x_X[ixN]^2 + space.y_Y[iyN]^2 > 1E-6 && xmin * xmax > 0 && ymin * ymax > 0
			for ixM in 1:sizeX
				for iyM in 1:sizeY
					# hcubature(phaseTerm, SVector(xmin, ymin), SVector(xmax, ymax), atol = 1E-5)[1]
					# t12[cart[ixM, iyM], iN] = 1 / 4π^2 * hcubature(phaseTerm, SVector(xmin, ymin), SVector(xmax, ymax), atol = 1E-5)[1]
					if bool
						f(x,y) = -k * (axicon.β * √(x^2 + y^2) + nsx_X[ixM] * x + nsy_Y[iyM] * y)
						aux = integral(f, xmin, xmax, ymin, ymax)
					else
						phaseTerm(r) = exp(-im * k * (axicon.β * √(r[1]^2 + r[2]^2) + nsx_X[ixM] * r[1] + nsy_Y[iyM] * r[2]))
					 	aux = hcubature(phaseTerm, SVector(xmin, ymin), SVector(xmax, ymax), atol = 1E-5)[1]
					end
					t12[cart[ixM, iyM], iN] = 1 / 4π^2 * aux
				end
			end
		end
	end

	if space.dir > 0
		fieldl = FieldSpace{T}(copy(space.x_X), copy(space.y_Y), space.e_SXY, space.λ, space.n, -1, space.ref)
		fieldr = FieldAngularSpectrum{T}(nsx_X, nsy_Y, copy(space.e_SXY), space.λ, space.n, 1, space.ref)
	else
		fieldr = FieldSpace{T}(copy(space.x_X), copy(space.y_Y), space.e_SXY, space.λ, space.n, 11, space.ref)
		fieldl = FieldAngularSpectrum{T}(nsx_X, nsy_Y, copy(space.e_SXY), space.λ, space.n, 1, space.ref)
	end
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), typeof(r), Matrix{Complex{T}}}(r, t12, r, t12, fieldl, fieldr)
end

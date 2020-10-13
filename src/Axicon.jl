struct Axicon{T} <: AbstractOpticalComponent{T}
    β::T
    ref::ReferenceFrame{T}
	function Axicon{T}(β, ref) where T
		return new{T}(β, ref)
	end
end

Axicon(β, ref) = Axicon{Float64}(β, ref)

@inline function axiconΔintegral(axicon::Axicon{T}, space::FieldSpace, nsx_A, nsy_B, iX, iY, iA, iB)::Tuple{Complex{T}, Complex{T}} where T
	(xmin, xmax) = integralExtremes(space.x_X, iX)
	(ymin, ymax) = integralExtremes(space.y_Y, iY)

	(nsxmin, nsxmax) = integralExtremes(nsx_A, iA)
	(nsymin, nsymax) = integralExtremes(nsy_B, iB)
	k = 2π / space.λ

	bool =  space.x_X[iX]^2 + space.y_Y[iY]^2 > 3E-6 && xmin * xmax > 0 && ymin * ymax > 0
	# Forward case
	if bool
		@inline f(x,y) = -k * (axicon.β * √(x^2 + y^2) + nsx_A[iA] * x + nsy_B[iB] * y)
		t12 = 1 / 4π^2 * integrate_exp_xy_x_y(f, xmin, xmax, ymin, ymax)
	else
		@inline phaseTerm(r) = exp(-im * k * (axicon.β * √(r[1]^2 + r[2]^2) + nsx_A[iA] * r[1] + nsy_B[iB] * r[2]))
	 	t12 = 1 / 4π^2 * hcubature(phaseTerm, SVector(xmin, ymin), SVector(xmax, ymax), rtol = 1E-2)[1]
	end
	t21 = k^2 * exp(-im * k * axicon.β * √(space.x_X[iX]^2 + space.y_Y[iY]^2)) * integrate_exp_x_y(k * space.x_X[iX], k * space.y_Y[iY], zero(T), nsxmin, nsxmax, nsymin, nsymax)
	return (t12, t21)
end

function coefficient_general(axicon::Axicon{T}, space::FieldSpace{T}, nsx_A = range(-2axicon.β, 2axicon.β, length = size(space.e_SXY)[2]), nsy_B = range(-2axicon.β, 2axicon.β, length = size(space.e_SXY)[3])) where T
	checkorientation(space.ref, axicon.ref) || tobedone()
	checkinplane(space.ref, axicon.ref) || tobedone()

	(sizeX, sizeY) = size(space.e_SXY)[2:3]
	sizeA, sizeB = length(nsx_A), length(nsy_B)
	r12 = UniformScaling(zero(Complex{T}))
	r21 = UniformScaling(zero(Complex{T}))
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)

	coordXY = LinearIndices((sizeX, sizeY))
	coordAB = LinearIndices((sizeA, sizeB))

	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let t12 = t12, t21 = t21, axicon = axicon, space = space
		@inbounds Threads.@threads for iX1 in 1:sizeX
	 	 	@simd for iY1 in 1:sizeY
				i1 = coordXY[iX1, iY1]
				for iA2 in 1:sizeA
					for iB2 in 1:sizeB
						i2 = coordAB[iA2, iB2]
						(t12[i2, i1], t21[i1,i2])  = axiconΔintegral(axicon, space, nsx_A, nsy_B, iX1, iY1, iA2, iB2)
					end
				end
			end
		end
	end
	if space.dir > 0
		fieldl = FieldSpace{T}(copy(space.x_X), copy(space.y_Y), space.e_SXY, space.λ, space.n, -1, copy(space.ref))
		fieldr = FieldAngularSpectrum{T}(nsx_A, nsy_B, zeros(Complex{T},1,sizeA,sizeB), space.λ, space.n, 1, copy(space.ref))
	else
		fieldl = FieldAngularSpectrum{T}(nsx_A, nsy_B, zeros(Complex{T},1,sizeA,sizeB), space.λ, space.n, -1, copy(space.ref))
		fieldr = FieldSpace{T}(copy(space.x_X), copy(space.y_Y), space.e_SXY, space.λ, space.n, 1, copy(space.ref))
		aux = t12
		t12 = t21
		t21 = aux
	end
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), typeof(r12), typeof(t12)}(r12, t12, r21, t21, fieldl, fieldr)
end

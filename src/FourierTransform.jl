struct FourierTransform{T} <: AbstractOpticalComponent{T}
end

FourierTransform() = FourierTransform{Float64}()

function coefficient_general(fourier::FourierTransform{T}, field::FieldSpace{T,X}, nsx_A::AbstractRange{<:Real}, nsy_B::AbstractRange{<:Real}) where {T, X <: AbstractRange}
	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	(sizeA, sizeB) = (length(nsx_A), length(nsy_B))
	r12 = spzeros(Complex{T}, sizeX*sizeY, sizeX * sizeY)
	r21 = spzeros(Complex{T}, sizeA*sizeB, sizeA * sizeB)
	t12 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)
	t21 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)

	k = 2π / field.λ
	imk = im * k
	Δx = field.x_X[2] - field.x_X[1]
	Δy = field.y_Y[2] - field.y_Y[1]
	ΔxΔy = Δx * Δy / 4π^2

	Δnsx = nsx_A[2] - nsx_A[1]
	Δnsy = nsy_B[2] - nsy_B[1]

	ΔkxΔky = k^2 * Δnsx * Δnsy
	coord = LinearIndices((sizeA, sizeB))
	@inbounds Threads.@threads for iX1 in 1:sizeX
	 	@simd for iY1 in 1:sizeY
			i1 = coord[iX1, iY1]
			for iA2 in 1:sizeA
				for iB2 in 1:sizeB
					i2 = coord[iA2, iB2]
					phaseterm = exp(-imk * (field.x_X[iX1] * nsx_A[iA2] + field.y_Y[iY1] * nsy_B[iB2]))
					t12[i2, i1] = phaseterm * ΔxΔy
					t21[i1, i2] = ΔkxΔky / phaseterm
				end
			end
		end
	end

	if field.dir > 0
		fieldl = FieldSpace{T}(copy(field.x_X), copy(field.y_Y), field.e_SXY, field.λ, field.n, -1, copy(field.ref))
		fieldr = FieldAngularSpectrum{T}(nsx_A, nsy_B, field.e_SXY, field.λ, field.n, 1, copy(field.ref))
	else
		fieldl = FieldAngularSpectrum{T}(nsx_A, nsy_B, field.e_SXY, field.λ, field.n, -1, copy(field.ref))
		fieldr = FieldSpace{T}(copy(field.x_X), copy(field.y_Y), field.e_SXY, field.λ, field.n, 1, copy(field.ref))
		aux = t12
		t12 = t21
		t21 = aux
		aux = r12
		r12 = r21
		r21 = aux
	end

	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), typeof(r12), typeof(t12)}(r12, t12, r21, t21, fieldl, fieldr)
end

@inline function coefficient_general(fourier::FourierTransform{T}, field::FieldSpace{T,X}) where {T, X <: AbstractRange}
	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	Δx = field.x_X[2] - field.x_X[1]
	Δy = field.y_Y[2] - field.y_Y[1]
	nsx_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δx)) * field.λ;
	nsy_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δy)) * field.λ;
	return coefficient_general(fourier, field, nsx_X, nsy_Y)
end

function coefficient_general(fourier::FourierTransform{T}, field::FieldAngularSpectrum{T,X}, x_A::AbstractRange{<:Real}, y_B::AbstractRange{<:Real}) where {T, X <: AbstractRange}
	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	(sizeA, sizeB) = (length(x_A), length(y_B))
	r12 = spzeros(Complex{T}, sizeX*sizeY, sizeX * sizeY)
	r21 = spzeros(Complex{T}, sizeA*sizeB, sizeA * sizeB)
	t12 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)
	t21 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)

	k = 2π / field.λ
	imk = im * k

	Δnsx = field.nsx_X[2] - field.nsx_X[1]
	Δnsy = field.nsy_Y[2] - field.nsy_Y[1]

	# x_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δnsx)) * field.λ;
	# y_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δnsy)) * field.λ;
	Δx = x_A[2] - x_A[1]
	Δy = y_B[2] - y_B[1]
	ΔxΔy = Δx * Δy / 4π^2

	ΔkxΔky = k^2 * Δnsx * Δnsy
	coord = LinearIndices((sizeX, sizeY))
	@inbounds Threads.@threads for iX1 in 1:sizeX
		@simd for iY1 in 1:sizeY
			i1 = coord[iX1, iY1]
			for iA2 in 1:sizeA
				for iB2 in 1:sizeB
					i2 = coord[iA2, iB2]
					phaseterm = exp(imk * (field.nsx_X[iX1] * x_A[iA2] + field.nsy_Y[iY1] * y_B[iB2]))
					t12[i2, i1] = phaseterm * ΔkxΔky
					t21[i1, i2] = ΔxΔy / phaseterm
				end
			end
		end
	end

	if field.dir > 0
		fieldl = FieldAngularSpectrum{T}(copy(field.nsx_X), copy(field.nsy_Y), field.e_SXY, field.λ, field.n, -1, copy(field.ref))
		fieldr = FieldSpace{T}(x_A, y_B, field.e_SXY, field.λ, field.n, 1, copy(field.ref))
	else
		fieldl = FieldSpace{T}(x_A, y_B, field.e_SXY, field.λ, field.n, -1, copy(field.ref))
		fieldr = FieldAngularSpectrum{T}(copy(field.nsx_X), copy(field.nsy_Y), field.e_SXY, field.λ, field.n, 1, copy(field.ref))
		aux = t12
		t12 = t21
		t21 = aux
		aux = r12
		r12 = r21
		r21 = aux
	end

	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), typeof(r12), typeof(t12)}(r12, t12, r21, t21, fieldl, fieldr)
end

@inline function coefficient_general(fourier::FourierTransform{T}, field::FieldAngularSpectrum{T,X}) where {T, X <: AbstractRange}
	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	Δnsx = field.nsx_X[2] - field.nsx_X[1]
	Δnsy = field.nsy_Y[2] - field.nsy_Y[1]

	x_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δnsx)) * field.λ;
	y_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δnsy)) * field.λ;
	return coefficient_general(fourier, field, x_X, y_Y)
end

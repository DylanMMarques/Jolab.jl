struct FourierTransform{T, X <:Union{AbstractVector{T}, Nothing}} <: AbstractOpticalComponent{T}
	x_X::X
	y_Y::X
	nsx_X::X
	nsy_Y::X
end

function FourierTransform(x_X::X1, y_Y::X2, nsx_X::X3, nsy_Y::X4) where {X1 <: AbstractVector{T}, X2 <: AbstractVector{T}, X3 <: AbstractVector{T}, X4 <: AbstractVector{T}} where T <: Real
	X = promote_type(X1, X2, X3, X4)
	return FourierTransform{T, X}(x_X, y_Y, nsx_X, nsy_Y)
end

FourierTransform{T}() where T = FourierTransform{T, Nothing}(nothing, nothing, nothing, nothing)
FourierTransform() = FourierTransform{Float64}()

@inline function checkapplicability(fourier::FourierTransform{T,X}, field::FieldSpace{T}) where {T, X <: AbstractVector}
	isapprox(fourier.x_X, field.x_X, atol = @tol) || return false
	isapprox(fourier.y_Y, field.y_Y, atol = @tol) || return false
	return true
end

@inline function checkapplicability(fourier::FourierTransform{T,X}, field::FieldAngularSpectrum{T}) where {T, X <: AbstractVector}
	isapprox(fourier.nsx_X, field.nsx_X, atol = @tol) || return false
	isapprox(fourier.nsy_Y, field.nsy_Y, atol = @tol) || return false
	return true
end

@inline checkapplicability(fourier::FourierTransform{T,X}, field::Union{FieldAngularSpectrum{T}, FieldSpace{T}}) where {T, X <:Nothing} = true

function coefficient_general(fourier::FourierTransform{T,X}, field::FieldSpace{T}) where {T, X <: AbstractVector}
	checkapplicability(fourier, field) || tobedone()

	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	(sizeA, sizeB) = (length(fourier.nsx_X), length(fourier.nsy_Y))
	r12 = UniformScaling(zero(Complex{T}))
	r21 = UniformScaling(zero(Complex{T}))
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)

	k = 2π / field.λ
	imk = im * k
	coordAB = LinearIndices((sizeA, sizeB))
	coordXY = LinearIndices((sizeX, sizeY))
	@inbounds Threads.@threads for iX1 in 1:sizeX
		(xmin, xmax) = integralExtremes(field.x_X, iX1)
	 	@simd for iY1 in 1:sizeY
			(ymin, ymax) = integralExtremes(field.y_Y, iY1)
			i1 = coordXY[iX1, iY1]
			for iA2 in 1:sizeA
				for iB2 in 1:sizeB
					t12[coordAB[iA2, iB2], i1] = integrate_exp_x_y(-k * fourier.nsx_X[iA2], -k * fourier.nsy_Y[iB2], zero(T), xmin, xmax, ymin, ymax) / 4π^2
				end
			end
		end
	end
	@inbounds Threads.@threads for iA1 in 1:sizeA
		(nsxmin, nsxmax) = integralExtremes(fourier.nsx_X, iA1)
	 	for iB1 in 1:sizeB
			(nsymin, nsymax) = integralExtremes(fourier.nsy_Y, iB1)
			i1 = coordAB[iA1, iB1]
			for iX2 in 1:sizeX
				for iY2 in 1:sizeY
					t21[coordXY[iX2, iY2], i1] = integrate_exp_x_y(k * field.x_X[iX2], k * field.y_Y[iY2], zero(T), nsxmin, nsxmax, nsymin, nsymax) * k^2
				end
			end
		end
	end

	if field.dir > 0
		fieldl = FieldSpace{T}(copy(field.x_X), copy(field.y_Y), field.e_SXY, field.λ, field.n, -1, copy(field.ref))
		fieldr = FieldAngularSpectrum{T}(fourier.nsx_X, fourier.nsy_Y, zeros(Complex{T},1,sizeA,sizeB), field.λ, field.n, 1, copy(field.ref))
	else
		fieldl = FieldAngularSpectrum{T}(fourier.nsx_X, fourier.nsy_Y, zeros(Complex{T},1,sizeA,sizeB), field.λ, field.n, -1, copy(field.ref))
		fieldr = FieldSpace{T}(copy(field.x_X), copy(field.y_Y), field.e_SXY, field.λ, field.n, 1, copy(field.ref))
		aux = t12
		t12 = t21
		t21 = aux
	end
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), typeof(r12), typeof(t12)}(r12, t12, r21, t21, fieldl, fieldr)
end

@inline function coefficient_general(fourier::FourierTransform{T,X}, field::FieldSpace{T,Y}) where {T, X<:Nothing, Y <: AbstractRange}
	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	Δx = Δvector(field.x_X, 1)
	Δy = Δvector(field.y_Y, 1)
	nsx_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δx)) * field.λ;
	nsy_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δy)) * field.λ;
	fourieraux = FourierTransform{T,Y}(field.x_X, field.y_Y, nsx_X, nsy_Y)
	return coefficient_general(fourieraux, field)
end

function coefficient_general(fourier::FourierTransform{T,X}, field::FieldAngularSpectrum{T}) where {T, X <: AbstractVector}
	checkapplicability(fourier, field) || tobedone()
	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	(sizeA, sizeB) = (length(fourier.x_X), length(fourier.y_Y))
	r12 = UniformScaling(zero(Complex{T}))
	r21 = UniformScaling(zero(Complex{T}))
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)

	k = 2π / field.λ
	imk = im * k

	coordXY = LinearIndices((sizeX, sizeY))
	coordAB = LinearIndices((sizeA, sizeB))
	@inbounds Threads.@threads for iX1 in 1:sizeX
		(nsxmin, nsxmax) = integralExtremes(field.nsx_X, iX1)
		@simd for iY1 in 1:sizeY
			(nsymin, nsymax) = integralExtremes(field.nsy_Y, iY1)
			i1 = coordXY[iX1, iY1]
			for iA2 in 1:sizeA
				for iB2 in 1:sizeB
					t12[coordAB[iA2, iB2], i1] = k^2 * integrate_exp_x_y(k * fourier.x_X[iA2], k*fourier.y_Y[iB2], zero(T), nsxmin, nsxmax, nsymin, nsymax)
				end
			end
		end
	end

	@inbounds Threads.@threads for iA1 in 1:sizeA
		(xmin, xmax) = integralExtremes(fourier.x_X, iA1)
		@simd for iB1 in 1:sizeB
			(ymin, ymax) = integralExtremes(fourier.y_Y, iB1)
			i1 = coordAB[iA1, iB1]
			for iX2 in 1:sizeX
				for iY2 in 1:sizeY
					t21[coordXY[iX2, iY2], i1] = integrate_exp_x_y(-k * field.nsx_X[iX2], -k * field.nsy_Y[iY2], zero(T), xmin, xmax, ymin, ymax)
				end
			end
		end
	end

	if field.dir > 0
		fieldl = FieldAngularSpectrum{T}(copy(field.nsx_X), copy(field.nsy_Y), field.e_SXY, field.λ, field.n, -1, copy(field.ref))
		fieldr = FieldSpace{T}(fourier.x_X, fourier.y_Y, zeros(Complex{T},1,sizeA,sizeB), field.λ, field.n, 1, copy(field.ref))
	else
		fieldl = FieldSpace{T}(fourier.x_X, fourier.y_Y, zeros(Complex{T}, 1,sizeA,sizeB), field.λ, field.n, -1, copy(field.ref))
		fieldr = FieldAngularSpectrum{T}(copy(field.nsx_X), copy(field.nsy_Y), field.e_SXY, field.λ, field.n, 1, copy(field.ref))
		aux = t12
		t12 = t21
		t21 = aux
	end

	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), typeof(r12), typeof(t12)}(r12, t12, r21, t21, fieldl, fieldr)
end

@inline function coefficient_general(fourier::FourierTransform{T,X}, field::FieldAngularSpectrum{T,Y}) where {T, X <: Nothing, Y <: AbstractRange}
	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	Δnsx = Δvector(field.nsx_X, 1)
	Δnsy = Δvector(field.nsy_Y, 1)
	x_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δnsx)) * field.λ;
	y_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δnsy)) * field.λ;
	fourieraux = FourierTransform{T,Y}(x_X, y_Y, field.nsx_X, field.nsy_Y)
	return coefficient_general(fourieraux, field)
end

function coefficient_specific(fourrier::FourierTransform{T,X}, field::F) where {T, X <: Nothing, F <: FieldSpace{T, <:AbstractRange}}
	planfft = FFTW.plan_fft(view(field.e_SXY,1,:,:))
	(sizeX, sizeY) = (length(field.x_X), length(field.y_Y))

	(Δx, Δy) = (Δvector(field.x_X, 1), Δvector(field.y_Y, 1))

	nsx_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δx)) * field.λ;
	nsy_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δy)) * field.λ;
	scalar = Δx * Δy / 4 / π^2

	(Δnsx, Δnsy) = (Δvector(nsx_X, 1), Δvector(nsy_Y, 1))
	i_scalar = sizeX * sizeY * Δnsx * Δnsy * 4 * π^2 / field.λ^2

	tmp = Array{Complex{T},2}(undef, sizeX, sizeY)
	if field.dir > 0
		fieldr = FieldAngularSpectrum{T}(nsx_X, nsy_Y, field.e_SXY, field.λ, field.n, 1, field.ref)
		return FFTCoefficient{T,F,typeof(fieldr),typeof(planfft)}(planfft, scalar, i_scalar, tmp, field, fieldr)
	else
		fieldl = FieldAngularSpectrum{T}(nsx_X, nsy_Y, field.e_SXY, field.λ, field.n, -1, field.ref)
		return FFTCoefficient{T,typeof(fieldl),F,typeof(planfft)}(planfft, scalar, i_scalar, tmp, fieldl, field)
	end
end

function coefficient_specific(fourrier::FourierTransform{T,X}, field::F) where {T, X <: Nothing, F <: FieldAngularSpectrum{T, <:AbstractRange}}
	planfft = FFTW.plan_fft(view(field.e_SXY,1,:,:))
	(sizeX, sizeY) = (length(field.nsx_X), length(field.nsy_Y))

	(Δnsx, Δnsy) = (Δvector(field.nsx_X, 1), Δvector(field.nsy_Y, 1))

	x_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δnsx)) * field.λ;
	y_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δnsy)) * field.λ;
	(Δx, Δy) = (Δvector(x_X, 1), Δvector(y_Y, 1))
	scalar = Δx * Δy / 4 / π^2

	i_scalar = sizeX * sizeY * Δnsx * Δnsy * 4 * π^2 / field.λ^2

	tmp = Array{Complex{T},2}(undef, sizeX, sizeY)
	if field.dir > 0
		fieldr = FieldSpace{T}(x_X, y_Y, field.e_SXY, field.λ, field.n, 1, field.ref)
		return FFTCoefficient{T,F,typeof(fieldr),typeof(planfft)}(planfft, scalar, i_scalar, tmp, field, fieldr)
	else
		fieldl = FieldSpace{T}(x_X, y_Y, field.e_SXY, field.λ, field.n, -1, field.ref)
		return FFTCoefficient{T,typeof(fieldl),F,typeof(planfft)}(planfft, scalar, i_scalar, tmp, fieldl, field)
	end
end

struct FFTCoefficient{T,L,R, F <: AbstractFFTs.Plan} <: AbstractCoefficient{T,L,R}
	planfft::F
	scalar::T
	i_scalar::T
	tmp::Array{Complex{T}, 2}
	fieldl::L
	fieldr::R
end

function lightinteraction!(fieldl::L, fieldr::R, coef::FFTCoefficient{T,L,R}, fieldi::Union{L,R}) where {T,L,R}
	if fieldi.dir > 0
		samedefinitions(fieldl, fieldi) || tobedone()
		if fieldi isa FieldSpace
			mul!(coef.tmp, coef.planfft, fftshift(view(fieldi.e_SXY,1,:,:)))
			vec(fieldr.e_SXY) .= vec(fftshift(coef.tmp)) * coef.scalar
		else
			ldiv!(coef.tmp, coef.planfft, fftshift(view(fieldi.e_SXY,1,:,:)))
			vec(fieldr.e_SXY) .= vec(fftshift(coef.tmp)) * coef.i_scalar
		end
		fieldl.e_SXY .= zero(Complex{T})
	else
		samedefinitions(fieldr, fieldi) || tobedone()
		if fieldi isa FieldSpace
			mul!(coef.tmp, coef.planfft, fftshift(view(fieldi.e_SXY,1,:,:)))
			vec(fieldl.e_SXY) .= vec(fftshift(coef.tmp)) * coef.scalar
		else
			ldiv!(coef.tmp, coef.planfft, fftshift(view(fieldi.e_SXY,1,:,:)))
			vec(fieldl.e_SXY) .= vec(fftshift(coef.tmp)) * coef.i_scalar
		end
		fieldr.e_SXY .= zero(Complex{T})
	end
end

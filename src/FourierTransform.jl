struct FourierTransform{T, X <:Union{AbstractVector{T}, Nothing}, Y<:Union{ReferenceFrame{T}, Nothing}} <: AbstractOpticalComponent{T}
	x_X::X
	y_Y::X
	nsx_X::X
	nsy_Y::X
	ref::Y
end

function FourierTransform(x_X::X1, y_Y::X2, nsx_X::X3, nsy_Y::X4) where {X1 <: AbstractVector{T}, X2 <: AbstractVector{T}, X3 <: AbstractVector{T}, X4 <: AbstractVector{T}} where T <: Real
	X = promote_type(X1, X2, X3, X4)
	return FourierTransform{T, X, Nothing}(x_X, y_Y, nsx_X, nsy_Y, nothing)
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

@inline function fourierΔintegral(fourier::FourierTransform{T,X,Y}, field::FieldSpace, iX, iY, iA, iB) where {T, X, Y<:Nothing}
	(xmin, xmax) = integralExtremes(field.x_X, iX)
	(ymin, ymax) = integralExtremes(field.y_Y, iY)

	(nsxmin, nsxmax) = integralExtremes(fourier.nsx_X, iA)
	(nsymin, nsymax) = integralExtremes(fourier.nsy_Y, iB)

	k = 2π / field.λ

	return (integrate_exp_x_y(-k * fourier.nsx_X[iA], -k * fourier.nsy_Y[iB], zero(T), xmin, xmax, ymin, ymax) / 4π^2, integrate_exp_x_y(k * field.x_X[iX], k * field.y_Y[iY], zero(T), nsxmin, nsxmax, nsymin, nsymax) * k^2)
end

@inline function fourierΔintegral(fourier::FourierTransform{T,X,Y}, field::FieldAngularSpectrum, iX, iY, iA, iB) where {T, X, Y<:Nothing}
	(nsxmin, nsxmax) = integralExtremes(field.nsx_X, iX)
	(nsymin, nsymax) = integralExtremes(field.nsy_Y, iY)

	(xmin, xmax) = integralExtremes(fourier.x_X, iA)
	(ymin, ymax) = integralExtremes(fourier.y_Y, iB)

	k = 2π / field.λ

	return (k^2 * integrate_exp_x_y(k * fourier.x_X[iA], k*fourier.y_Y[iB], zero(T), nsxmin, nsxmax, nsymin, nsymax), integrate_exp_x_y(-k * field.nsx_X[iX], -k * field.nsy_Y[iY], zero(T), xmin, xmax, ymin, ymax) / 4π^2)
end

@inline function fourierΔintegral(fourier::FourierTransform{T,X,Y}, field::FieldAngularSpectrum, iX, iY, iA, iB) where {T, X, Y<:ReferenceFrame}
	(nsxmin, nsxmax) = integralExtremes(field.nsx_X, iX)
	(nsymin, nsymax) = integralExtremes(field.nsy_Y, iY)

	(xmin, xmax) = integralExtremes(fourier.x_X, iA)
	(ymin, ymax) = integralExtremes(fourier.y_Y, iB)

	k = 2π / field.λ

	(refΔx, refΔy, refΔz) = (fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z)
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ);

	imagina_waves = (imag(field.n) > @tol) || nsxmin^2 + nsymin^2 > real(field.n)^2 || nsxmin^2 + nsymax^2 > real(field.n)^2 || nsxmax^2 + nsymin^2 > real(field.n)^2 || nsxmax^2 + nsymax^2 > real(field.n)^2

	if imagina_waves
		@inline f(nsr) = k * ((fourier.x_X[iA] + refΔx) * nsr[1] + (fourier.y_Y[iB] + refΔy) * nsr[2] + √(complex(field.n^2 - nsr[1]^2 - nsr[2]^2)) * refΔz)
		t12 = k^2 * hcubature(f, SVector(nsxmin, nsymin), SVector(nsxmax, nsymax))
	else
		@inline g(nsx, nsy) = k * ((fourier.x_X[iA] + refΔx) * nsx + (fourier.y_Y[iB] + refΔy) * nsy + √(real(field.n)^2 - nsx^2 - nsy^2) * refΔz)
		t12 = k^2 * integrate_exp_xy_x_y(g, nsxmin, nsxmax, nsymin, nsymax)
	end

	aux = exp(-im * k * refΔz * √(field.n^2 - fourier.nsx_X[iX]^2 - fourier.nsy_Y[iY]^2))

	return (t12, aux * integrate_exp_x_y(-k * field.nsx_X[iX], -k * field.nsy_Y[iY], zero(T), xmin, xmax, ymin, ymax) / 4π^2)
end

@inline function fourierΔintegral(fourier::FourierTransform{T,X,Y}, field::FieldSpace, iX, iY, iA, iB) where {T, X, Y<:ReferenceFrame}
	(xmin, xmax) = integralExtremes(field.x_X, iX)
	(ymin, ymax) = integralExtremes(field.y_Y, iY)

	(nsxmin, nsxmax) = integralExtremes(fourier.nsx_X, iA)
	(nsymin, nsymax) = integralExtremes(fourier.nsy_Y, iB)

	k = 2π / field.λ

	(refΔx, refΔy, refΔz) = (fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z)
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ);

	imagina_waves = (imag(field.n) > @tol) || nsxmin^2 + nsymin^2 > real(field.n)^2 || nsxmin^2 + nsymax^2 > real(field.n)^2 || nsxmax^2 + nsymin^2 > real(field.n)^2 || nsxmax^2 + nsymax^2 > real(field.n)^2
	if imagina_waves
		f(nsr) = k * ((fourier.x_X[iX] + refΔx) * nsr[1] + (fourier.y_Y[iY] + refΔy) * nsr[2] - √(complex(field.n^2 - nsr[1]^2 - nsr[2]^2)) * refΔz)
		t21 = k^2 * hcubature(f, SVector(nsxmin, nsymin), SVector(nsxmax, nsymax))
	else
		g(nsx, nsy) = k * ((fourier.x_X[iA] + refΔx) * nsx + (fourier.y_Y[iB] + refΔy) * nsy - √(real(field.n)^2 - nsx^2 - nsy^2) * refΔz)
		t21 = k^2 * integrate_exp_xy_x_y(g, nsxmin, nsxmax, nsymin, nsymax)
	end

	aux = exp(im * k * refΔz * √(field.n^2 - fourier.nsx_X[iA]^2 - fourier.nsy_Y[iB]^2))

	return (aux * integrate_exp_x_y(-k * field.nsx_X[iA], -k * field.nsy_Y[iB], zero(T), xmin, xmax, ymin, ymax) / 4π^2, t21)
end


function coefficient_general(fourier::FourierTransform{T,X,Y}, field::FieldSpace{T}) where {T, X <: AbstractVector, Y<:Nothing}
	checkapplicability(fourier, field) || tobedone()

	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	(sizeA, sizeB) = (length(fourier.nsx_X), length(fourier.nsy_Y))
	r12 = UniformScaling(zero(Complex{T}))
	r21 = UniformScaling(zero(Complex{T}))
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)

	coordAB = LinearIndices((sizeA, sizeB))
	coordXY = LinearIndices((sizeX, sizeY))
	@inbounds Threads.@threads for iX1 in 1:sizeX
	 	@simd for iY1 in 1:sizeY
			i1 = coordXY[iX1, iY1]
			for iA2 in 1:sizeA
				for iB2 in 1:sizeB
					i2 = coordAB[iA2, iB2]
					(t12[i2, i1], t21[i1, i2]) = fourierΔintegral(fourier, field, iX1, iY1, iA2, iB2)
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

	(r12, r21) = (UniformScaling(zero(Complex{T})), UniformScaling(zero(Complex{T})))
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)

	coordXY = LinearIndices((sizeX, sizeY))
	coordAB = LinearIndices((sizeA, sizeB))
	@inbounds Threads.@threads for iX1 in 1:sizeX
		@simd for iY1 in 1:sizeY
			i1 = coordXY[iX1, iY1]
			for iA2 in 1:sizeA
				for iB2 in 1:sizeB
					i2 = coordAB[iA2, iB2]
					(t12[i2, i1], t21[i1, i2]) = fourierΔintegral(fourier, field, iX1, iY1, iA2, iB2)
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

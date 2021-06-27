function fourriertransformfft(x_X::AbstractRange{<:Real}, y_Y::AbstractRange{<:Real}, z_XY::AbstractArray{T,2}, λ::Real; padding = 0::Integer) where {T<:Number}
	(sizeX, sizeY) = size(z_XY)

	if padding > 0
		tmp_z_XY = zeros(T, sizeX + 2 * padding, sizeY + 2 * padding)
		@inbounds @simd for iY in 1: sizeY
			for iX in 1:sizeX
				tmp_z_XY[iX+padding,iY+padding] = z_XY[iX,iY];
			end
		end
	else
		tmp_z_XY = z_XY;
	end
	dX = x_X[2] - x_X[1]
	dY = y_Y[2] - y_Y[1]
	nsx = fftshift(FFTW.fftfreq(sizeX + 2 * padding, 1 / dX)) * λ;
	nsy = fftshift(FFTW.fftfreq(sizeY + 2 * padding, 1 / dY)) * λ;
	if padding > 0
		nsx = nsx[padding+1:padding+sizeX];
		nsy = nsy[padding+1:padding+sizeX];
	end
	# Normalization factors makes the fft with the same amplitude as the analytical fourier transform
	tmp_fz_XY = fftshift(fft(fftshift(tmp_z_XY))) * (dX * dY / 4 / π^2);
	if padding > 0
		fz_XY = tmp_fz_XY[padding+1:padding+sizeX, padding+1:padding+sizeY]
	else
		fz_XY = tmp_fz_XY;
	end
	return (nsx, nsy, fz_XY);
end

function fourriertransformfft(x_X::AbstractRange{<:Real}, y_Y::AbstractRange{<:Real}, z_SXY::AbstractArray{T,3}, λ::Real; padding = 0::Integer) where {T<:Number}
	(sizeS, sizeX, sizeY) = size(z_SXY);

	dX = x_X[2] - x_X[1]
	dY = y_Y[2] - y_Y[1]

	nsx = fftshift(FFTW.fftfreq(sizeX + 2 * padding, 1 / dX)) * λ;
	nsy = fftshift(FFTW.fftfreq(sizeY + 2 * padding, 1 / dY)) * λ;
	if padding > 0
		nsx = nsx[padding+1:padding+sizeX];
		nsy = nsy[padding+1:padding+sizeY];
	end
	cons = (dX * dY / 4 / π^2);

	if padding > 0
		fz_SXY = Array{T,3}(undef, sizeS, sizeX, sizeY);
		tmp_z_XY = zeros(T, sizeX + 2 * padding, sizeY + 2 * padding)
	else
		tmp_z_XY = @view z_SXY[1,:,:];
	end

	fftOperator = plan_fft(tmp_z_XY);

	@inbounds for iS in 1:sizeS
		if padding > 0
			for iY in 1:sizeY
				for iX in 1:sizeX
					tmp_z_XY[iX+padding,iY+padding] = z_SXY[iS,iX,iY];
				end
			end
		else
			tmp_z_XY = @view z_SXY[iS,:,:];
		end
		tmp_fz_XY = fftshift(fftOperator * fftshift(tmp_z_XY));
		tmp_fz_XY .*= cons;
		if padding > 0
			fz_SXY[iS,:,:] .= tmp_fz_XY[padding+1:padding+sizeX, padding+1:padding+sizeY]
		else
			fz_SXY = reshape(tmp_fz_XY, sizeS, sizeX, sizeY)
		end
	end
	return (nsx, nsy, fz_SXY);
end

function inversefourriertransformfft(nsx_X::AbstractRange, dnsy_Y::AbstractRange, fz_XY::AbstractArray{T,2}, λ::Real; padding = 0::Integer)::Tuple{AbstractVector{<:Number}, AbstractVector{<:Number}, AbstractArray{<:Number, 2}} where {T<:Number}
	(sizeX, sizeY) = size(fz_XY);

	dnsx = nsx_X[2] - nsx_X[1]
	dnsy = nsy_Y[2] - nsy_Y[1]

	if padding > 0
		tmp_fz_XY = zeros(T, sizeX + 2 * padding, sizeY + 2 * padding)
		@inbounds for iY in 1: sizeY
			for iX in 1:sizeX
				tmp_fz_XY[iX+padding,iY+padding] = fz_XY[iX,iY];
			end
		end
	else
		tmp_fz_XY = fz_XY;
	end

	x = fftshift(FFTW.fftfreq(sizeX + 2 * padding, 1 / dnsx)) * λ
	y = fftshift(FFTW.fftfreq(sizeY + 2 * padding, 1 / dnsy)) * λ
	if padding > 0
		x = x[padding+1:padding+sizeX]
		y = y[padding+1:padding+sizeY]
	end
	tmp_z_XY = fftshift(ifft(fftshift(tmp_fz_XY))) .* ((sizeX + 2 * padding) * (sizeY + 2 * padding) * dnsx * dnsy * 4 * π^2 / λ^2); # Normalization factors makes the fft with the same amplitude as the analytical fourier transform
	if padding > 0
		z_XY = tmp_z_XY[padding+1:padding+sizeX, padding+1:padding+sizeY]
	else
		z_XY = tmp_z_XY;
	end

	return (x, y, z_XY);
end

function inversefourriertransformfft(nsx_X::AbstractRange{<:Real}, nsy_Y::AbstractRange{<:Real}, fz_SXY::AbstractArray{Complex{T}, 3}, λ::Real; padding = 0::Integer) where {T<:Real}
	(sizeS, sizeX, sizeY) = size(fz_SXY);

	dnsx = nsx_X[2] - nsx_X[1]
	dnsy = nsy_Y[2] - nsy_Y[1]

	x = fftshift(FFTW.fftfreq(sizeX + 2 * padding, 1 / dnsx)) * λ;
	y = fftshift(FFTW.fftfreq(sizeY + 2 * padding, 1 / dnsy)) * λ;
	if padding > 0
		x = x[padding+1:padding+sizeX]
		y = y[padding+1:padding+sizeY]
	end

	z_SXY = Array{Complex{T},3}(undef, sizeS, sizeX, sizeY);
	if padding > 0
		tmp_fz_XY = zeros(Complex{T}, sizeX + 2 * padding, sizeY + 2 * padding)
	else
		tmp_fz_XY = @view fz_SXY[1,:,:];
	end

	fftOperator = plan_ifft(tmp_fz_XY);
	cons = ((sizeX + 2 * padding) * (sizeY + 2 * padding) * dnsy * dnsx * 4 * π^2/ λ^2);
	@inbounds for iS in 1:sizeS
		if padding > 0
			for iY in 1:sizeY
				for iX in 1:sizeX
					tmp_fz_XY[iX+padding,iY+padding] = fz_SXY[iS,iX,iY];
				end
			end
		else
			tmp_fz_XY = @view fz_SXY[iS,:,:];
		end
		tmp_z_XY = fftshift(fftOperator * fftshift(tmp_fz_XY)) * cons;
		if padding > 0
			z_SXY[iS,:,:] .= tmp_z_XY[padding+1:sizeX+padding, padding+1:sizeY+padding]
		else
			view(z_SXY,iS,:,:) .= tmp_z_XY;
		end
	end
	return (x, y, z_SXY);
end

function ΔIntegrationTrap(x_X::AbstractVector{T}) where T<:Number
	sizeX = length(x_X);
	Δx_X = Vector{T}(undef, sizeX);
	Δx_X[1] = (x_X[2] - x_X[1]) / 2;
	@inbounds for iX in 2:sizeX-1
		Δx_X[iX] = (x_X[iX+1] - x_X[iX-1]) / 2;
	end
	Δx_X[sizeX] = (x_X[sizeX] - x_X[sizeX-1]) / 2;
	return Δx_X;
end

function dotdim(a, b, dim)
	return real.(dropdims(sum(a .* conj(b), dims = dim), dims = dim));
end

function ∫(f_X::AbstractVector{T}, x_X::AbstractVector{<:Number}) where T <: Number
	Δx_X = ΔIntegrationTrap(x_X);
	sizeX = length(x_X);
	res = convert(T, 0);
	sizeX == length(f_X) || error("Size of vectors for integration must be the same")
	@inbounds for iX in 1:sizeX
		res += f_X[iX] * Δx_X[iX];
	end
	return res;
end

function ∫∫(f_XY::AbstractArray{T,2}, x_X::AbstractVector{<:Number}, y_Y::AbstractVector{<:Number}) where T <: Number
	Δx_X = ΔIntegrationTrap(x_X);
	Δy_Y = ΔIntegrationTrap(y_Y);
	sizeX = length(x_X);
	sizeY = length(y_Y);
	res = zero(T);
	(sizeX == size(f_XY, 1) && sizeY == size(f_XY, 2)) || error("Size of vector for integration must be the same");

	@inbounds @simd for iX in 1:sizeX
		aux = zero(T);
		@simd for iY in 1:sizeY
			aux += f_XY[iX,iY] * Δy_Y[iY];
		end
		res += aux * Δx_X[iX];
	end
	return res;
end

function ∫∫(f_XY::AbstractArray{T,2}, x_X::AbstractRange{<:Number}, y_Y::AbstractRange{<:Number}) where T <: Number
	(length(x_X) == size(f_XY, 1) && length(y_Y) == size(f_XY, 2)) || error("Size of vector for integration must be the same");
	res = zero(T)
	@inbounds @simd for i in eachindex(f_XY)
		res += f_XY[i]
	end
	return res * (x_X[2] - x_X[1]) * (y_Y[2] - y_Y[1]);
end

function anglebetweenvectors(v1::Vector{<:Number}, v2::Vector{<:Number})::Number
	return acos(v1 ⋅ v2 ./ norm(v1) ./ norm(v2));
end

function checkrange(x_X::AbstractVector{<:Number})
	Δx = x_X[2] - x_X[1]
	@inbounds @simd for i in 2:length(x_X)
		(abs(Δx - (x_X[i] - x_X[i-1])) < @tol) || return false
	end
	return true
end

function rextrema(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
	min = x[1]^2 + y[1]^2
	max = x[1]^2 + y[1]^2
	@inbounds @simd for xi in x
		for yi in y
			r = xi^2 + yi^2
			(min > r) && (min = r)
			(max < r) && (max = r)
		end
	end
	return (√min, √max)
end

@inline function Δvector(x_X::Vector{T}, i::Integer) where T
	sizeX = length(x_X)
	if 1 < i < sizeX
		return (x_X[i + 1] - x_X[i - 1]) / 2
	elseif i == 1
		return (x_X[2] - x_X[1]) / 2
	elseif i == sizeX
		return (x_X[sizeX] - x_X[sizeX - 1]) / 2
	else
		return zero(T)
	end
end

@inline function Δvector(x_X::AbstractRange{<:Number}, i::Integer)
	return x_X.step.hi
end

@inline function integralExtremes(x_X::AbstractVector{<:Number}, i::Integer)
	sizeX = length(x_X)
	if i == 1
		return (x_X[1], x_X[1] + (x_X[2] - x_X[1]) / 2)
	elseif i >= sizeX
		i == sizeX && return (x_X[sizeX] - (x_X[sizeX] - x_X[sizeX-1]) / 2, x_X[sizeX])
		error("Programming error")
	else
		return (x_X[i] - (x_X[i] - x_X[i-1]) / 2, x_X[i] + (x_X[i+1] - x_X[i]) / 2)
	end
end

@inline function integralExtremes(x_X::AbstractRange{T}, i::Integer)::Tuple{T,T} where T
	sizeX = length(x_X)
	if i == 1
		return  (x_X[1], x_X[1] + x_X.step.hi / 2)
	elseif i >= sizeX
		i == sizeX && return (x_X[i] - x_X.step.hi / 2, x_X[i])
		error("Programming error")
	else
		return (x_X[i] - x_X.step.hi / 2, x_X[i] + x_X.step.hi / 2)
	end
end

@inline function integrate_exp_x_y(β::T, γ::T, δ::T, xmin::T, xmax::T, ymin::T, ymax::T)::Complex{T} where T
	if β == 0 && γ == 0
		return exp(im * δ) * (xmax - xmin) * (ymax - ymin)
	elseif γ == 0
		return im * exp(im * δ) * (-exp(im * β * xmax) + exp(im * β * xmin)) * (ymax - ymin) / β
	elseif β == 0
		return im * exp(im * δ) * (-exp(im * γ * ymax) + exp(im * γ * ymin)) * (xmax - xmin) / γ
	else
		return exp(im * δ) * (exp(im * β * xmax) - exp(im * β * xmin)) * (-exp(im * γ * ymax) + exp(im * γ * ymin)) / (γ * β)
	end
end

@inline function integrate_exp_xy_x_y(α::T, β::T, γ::T, δ::T, xmin::T, xmax::T, ymin::T, ymax::T)::Complex{T} where T
	# g(x,y) = (α * x * y + β * x + γ * y + δ)
	if abs(α) <= 1E-7
		return integrate_exp_x_y(β, γ, δ, xmin, xmax, ymin, ymax)
	else
		if (0 <= (γ + α * xmin) / (α * xmin - α * xmax) <= 1) && (0 <= (β + α * ymin) / (α * ymin - α * ymax) <= 1)
			b(x,y) = exp(im * (α * x * y + β * x + γ * y + δ))
			return (xmax - xmin) * (ymax - ymin) * (b(xmax,ymax) + b(xmax, ymin) + b(xmin, ymax) + b(xmin, ymin)) / 4
		else
			@inline Pf(x::T, y::T)::Complex{T} = begin
				aux = (γ + α * x) * (β + α * y) / α
				if aux < 0
					return ((-conj(Complex{T}(expint(im * aux)))) - π * im)
				else
					return ((-conj(Complex{T}(expint(im * aux)))) + π * im)
				end
 			end
			return - im / α * exp(im * (δ - β * γ / α)) * (Pf(xmax, ymax) - Pf(xmax, ymin) - Pf(xmin, ymax) + Pf(xmin, ymin))
		end
	end
	return zero(Complex{T})
end

@inline function integrate_exp_xy_x_y(f::Function, xmin::T, xmax::T, ymin::T, ymax::T)::Complex{T} where T
	f11, f12, f21, f22 = f(xmin, ymin), f(xmin, ymax), f(xmax, ymin), f(xmax, ymax)
	(α, β, γ, δ) = bilinearinterpolation(f11, f12, f21, f22, xmin, xmax, ymin, ymax)
	return integrate_exp_xy_x_y(α, β, γ, δ, xmin, xmax, ymin, ymax)
end

function bilinearinterpolation(f11::T, f12, f21, f22, xmin, xmax, ymin, ymax)::Tuple{T,T,T,T} where T
	# fits the 4 data points to a xy + b x + c y + d
	a = (f11 - f12 - f21 + f22) / (xmax - xmin) / (ymax - ymin)
	b = (-f11 * ymax + f21 * ymax + f12 * ymin - f22 * ymin) / (xmax - xmin) / (ymax - ymin)
	c = (-f11 * xmax + f21 * xmin + f12 * xmax - f22 * xmin) / (xmax - xmin) / (ymax - ymin)
	d = (f11 * xmax * ymax - f21 * xmin * ymax - f12 * xmax * ymin + f22 * xmin * ymin) / (xmax - xmin) / (ymax - ymin)
	return (a,b,c,d)
end

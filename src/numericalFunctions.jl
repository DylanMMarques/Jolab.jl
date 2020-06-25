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

function fourriertransform!(o_SAB::AbstractArray{Complex{T}, 3}, x_X::AbstractVector{<:Real}, y_Y::AbstractVector{<:Real}, z_SXY::AbstractArray{<:Number, 3}, λ::Real, n::Number, nsx_A::AbstractVector{<:Real}, nsy_B::AbstractVector{<:Real}) where {T<:Real}
	sizeS, sizeX, sizeY, sizeA, sizeB = size(z_SXY, 1), length(x_X), length(y_Y), length(nsx_A), length(nsy_B);

	sizeX == size(z_SXY, 2) && sizeY == size(z_SXY, 3) || error("Wrong sizes");
	sizeA == size(o_SAB, 2) && sizeB == size(o_SAB, 3) || error("Wrong sizes");
	sizeS == size(o_SAB, 1) || error("Wrong sizes")

	k = 2π / λ;
	Δx = ΔIntegrationTrap(x_X);
	Δy = 1 / (4π^2) .* ΔIntegrationTrap(y_Y);
	o_SAB .= zero(Complex{T})

	@inbounds for iX in 1:sizeX
		for iY in 1:sizeY
			ΔxΔy = Δx[iX] * Δy[iY]
			@simd for iA in 1:sizeA
				for iB in 1:sizeB
					phaseterm = exp(-im * k * (x_X[iX] * nsx_A[iA] + y_Y[iY] * nsy_B[iB])) * ΔxΔy
					for iS in 1:sizeS
						o_SAB[iS,iA,iB] += z_SXY[iS,iX,iY] * phaseterm
					end
				end
			end
		end
	end
end

@inline function fourriertransform(x_X::AbstractVector{<:Real}, y_Y::AbstractVector{<:Real}, z_SXY::AbstractArray{T, 3}, λ::Real, n::Number, nsx_A::AbstractVector{<:Real}, nsy_B::AbstractVector{<:Real}) where T<:Number
	o_SAB = zeros(T, size(z_SXY, 1), length(nsx_A), length(nsy_B));
	fourriertransform!(o_SAB, x_X, y_Y, z_SXY, λ, n, nsx_A, nsy_B);
	return o_SAB
end

function inversefourriertransform!(o_SAB::AbstractArray{Complex{T}, 3}, nsx_X::AbstractVector{<:Real}, nsy_Y::AbstractVector{<:Real}, z_SXY::AbstractArray{<:Number, 3}, λ::Real, n::Number, x_A::AbstractVector{<:Real}, y_B::AbstractVector{<:Real}) where {T<:Real}
	sizeS, sizeX, sizeY, sizeA, sizeB = size(z_SXY, 1), length(nsx_X), length(nsy_Y), length(x_A), length(y_B);

	sizeX == size(z_SXY, 2) && sizeY == size(z_SXY, 3) || error("Wrong sizes");
	sizeA == size(o_SAB, 2) && sizeB == size(o_SAB, 3) || error("Wrong sizes");
	sizeS == size(o_SAB, 1) || error("Wrong sizes")

	k = 2π / λ;
	Δkx = ΔIntegrationTrap(nsx_X) .* k
	Δky = ΔIntegrationTrap(nsy_Y) .* k
	o_SAB .= zero(Complex{T})

	@inbounds for iX in 1:sizeX
		for iY in 1:sizeY
			ΔkxΔky = Δkx[iX] * Δky[iY]
			@simd for iA in 1:sizeA
				for iB in 1:sizeB
					phaseterm = exp(im * k * (x_A[iA] * nsx_X[iX] + y_B[iB] * nsy_Y[iY])) * ΔkxΔky
					for iS in 1:sizeS
						o_SAB[iS,iA,iB] += z_SXY[iS,iX,iY] * phaseterm
					end
				end
			end
		end
	end
end

@inline function inversefourriertransform(nsx_X::AbstractVector{<:Real}, nsy_Y::AbstractVector{<:Real}, z_SXY::AbstractArray{<:Number, 3}, λ::Real, n::Number, x_A::AbstractVector{<:Real}, y_B::AbstractVector{<:Real})
	o_SAB = zeros(eltype(z_SXY), size(z_SXY, 1), length(x_A), length(y_B));
	inversefourriertransform!(o_SAB, nsx_X, nsy_Y, z_SXY, λ, n, x_A, y_B);
	return o_SAB
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

function ndgrid(v1::AbstractVector{<:Number}, v2::AbstractVector{<:Number})
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repeat(v1, 1, n), repeat(v2, m, 1))
end

function adddims!(v1::Array{T,M}, dims::NTuple{N,Int}) where {T, N, M}
	siz = collect(size(v1));
	dims = collect(dims);

	for i in dims
		insert!(siz, i, 1);
	end
	v1 = reshape(v1, Tuple(siz));
end

function adddims(v1::Array{T,M}, dims::NTuple{N,Int}) where {T, N, M}
	siz = collect(size(v1));
	dims = collect(dims);

	for i in dims
		insert!(siz, i, 1);
	end
	addims = reshape(v1, Tuple(siz));
	return addims;
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

@inline function intensity(e_SXY::AbstractArray{Complex{T}}) where {T<:Number}
	val = zero(T)
	@inbounds @simd for i in eachindex(e_SXY)
		val += abs2(e_SXY[i])
	end
	return val
end

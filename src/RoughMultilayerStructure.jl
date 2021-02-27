struct RoughMultilayerStructure{T<:Real} <: AbstractOpticalComponent{T}
	n::Vector{Material{T}}
	h::Vector{T}
	Δz::Vector{JolabFunction2D{T}}
	ref::ReferenceFrame{T}
	function RoughMultilayerStructure{T}(n, h, Δz, ref) where T
		length(n) == length(h) + 2 == length(Δz) + 1 || error("wrong sizes")
		return new{T}(n, h, Δz, ref)
	end
end
RoughMultilayerStructure(n1, n2, Δz, ref) = RoughMultilayerStructure{Float64}(n1, n2, Δz, ref)

@inline function r12_i(n_M, h_M, nsr::T, λ, iM) where T
	sizeM = length(n_M)

	k = 2π / λ
	(rm_m1, tm_m1) =  rtss₁₂(nsr, view(n_M, iM:sizeM), view(h_M, iM:sizeM-2),λ)
	if iM > 1
		(tmp, t0_m) =  rtss₁₂(nsr, view(n_M, 1:iM), view(h_M, 1:iM-2),λ)
		nsz_m = √(n_M[iM]^2 - nsr^2)
		phaseTerm_m = exp(im * k * nsz_m * h_M[iM-1])
		t0_m *= phaseTerm_m
	else
		t0_m = one(Complex{T})
	end
	return -im * k / 2 * (n_M[iM+1]^2 - n_M[iM]^2) * t0_m * (1 + rm_m1)
end

@inline function r12_s(n_M, h_M, nsr::T, λ, iM) where T
	sizeM = length(n_M)
	k = 2π / λ
	nsz_m = √(n_M[iM]^2 - nsr^2)
	(rm_m1, tm_m1) =  rtss₁₂(nsr, view(n_M, iM:sizeM), view(h_M, iM:sizeM-2),λ)
	if iM > 1
		(rm_ml1, tm_0) =  rtss₁₂(nsr, view(n_M, iM:-1:1), view(h_M, iM-2:-1:1),λ)
		phaseTerm_m = exp(im * k * nsz_m * h_M[iM-1])
		tm_0 *= phaseTerm_m
		M_m = 1 - rm_ml1 * rm_m1 * phaseTerm_m^2
	else
		t0_m = one(Complex{T})
		tm_0 = one(Complex{T})
		M_m = one(Complex{T})
	end
	return tm_0 * (1 + rm_m1) / nsz_m / M_m
end

@inline function t12_i(n_M, h_M, nsr::T, λ, iM) where T
	sizeM = length(n_M)

	k = 2π / λ
	(rm_m1, tm_m1) =  rtss₁₂(nsr, view(n_M, iM:sizeM), view(h_M, iM:sizeM-2),λ)
	if iM > 1
		nsz_m = √(n_M[iM]^2 - nsr^2)
		(tmp, t0_m) =  rtss₁₂(nsr, view(n_M, 1:iM), view(h_M, 1:iM-2),λ)
		phaseTerm_m = exp(im * k * nsz_m * h_M[iM-1])
		t0_m *= phaseTerm_m
	else
		t0_m = one(Complex{T})
	end
	return -im * k / 2 * (n_M[iM+1]^2 - n_M[iM]^2) * t0_m * (1 + rm_m1)
end

@inline function t12_s(n_M, h_M, nsr::T, λ, iM) where T
	sizeM = length(n_M)

	k = 2π / λ
	nsz_m1 = √(n_M[iM+1]^2 - nsr^2)

	(rm1_0, tm1_0) = rtss₁₂(nsr, view(n_M, iM+1:-1:1), view(h_M, iM-1:-1:1),λ)
	if iM < sizeM - 1
		(tmp, tm1_M) = rtss₁₂(nsr, view(n_M, iM+1:sizeM), view(h_M, iM+1:sizeM-2),λ)
		tm1_M *= exp(im * k * nsz_m1 * h_M[iM])
	else
		tm1_M = one(Complex{T})
	end
	return tm1_M * (1 + rm1_0) / nsz_m1
end

function coefficient_specific(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrum{T,D,X}) where {T,D, X<:AbstractRange}
	checkapplicability(rmls, fieldi)
	(sizeX, sizeY, sizeM) = (length(fieldi.nsx_X), length(fieldi.nsy_Y), length(rmls.n))

	k = 2π / fieldi.λ
	n_M = Vector{Complex{T}}(undef, sizeM)
	for i in 1:sizeM
		n_M[i] = rmls.n[i](fieldi.λ)
	end
	Δz_M = rmls.Δz
	h_M = rmls.h

	scat = get_coefficientspecifictype(rmls, fieldi)

	x_X = FFTW.fftfreq(sizeX, 1 / (fieldi.nsx_X[2] - fieldi.nsx_X[1])) * fieldi.λ
	y_Y = FFTW.fftfreq(sizeY, 1 / (fieldi.nsy_Y[2] - fieldi.nsy_Y[1])) * fieldi.λ

	#Constant to convert fft values to fourrier transform to apply the convolution theorem
	fftconst = (x_X[2] - x_X[1]) * (y_Y[2] - y_Y[1]) / 4π^2 * sizeX * sizeY * (fieldi.nsx_X[2] - fieldi.nsx_X[1]) * (fieldi.nsy_Y[2] - fieldi.nsy_Y[1]) * k^2

	in_M = n_M[sizeM:-1:1]
	ih_M = h_M[sizeM-2:-1:1]
	for iM in eachindex(Δz_M)
		for iY in eachindex(fieldi.nsy_Y)
			for iX in eachindex(fieldi.nsx_X)
				nsr = √(fieldi.nsx_X[iX]^2 + fieldi.nsy_Y[iY]^2)
				if iM == 1
					(scat.r₁₂[iX, iY], scat.t₁₂[iX,iY]) = rtss₁₂(nsr, n_M, h_M, fieldi.λ)
					(scat.r₂₁[iX, iY], scat.t₂₁[iX,iY]) = rtss₁₂(nsr, in_M, ih_M, fieldi.λ)
				end
				if nsr >= real(n_M[iM]) || nsr >= real(n_M[iM+1])
					scat.ir₁₂[iX,iY,iM] = zero(Complex{T})
					scat.sr₁₂[iX,iY,iM] = zero(Complex{T})
					scat.it₁₂[iX,iY,iM] = zero(Complex{T})
					scat.st₁₂[iX,iY,iM] = zero(Complex{T})

					scat.ir₂₁[iX,iY,iM] = zero(Complex{T})
					scat.sr₂₁[iX,iY,iM] = zero(Complex{T})
					scat.it₂₁[iX,iY,iM] = zero(Complex{T})
					scat.st₂₁[iX,iY,iM] = zero(Complex{T})

					scat.Δ[iX,iY,iM] = zero(Complex{T})
					continue
				end
				scat.ir₁₂[iX,iY,iM] = r12_i(n_M, h_M, nsr, fieldi.λ, iM)
				scat.sr₁₂[iX,iY,iM] = r12_s(n_M, h_M, nsr, fieldi.λ, iM)
				scat.it₁₂[iX,iY,iM] = t12_i(n_M, h_M, nsr, fieldi.λ, iM)
				scat.st₁₂[iX,iY,iM] = t12_s(n_M, h_M, nsr, fieldi.λ, iM)

				scat.ir₂₁[iX,iY,iM] = r12_i(in_M, ih_M, nsr, fieldi.λ, iM)
				scat.sr₂₁[iX,iY,iM] = r12_s(in_M, ih_M, nsr, fieldi.λ, iM)
				scat.it₂₁[iX,iY,iM] = t12_i(in_M, ih_M, nsr, fieldi.λ, iM)
				scat.st₂₁[iX,iY,iM] = t12_s(in_M, ih_M, nsr, fieldi.λ, iM)

				scat.Δ[iX,iY,iM] = rmls.Δz[iM](x_X[iX], y_Y[iY]) * fftconst
			end
		end
	end
	scat.ir₂₁ .*= -1
	scat.it₂₁ .*= -1

	correctscatteringmatrix_referenceframes!(scat, rmls, fieldi)
	return scat
end

function get_coefficientspecifictype(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrum{T,D,X}) where {T,D,X}
	(sizeX, sizeY, sizeM) = (length(fieldi.nsx_X), length(fieldi.nsy_Y), length(rmls.n))

	r12 = Matrix{Complex{T}}(undef,	sizeX, sizeY)
	t12 = Matrix{Complex{T}}(undef,	sizeX, sizeY)
	r21 = Matrix{Complex{T}}(undef,	sizeX, sizeY)
	t21 = Matrix{Complex{T}}(undef,	sizeX, sizeY)
	ir12 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	sr12 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	it12 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	st12 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	ir21 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	sr21 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	it21 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	st21 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	Δ = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)

	tmp2 = Array{Complex{T},3}(undef, sizeX, sizeY, sizeM-1) # preallocate matrix to avoid allocation after
	tmp1 = Array{Complex{T},3}(undef, sizeX, sizeY, sizeM-1) # preallocate matrix to avoid allocation after
	planfft = plan_fft(r12)  # precalculates the fft plan
	inv(planfft) # precalculates the inverse fft plan

	fieldl = FieldAngularSpectrum{T,-1,X}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, first(rmls.n)(fieldi.λ), dir(fieldi) > 0 ? fieldi.ref : ref1(rmls))
	fieldr = FieldAngularSpectrum{T,1,X}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, last(rmls.n)(fieldi.λ), dir(fieldi) > 0 ? ref2(rmls) : fieldi.ref)
	return RoughInterfaceConvolutionCoefficient{T, FieldAngularSpectrum{T,-1,X}, FieldAngularSpectrum{T,1,X}, Matrix{Complex{T}}, Array{Complex{T},3}, typeof(planfft)}(r12, t12, r21, t21, ir12, sr12, it12, st12, ir21, sr21, it21, st21, Δ, fieldl, fieldr, planfft, tmp1, tmp2)
end

function coefficient_general(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrum{T,D,X}) where {T,D, X<:AbstractRange}
	checkapplicability(rmls, fieldi)

	(sizeX, sizeY, sizeM) = (length(fieldi.nsx_X), length(fieldi.nsy_Y), length(rmls.n))
	n_M = Vector{Complex{T}}(undef, sizeM)
	for i in 1:sizeM
		n_M[i] = rmls.n[i](fieldi.λ)
	end

	k = 2π / fieldi.λ
	x = fftshift(FFTW.fftfreq(sizeX, 1 / (fieldi.nsx_X[2] - fieldi.nsx_X[1]))) * fieldi.λ
	y = fftshift(FFTW.fftfreq(sizeY, 1 / (fieldi.nsy_Y[2] - fieldi.nsy_Y[1]))) * fieldi.λ

	cons = (fieldi.nsx_X[2]- fieldi.nsx_X[1]) * (fieldi.nsy_Y[2] - fieldi.nsy_Y[1]) * k^2

	ir_12 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	sr_12 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	it_12 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	st_12 = Matrix{Complex{T}}(undef, sizeX, sizeY)

	ir_21 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	sr_21 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	it_21 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	st_21 = Matrix{Complex{T}}(undef, sizeX, sizeY)

	scat = get_scatteringmatrixtype(rmls, fieldi)

	in_M = n_M[sizeM:-1:1]
	ih_M = rmls.h[sizeM-2:-1:1]

	ind = LinearIndices((sizeX, sizeY))

	@inbounds for iM in 1:sizeM-1
		Threads.@threads for iX in 1:sizeX
			for iY in 1:sizeY
				nsr = √(fieldi.nsx_X[iX]^2 + fieldi.nsy_Y[iY]^2)
				ir_12[iX,iY] = r12_i(n_M, rmls.h, nsr, fieldi.λ, iM)
				sr_12[iX,iY] = r12_s(n_M, rmls.h, nsr, fieldi.λ, iM)
				it_12[iX,iY] = t12_i(n_M, rmls.h, nsr, fieldi.λ, iM)
				st_12[iX,iY] = t12_s(n_M, rmls.h, nsr, fieldi.λ, iM)

				ir_21[iX,iY] = r12_i(in_M, ih_M, nsr, fieldi.λ, iM)
				sr_21[iX,iY] = r12_s(in_M, ih_M, nsr, fieldi.λ, iM)
				it_21[iX,iY] = t12_i(in_M, ih_M, nsr, fieldi.λ, iM)
				st_21[iX,iY] = t12_s(in_M, ih_M, nsr, fieldi.λ, iM)

				if iM == 1
					i = ind[iX,iY]
					(scat.r₁₂[i,i], scat.t₁₂[i,i]) = rtss₁₂(nsr, n_M, rmls.h, fieldi.λ)
					(scat.r₂₁[i,i], scat.t₂₁[i,i]) = rtss₁₂(nsr, in_M, ih_M, fieldi.λ)
				end
			end
		end

		z = rmls.Δz[iM].(x, y')
		(nsxfft, nsyfft, fftz) = fourriertransformfft(x, y, z, fieldi.λ)
		fftz_itp = LinearInterpolation((nsxfft, nsyfft), fftz)

		Threads.@threads for iX1 in 1:sizeX
			for iY1 in 1:sizeY
				iI = ind[iX1,iY1]
				for iX2 in 1:sizeX
					for iY2 in 1:sizeY
						if (nsxfft[1] < fieldi.nsx_X[iX2] - fieldi.nsx_X[iX1] < nsxfft[sizeX] && nsyfft[1] < fieldi.nsy_Y[iY2] - fieldi.nsy_Y[iY1] < nsyfft[sizeY])
							ζ = fftz_itp(fieldi.nsx_X[iX2] - fieldi.nsx_X[iX1], fieldi.nsy_Y[iY2] - fieldi.nsy_Y[iY1])
							iS = ind[iX2,iY2]

							scat.r₁₂[iS,iI] += ζ * ir_12[iX1,iY1] * sr_12[iX2,iY2] * cons
							scat.t₁₂[iS,iI] += ζ * it_12[iX1,iY1] * st_12[iX2,iY2] * cons

							scat.r₂₁[iS,iI] -= ζ * ir_21[iX1,iY1] * sr_21[iX2,iY2] * cons
							scat.t₂₁[iS,iI] -= ζ * it_21[iX1,iY1] * st_21[iX2,iY2] * cons
						end
					end
				end
			end
		end
	end

	correctscatteringmatrix_referenceframes!(scat, rmls, fieldi)
	return scat
end

function get_scatteringmatrixtype(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrum{T,D,X}) where {T,X,D}
	(sizeX, sizeY) = (length(fieldi.nsx_X), length(fieldi.nsy_Y))
	r12 = zeros(Complex{T}, sizeX * sizeY, sizeX * sizeY)
	r21 = zeros(Complex{T}, sizeX * sizeY, sizeX * sizeY)
	t12 = zeros(Complex{T}, sizeX * sizeY, sizeX * sizeY)
	t21 = zeros(Complex{T}, sizeX * sizeY, sizeX * sizeY)

	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrum{T,-1,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[1](fieldi.λ), fieldi.ref)
		fieldr = FieldAngularSpectrum{T,1,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[end](fieldi.λ), rmls.ref)
	else
		fieldl = FieldAngularSpectrum{T,-1,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[1](fieldi.λ), rmls.ref)
		fieldr = FieldAngularSpectrum{T,1,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[end](fieldi.λ), fieldi.ref)
	end
	return ScatteringMatrix{T, FieldAngularSpectrum{T,-1,X}, FieldAngularSpectrum{T,1,X}, Matrix{Complex{T}}, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
end

@inline function ref2(rmls::RoughMultilayerStructure)
	totalh = sum(rmls.h);
	return rmls.ref + ReferenceFrame(sin(rmls.ref.θ) * cos(rmls.ref.ϕ) * totalh, sin(rmls.ref.θ) * sin(rmls.ref.ϕ) * totalh, cos(rmls.ref.θ) * totalh, rmls.ref.θ, rmls.ref.ϕ);
end

function roughfft(rmls::RoughMultilayerStructure{T}, nsx::AbstractRange, nsy::AbstractRange, λ) where T
	x = fftshift(FFTW.fftfreq(length(nsx), 1 / (nsx[2] - nsx[1]))) * λ
	y = fftshift(FFTW.fftfreq(length(nsy), 1 / (nsy[2] - nsy[1]))) * λ
	z = rmls.Δz[1].(x, y')
	(nsxfft, nsyfft, fftz) = fourriertransformfft(x, y, z, λ)
	return (x, y, z)
end

function plotroughsampling(rmls::RoughMultilayerStructure{T}, field::FieldAngularSpectrum) where T
	nsx = field.nsx_X
	nsy = field.nsy_Y
	x = fftshift(FFTW.fftfreq(length(nsx), 1 / (nsx[2] - nsx[1]))) * field.λ
	y = fftshift(FFTW.fftfreq(length(nsy), 1 / (nsy[2] - nsy[1]))) * field.λ
	z = rmls.Δz[1].(x, y')
	(nsxfft, nsyfft, fftz) = fourriertransformfft(x, y, z, field.λ)
	return (x,y,z)
end

function getfields_lr(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrum{T,A,X}) where {T,X,A}
	fieldl = FieldAngularSpectrum{T,-1,X}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, first(rmls.n_A)(fieldi.λ), ref1(rmls))
	fieldr = FieldAngularSpectrum{T,1,X}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, last(rmls.n_A)(fieldi.λ), ref2(rmls))
	return (fieldl, fieldr)
end

function checkapplicability(rmls::RoughMultilayerStructure, fieldi::FieldAngularSpectrum)
	isapprox(fieldi.n, dir(fieldi) > 0 ? rmls.n[1](fieldi.λ) : rmls.n[end](fieldi.λ), atol = @tol) || error("Field medium and rmls incident medium are different")
	checkorientation(fieldi.ref, rmls.ref) || errorToDo()
end

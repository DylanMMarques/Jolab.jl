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

n(rmls::RoughMultilayerStructure, λ) = [ni(λ) for ni in rmls.n]

function r12_i(n_M, h_M, nsr::T, λ, iM) where T
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

function r12_s(n_M, h_M, nsr::T, λ, iM) where T
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

function t12_i(n_M, h_M, nsr::T, λ, iM) where T
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

function t12_s(n_M, h_M, nsr::T, λ, iM) where T
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

function coefficient_specific(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,D, X<:AbstractRange,B}
	checkapplicability(rmls, fieldi)

	k = 2π / fieldi.λ
	n_M = n(rmls, fieldi.λ)
	scat = get_coefficientspecifictype(rmls, fieldi)

	(sizeX, sizeY) = (length(fieldi.nsx_X), length(fieldi.nsy_Y))
	x_X = FFTW.fftfreq(sizeX, 1 / (fieldi.nsx_X[2] - fieldi.nsx_X[1])) * fieldi.λ
	y_Y = FFTW.fftfreq(sizeY, 1 / (fieldi.nsy_Y[2] - fieldi.nsy_Y[1])) * fieldi.λ

	#Constant to convert fft values to fourrier transform to apply the convolution theorem
	fftconst = (x_X[2] - x_X[1]) * (y_Y[2] - y_Y[1]) / 4π^2 * sizeX * sizeY * (fieldi.nsx_X[2] - fieldi.nsx_X[1]) * (fieldi.nsy_Y[2] - fieldi.nsy_Y[1]) * k^2

	in_M = n_M[end:-1:1]
	ih_M = rmls.h[end-2:-1:1]
	cart = CartesianIndices(fieldi)
	for iM in eachindex(rmls.Δz)
		for i in eachindex(fieldi.e_SXY)
			nsr = √(fieldi.nsx_X[cart[i][2]]^2 + fieldi.nsy_Y[cart[i][3]]^2)
			if iM == 1
				(scat.r₁₂[i], scat.t₁₂[i]) = rtss₁₂(nsr, n_M, rmls.h, fieldi.λ)
				(scat.r₂₁[i], scat.t₂₁[i]) = rtss₁₂(nsr, in_M, ih_M, fieldi.λ)
			end
			scat.ir₁₂[i,iM] = r12_i(n_M, rmls.h, nsr, fieldi.λ, iM)
			scat.sr₁₂[i,iM] = r12_s(n_M, rmls.h, nsr, fieldi.λ, iM)
			scat.it₁₂[i,iM] = t12_i(n_M, rmls.h, nsr, fieldi.λ, iM)
			scat.st₁₂[i,iM] = t12_s(n_M, rmls.h, nsr, fieldi.λ, iM)

			scat.ir₂₁[i,iM] = -r12_i(in_M, ih_M, nsr, fieldi.λ, iM) # MINUS BECAUSE PERTURBATION LOOKS INVERTED
			scat.sr₂₁[i,iM] = r12_s(in_M, ih_M, nsr, fieldi.λ, iM)
			scat.it₂₁[i,iM] = -t12_i(in_M, ih_M, nsr, fieldi.λ, iM) # MINUS BECAUSE PERTURBATION LOOKS INVERTED
			scat.st₂₁[i,iM] = t12_s(in_M, ih_M, nsr, fieldi.λ, iM)

			scat.Δ[i,iM] = rmls.Δz[iM](x_X[cart[i][2]], y_Y[cart[i][2]]) * fftconst
		end
	end

	correctscatteringmatrix_referenceframes!(scat, rmls, fieldi)
	return scat
end

function get_coefficientspecifictype(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,D,X,B}
	sizeM = length(rmls.n)
	sizeI = length(fieldi.e_SXY)

	r12 = Vector{Complex{T}}(undef,	sizeI)
	t12 = Vector{Complex{T}}(undef,	sizeI)
	r21 = Vector{Complex{T}}(undef,	sizeI)
	t21 = Vector{Complex{T}}(undef,	sizeI)
	ir12 = Array{Complex{T}, 2}(undef, sizeI, sizeM-1)
	sr12 = Array{Complex{T}, 2}(undef, sizeI, sizeM-1)
	it12 = Array{Complex{T}, 2}(undef, sizeI, sizeM-1)
	st12 = Array{Complex{T}, 2}(undef, sizeI, sizeM-1)
	ir21 = Array{Complex{T}, 2}(undef, sizeI, sizeM-1)
	sr21 = Array{Complex{T}, 2}(undef, sizeI, sizeM-1)
	it21 = Array{Complex{T}, 2}(undef, sizeI, sizeM-1)
	st21 = Array{Complex{T}, 2}(undef, sizeI, sizeM-1)
	Δ = Array{Complex{T}, 2}(undef, sizeI, sizeM-1)

	tmp2 = Array{Complex{T},2}(undef, sizeI, sizeM-1) # preallocate matrix to avoid allocation after
	tmp1 = Array{Complex{T},2}(undef, sizeI, sizeM-1) # preallocate matrix to avoid allocation after
	planfft = plan_fft(reshape(r12, length(fieldi.nsx_X), length(fieldi.nsy_Y)))  # precalculates the fft plan
	inv(planfft) # precalculates the inverse fft plan

	fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, first(rmls.n)(fieldi.λ), dir(fieldi) > 0 ? fieldi.ref : ref1(rmls))
	fieldr = FieldAngularSpectrumScalar{T,1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, last(rmls.n)(fieldi.λ), dir(fieldi) > 0 ? ref2(rmls) : fieldi.ref)
	return RoughInterfaceConvolutionCoefficient{T, FieldAngularSpectrumScalar{T,-1,X,B}, FieldAngularSpectrumScalar{T,1,X,B}, Vector{Complex{T}}, Matrix{Complex{T}}, typeof(planfft)}(r12, t12, r21, t21, ir12, sr12, it12, st12, ir21, sr21, it21, st21, Δ, fieldl, fieldr, planfft, tmp1, tmp2)
end

function coefficient_general(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,D, X<:AbstractRange,B}
	checkapplicability(rmls, fieldi)

	sizeI = length(fieldi.e_SXY)
	n_M = n(rmls, fieldi.λ)
	in_M = n_M[end:-1:1]
	ih_M = rmls.h[end-2:-1:1]

	k = 2π / fieldi.λ
	cons = (fieldi.nsx_X[2]- fieldi.nsx_X[1]) * (fieldi.nsy_Y[2] - fieldi.nsy_Y[1]) * k^2

	ir_12 = Vector{Complex{T}}(undef, sizeI)
	sr_12 = Vector{Complex{T}}(undef, sizeI)
	it_12 = Vector{Complex{T}}(undef, sizeI)
	st_12 = Vector{Complex{T}}(undef, sizeI)

	ir_21 = Vector{Complex{T}}(undef, sizeI)
	sr_21 = Vector{Complex{T}}(undef, sizeI)
	it_21 = Vector{Complex{T}}(undef, sizeI)
	st_21 = Vector{Complex{T}}(undef, sizeI)

	scat = get_scatteringmatrixtype(rmls, fieldi)

	cart_i = CartesianIndices(fieldi)
	#MISSING stuff
	for iM in eachindex(rmls.Δz)
		for i in 1:sizeI
			nsr = √(fieldi.nsx_X[cart_i[i][2]]^2 + fieldi.nsy_Y[cart_i[i][3]]^2)
			ir_12[i] = r12_i(n_M, rmls.h, nsr, fieldi.λ, iM)
			sr_12[i] = r12_s(n_M, rmls.h, nsr, fieldi.λ, iM)
			it_12[i] = t12_i(n_M, rmls.h, nsr, fieldi.λ, iM)
			st_12[i] = t12_s(n_M, rmls.h, nsr, fieldi.λ, iM)
			ir_21[i] = r12_i(in_M, ih_M, nsr, fieldi.λ, iM)
			sr_21[i] = r12_s(in_M, ih_M, nsr, fieldi.λ, iM)
			it_21[i] = t12_i(in_M, ih_M, nsr, fieldi.λ, iM)
			st_21[i] = t12_s(in_M, ih_M, nsr, fieldi.λ, iM)
			if iM == 1
				(scat.r₁₂[i,i], scat.t₁₂[i,i]) = rtss₁₂(nsr, n_M, rmls.h, fieldi.λ)
				(scat.r₂₁[i,i], scat.t₂₁[i,i]) = rtss₁₂(nsr, in_M, ih_M, fieldi.λ)
			end
		end

		(sizeX, sizeY) = size(cart_i)[2:3]
		x = fftshift(FFTW.fftfreq(sizeX, 1 / (fieldi.nsx_X[2] - fieldi.nsx_X[1]))) * fieldi.λ
		y = fftshift(FFTW.fftfreq(sizeY, 1 / (fieldi.nsy_Y[2] - fieldi.nsy_Y[1]))) * fieldi.λ
		z = rmls.Δz[iM].(x, y')
		(nsxfft, nsyfft, fftz) = fourriertransformfft(x, y, z, fieldi.λ)
		fftz_itp = LinearInterpolation((nsxfft, nsyfft), fftz)

		Threads.@threads for i_i in iterator_index(fieldi)
			iX1, iY1 = cart_i[i_i][2], cart_i[i_i][3]
			for i_s in iterator_index(fieldi)
				iX2, iY2 = cart_i[i_s][2], cart_i[i_s][3]
				if (nsxfft[1] < fieldi.nsx_X[iX2] - fieldi.nsx_X[iX1] < nsxfft[sizeX] && nsyfft[1] < fieldi.nsy_Y[iY2] - fieldi.nsy_Y[iY1] < nsyfft[sizeY])
					ζ = fftz_itp(fieldi.nsx_X[iX2] - fieldi.nsx_X[iX1], fieldi.nsy_Y[iY2] - fieldi.nsy_Y[iY1])
					scat.r₁₂[i_s,i_i] += ζ * ir_12[i_i] * sr_12[i_s] * cons
					scat.t₁₂[i_s,i_i] += ζ * it_12[i_i] * st_12[i_s] * cons
					scat.r₂₁[i_s,i_i] -= ζ * ir_21[i_i] * sr_21[i_s] * cons
					scat.t₂₁[i_s,i_i] -= ζ * it_21[i_i] * st_21[i_s] * cons
				end
			end
		end
	end

	correctscatteringmatrix_referenceframes!(scat, rmls, fieldi)
	return scat
end

function get_scatteringmatrixtype(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,X,D,B}
	sizeI = length(fieldi.e_SXY)
	r12 = zeros(Complex{T}, sizeI, sizeI)
	r21 = zeros(Complex{T}, sizeI, sizeI)
	t12 = zeros(Complex{T}, sizeI, sizeI)
	t21 = zeros(Complex{T}, sizeI, sizeI)

	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[1](fieldi.λ), fieldi.ref)
		fieldr = FieldAngularSpectrumScalar{T,1,X,B}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[end](fieldi.λ), rmls.ref)
	else
		fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[1](fieldi.λ), rmls.ref)
		fieldr = FieldAngularSpectrumScalar{T,1,X,B}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[end](fieldi.λ), fieldi.ref)
	end
	return ScatteringMatrix{T, FieldAngularSpectrumScalar{T,-1,X,B}, FieldAngularSpectrumScalar{T,1,X,B}, Matrix{Complex{T}}, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
end

function ref2(rmls::RoughMultilayerStructure)
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

function plotroughsampling(rmls::RoughMultilayerStructure{T}, field::FieldAngularSpectrumScalar) where T
	nsx = field.nsx_X
	nsy = field.nsy_Y
	x = fftshift(FFTW.fftfreq(length(nsx), 1 / (nsx[2] - nsx[1]))) * field.λ
	y = fftshift(FFTW.fftfreq(length(nsy), 1 / (nsy[2] - nsy[1]))) * field.λ
	z = rmls.Δz[1].(x, y')
	(nsxfft, nsyfft, fftz) = fourriertransformfft(x, y, z, field.λ)
	return (x,y,z)
end

function getfields_lr(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrumScalar{T,A,X,B}) where {T,X,A,B}
	fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, first(rmls.n_A)(fieldi.λ), ref1(rmls))
	fieldr = FieldAngularSpectrumScalar{T,1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, last(rmls.n_A)(fieldi.λ), ref2(rmls))
	return (fieldl, fieldr)
end

function checkapplicability(rmls::RoughMultilayerStructure, fieldi::AbstractFieldAngularSpectrum)
	isapprox(fieldi.n, dir(fieldi) > 0 ? rmls.n[1](fieldi.λ) : rmls.n[end](fieldi.λ), atol = @tol) || error("Field medium and rmls incident medium are different")
	checkorientation(fieldi.ref, rmls.ref) || errorToDo()
end

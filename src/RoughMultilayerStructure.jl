struct RoughMultilayerStructure{T<:Real} <: AbstractOpticalComponent{T}
	n::Vector{JolabFunction1D{T,Complex{T}}}
	h::Vector{T}
	Δz::Vector{JolabFunction2D{T,T}}
	ref::ReferenceFrame{T}
	function RoughMultilayerStructure{T}(n, h, Δz, ref) where T
		length(n) == length(h) + 2 == length(Δz) + 1 || error("wrong sizes")
		return new{T}(n, h, Δz, ref)
	end
end
RoughMultilayerStructure(n1, n2, Δz, ref) = RoughMultilayerStructure{Float64}(n1, n2, Δz, ref)

function coefficient_specific(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrum{T,X}) where {T, X<:AbstractRange}
	checkorientation(rmls.ref, fieldi.ref) || tobedone()
	(sizeX, sizeY, sizeM) = (length(fieldi.nsx_X), length(fieldi.nsy_Y), length(rmls.n))
	isapprox(fieldi.n, fieldi.dir > 0 ? rmls.n[1](fieldi.λ) : rmls.n[sizeM](fieldi.λ), atol = @tol) || error("Refractive index missmatch")

	k = 2π / fieldi.λ
	n_M = Vector{Complex{T}}(undef, sizeM)
	for i in 1:sizeM
		n_M[i] = rmls.n[i](fieldi.λ)
	end
	Δz_M = rmls.Δz
	h_M = rmls.h

	r12 = Matrix{Complex{T}}(undef,sizeX, sizeY)
	t12 = Matrix{Complex{T}}(undef,sizeX, sizeY)
	r21 = Matrix{Complex{T}}(undef,sizeX, sizeY)
	t21 = Matrix{Complex{T}}(undef,sizeX, sizeY)
	ir12 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	sr12 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	it12 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	st12 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	ir21 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	sr21 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	it21 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	st21 = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)
	Δ = Array{Complex{T}, 3}(undef, sizeX, sizeY, sizeM-1)

	x_X = FFTW.fftfreq(sizeX, 1 / (fieldi.nsx_X[2] - fieldi.nsx_X[1])) * fieldi.λ
	y_Y = FFTW.fftfreq(sizeY, 1 / (fieldi.nsy_Y[2] - fieldi.nsy_Y[1])) * fieldi.λ

	#Constant to convert fft values to fourrier transform to apply the convolution theorem
	fftconst = (x_X[2] - x_X[1]) * (y_Y[2] - y_Y[1]) / 4π^2 * sizeX * sizeY * (fieldi.nsx_X[2] - fieldi.nsx_X[1]) * (fieldi.nsy_Y[2] - fieldi.nsy_Y[1]) * k^2

	for aux in 1:2
		if aux == 2
			n_M = n_M[sizeM:-1:1]
			Δz_M = Δz_M[sizeM-1:-1:1]
			h_M = rmls.h[sizeM-2:-1:1]

			it = it21
			st = st21
			ir = ir21
			sr = sr21
			r = r21
			t = t21
		else
			it = it12
			st = st12
			ir = ir12
			sr = sr12
			r = r12
			t = t12
		end
		Threads.@threads for iM in eachindex(Δz_M)
			for iY in eachindex(fieldi.nsy_Y)
				for iX in eachindex(fieldi.nsx_X)
					nsr = √(fieldi.nsx_X[iX]^2 + fieldi.nsy_Y[iY]^2)
					sz_m = √(1 - nsr^2 / n_M[iM]^2)
					sz_m1 = √(1 - nsr^2 / n_M[iM+1]^2)
					(rm_m1, tm_m1) =  rtss₁₂(nsr, view(n_M, iM:sizeM), view(h_M, iM:sizeM-2), fieldi.λ)
					if iM > 1
						(tmp, t0_m) =  rtss₁₂(nsr, view(n_M, 1:iM), view(h_M, 1:iM-2), fieldi.λ)
						(rm_ml1, tm_0) =  rtss₁₂(nsr, view(n_M, iM:-1:1), view(h_M, iM-2:-1:1), fieldi.λ)
						phaseTerm_m = exp(im * k * sz_m * h_M[iM-1])
						t0_m *= phaseTerm_m
						tm_0 *= phaseTerm_m
						M_m = 1 - rm_ml1 * rm_m1 * phaseTerm_m^2
					else
						t0_m = one(Complex{T})
						tm_0 = one(Complex{T})
						M_m = one(Complex{T})
					end
					ir[iX,iY,iM] = -im * k / 2 * (n_M[iM+1]^2 - n_M[iM]^2) * t0_m * (1 + rm_m1)
					sr[iX,iY,iM] = tm_0 * (1 + rm_m1) / n_M[iM] / sz_m / M_m

					(rm1_0, tm1_0) = rtss₁₂(nsr, view(n_M, iM+1:-1:1), view(h_M, iM-1:-1:1), fieldi.λ)
					if iM < sizeM - 1
						(tmp, tm1_M) = rtss₁₂(nsr, view(n_M, iM+1:sizeM), view(h_M, iM+1:sizeM-2), fieldi.λ)
						tm1_M *= exp(im * k * sz_m1 * h_M[iM])
					else
						tm1_M = one(Complex{T})
					end
					it[iX,iY,iM] = -im * k / 2 * (n_M[iM+1]^2 - n_M[iM]^2) * t0_m * (1 + rm_m1)
					st[iX,iY,iM] = tm1_M * (1 + rm1_0)/ n_M[iM+1] / sz_m1

					aux == 1 && (Δ[iX,iY,iM] = rmls.Δz[iM](x_X[iX], y_Y[iY]) * fftconst)

					if iM == 1
						(r[iX, iY], t[iX,iY]) = rtss₁₂(nsr, n_M, h_M, fieldi.λ)
					end
				end
			end
		end
	end
	ir21 .*= -1
	it21 .*= -1

	ref = fieldi.dir > 0 ? ref1(rmls) : ref2(rmls)
	if !checkposition(fieldi.ref, ref)
		propM = propagationmatrix(fieldi, ref)
		ind = LinearIndices((sizeX, sizeY))
		if fieldi.dir > 0
			vec(r12) .*= propM.diag .* propM.diag
			vec(t12) .*= propM.diag
			vec(t21) .*= propM.diag
			@inbounds @simd for iM in 1:sizeM-1
				for iY in 1:sizeY
					for iX in 1:sizeX
						i = ind[iX,iY]
						ir12[iX,iY,iM] *= propM.diag[i]
						sr12[iX,iY,iM] *= propM.diag[i]
						it12[iX,iY,iM] *= propM.diag[i]
						st21[iX,iY,iM] *= propM.diag[i]
					end
				end
			end
		else
			conj!(propM.diag)
			vec(r21) .*= propM.diag .* propM.diag
			vec(t12) .*= propM.diag
			vec(t21) .*= propM.diag
			@inbounds @simd for iM in 1:sizeM-1
				for iY in 1:sizeY
					for iX in 1:sizeX
						i = ind[iX,iY]
						ir21[iX,iY,iM] *= propM.diag[i]
						sr21[iX,iY,iM] *= propM.diag[i]
						st12[iX,iY,iM] *= propM.diag[i]
						it21[iX,iY,iM] *= propM.diag[i]
					end
				end
			end
		end
	end
	if fieldi.dir > 0
		fieldl = FieldAngularSpectrum{T,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, -1, fieldi.ref)
		fieldr = FieldAngularSpectrum{T,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[sizeM](fieldi.λ), 1, ref2(rmls))
	else
		fieldl = FieldAngularSpectrum{T,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[1](fieldi.λ), -1, ref1(rmls))
		fieldr = FieldAngularSpectrum{T,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, 1, fieldi.ref)
	end

	tmp2 = Array{Complex{T},3}(undef, sizeX, sizeY, sizeM-1) # preallocate matrix to avoid allocation after
	tmp1 = Array{Complex{T},3}(undef, sizeX, sizeY, sizeM-1) # preallocate matrix to avoid allocation after
	planfft = plan_fft(r12)  # precalculates the fft plan
	inv(planfft) # precalculates the inverse fft plan
	return RoughInterfaceConvolutionCoefficient{T, FieldAngularSpectrum{T,X}, FieldAngularSpectrum{T,X}, Matrix{Complex{T}}, Array{Complex{T},3}, typeof(planfft)}(r12, t12, r21, t21, ir12, sr12, it12, st12, ir21, sr21, it21, st21, Δ, fieldl, fieldr, planfft, tmp1, tmp2)
end

function coefficient_general(rmls::RoughMultilayerStructure{T}, fieldi::FieldAngularSpectrum{T,X}) where {T, X<:AbstractRange}
	isapprox(fieldi.n, fieldi.dir > 0 ? rmls.n[1](fieldi.λ) : rmls.n[end](fieldi.λ), atol = @tol) || error("Field medium and rmls incident medium are different")
	checkorientation(fieldi.ref, rmls.ref) || errorToDo()

	errorToDo()

	sizeX = length(fieldi.nsx_X)
	sizeY = length(fieldi.nsy_Y)
	k = 2π / fieldi.λ
	n1 = rmls.n[1](fieldi.λ)
	n2 = rmls.n2(fieldi.λ)
	x = fftshift(FFTW.fftfreq(sizeX, 1 / (fieldi.nsx_X[2] - fieldi.nsx_X[1]))) * fieldi.λ
	y = fftshift(FFTW.fftfreq(sizeY, 1 / (fieldi.nsy_Y[2] - fieldi.nsy_Y[1]))) * fieldi.λ
	z = rmls.Δz.(x, y')
	(nsxfft, nsyfft, fftz) = fourriertransformfft(x, y, z, fieldi.λ)

	fftz_itp = LinearInterpolation((nsxfft, nsyfft), fftz)

	r12 = zeros(Complex{T}, sizeX, sizeY, sizeX, sizeY)
	t12 = zeros(Complex{T}, sizeX, sizeY, sizeX, sizeY)
	r21 = zeros(Complex{T}, sizeX, sizeY, sizeX, sizeY)
	t21 = zeros(Complex{T}, sizeX, sizeY, sizeX, sizeY)
	cons = (fieldi.nsx_X[2]- fieldi.nsx_X[1]) * (fieldi.nsy_Y[2] - fieldi.nsy_Y[1]) * k^2
	@inbounds Threads.@threads for iX1 in 1:sizeX
		@simd for iY1 in 1:sizeY
			sz1_i = √(1 - (fieldi.nsx_X[iX1]^2 + fieldi.nsy_Y[iY1]^2) / n1^2)
			sz2_i = √(1 - (fieldi.nsx_X[iX1]^2 + fieldi.nsy_Y[iY1]^2) / n2^2)
			ri12 = reflectioncoefficientinterfaces(n1, sz1_i, n2, sz2_i)
			for iX2 in 1:sizeX
				for iY2 in 1:sizeY
					for iM in 1:sizeM
						if (nsxfft[1] < fieldi.nsx_X[iX2] - fieldi.nsx_X[iX1] < nsxfft[sizeX] && nsyfft[1] < fieldi.nsy_Y[iY2] - fieldi.nsy_Y[iY1] < nsyfft[sizeY])
							aux = fftz_itp(fieldi.nsx_X[iX2] - fieldi.nsx_X[iX1], fieldi.nsy_Y[iY2] - fieldi.nsy_Y[iY1])
							sz1_s = √(1 - (fieldi.nsx_X[iX2]^2 + fieldi.nsy_Y[iY2]^2) / n1^2)
							sz2_s = √(1 - (fieldi.nsx_X[iX2]^2 + fieldi.nsy_Y[iY2]^2) / n2^2)
							rs12 = reflectioncoefficientinterfaces(n1, sz1_s, n2, sz2_s)
							ts21 = transmissioncoefficientinterfaces(n2, sz2_s, n1, sz1_s)
							ts12 = transmissioncoefficientinterfaces(n1, sz1_s, n2, sz2_s)
							# add the minus sign because the z is upside down on the paper
							r12[iX2, iY2, iX1, iY1] = -im * aux * (n2^2 - n1^2) / 2 * k / n1 / sz1_s * (1 + ri12) * (1 + rs12) * cons
							t12[iX2, iY2, iX1, iY1] = -im * k / 2 / n2 / sz2_s * (n2^2 - n1^2) * aux * (1 + ri12) * ts21 * cons
							r21[iX2, iY2, iX1, iY1] = im * aux * (n1^2 - n2^2) / 2 * k / n2 / sz2_s * (1 - ri12) * (1 - rs12) * cons
							t21[iX2, iY2, iX1, iY1] = im * k / 2 / n1 / sz1_s * (n1^2 - n2^2) * aux * (1 - ri12) * ts12 * cons
						end
					end
					if (iX1 == iX2 && iY1 == iY2)
						r12[iX2, iY2, iX1, iY1] += ri12
						t12[iX2, iY2, iX1, iY1] += ts12
						r21[iX2, iY2, iX1, iY1] -= ri12
						t21[iX2, iY2, iX1, iY1] += ts21
					end
				end
			end
		end
	end
	r12 = reshape(r12, sizeX * sizeY, sizeX * sizeY)
	t12 = reshape(t12, sizeX * sizeY, sizeX * sizeY)
	r21 = reshape(r21, sizeX * sizeY, sizeX * sizeY)
	t21 = reshape(t21, sizeX * sizeY, sizeX * sizeY)

	if !checkposition(fieldi.ref, rmls.ref)
		propM = propagationmatrix(fieldi, rmls.ref)
		if fieldi.dir > 0
			rmul!(r12, propM)
			lmul!(propM, r12)
			rmul!(t12, propM)
			lmul!(propM, t21)
		else
			conj!(propM.diag)
			rmul!(r21, propM)
			lmul!(propM, r21)
			rmul!(t21, propM)
			lmul!(propM, t12)
		end
	end

	if fieldi.dir > 0
		fieldl = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[1](fieldi.λ), -1, fieldi.ref)
		fieldr = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[end](fieldi.λ), 1, rmls.ref)
	else
		fieldl = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[1](fieldi.λ), -1, rmls.ref)
		fieldr = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n[end](fieldi.λ), 1, fieldi.ref)
	end
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), Matrix{Complex{T}}, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
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

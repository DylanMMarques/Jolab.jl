struct RoughInterface{T<:Real} <: AbstractOpticalComponent{T}
	n1::JolabFunction1D{T,Complex{T}}
	n2::JolabFunction1D{T,Complex{T}}
	Δz::JolabFunction2D{T,T}
	ref::ReferenceFrame{T}
	RoughInterface{T}(n1, n2, Δz, ref) where T = new{T}(n1, n2, Δz, ref)
end
RoughInterface(n1, n2, Δz, ref) = RoughInterface{Float64}(n1, n2, Δz, ref)

function coefficient_specific(rmls::RoughInterface{T}, fieldi::FieldAngularSpectrum{T,X}) where {T, X<:AbstractRange}
	checkorientation(rmls.ref, fieldi.ref) || tobedone()
	isapprox(fieldi.n, fieldi.dir > 0 ? rmls.n1(fieldi.λ) : rmls.n2(fieldi.λ), atol = @tol) || error("Refractive index missmatch")

	k = 2π / fieldi.λ
	n1 = rmls.n1(fieldi.λ)
	n2 = rmls.n2(fieldi.λ)
	(sizeX, sizeY) = (length(fieldi.nsx_X), length(fieldi.nsy_Y))

	r12 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	t12 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	r21 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	ir12 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	sr12 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	it12 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	st12 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	ir21 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	sr21 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	it21 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	st21 = Matrix{Complex{T}}(undef, sizeX, sizeY)
	Δ = Matrix{Complex{T}}(undef, sizeX, sizeY)

	x_X = FFTW.fftfreq(sizeX, 1 / (fieldi.nsx_X[2] - fieldi.nsx_X[1])) * fieldi.λ
	y_Y = FFTW.fftfreq(sizeY, 1 / (fieldi.nsy_Y[2] - fieldi.nsy_Y[1])) * fieldi.λ

	#Constant to convert fft values to fourrier transform to apply the convolution theorem
	fftconst = (x_X[2] - x_X[1]) * (y_Y[2] - y_Y[1]) / 4π^2 * sizeX * sizeY * (fieldi.nsx_X[2] - fieldi.nsx_X[1]) * (fieldi.nsy_Y[2] - fieldi.nsy_Y[1]) * k^2
	@inbounds Threads.@threads for iY in eachindex(fieldi.nsy_Y)
		@simd for iX in eachindex(fieldi.nsx_X)
			nsr² = fieldi.nsx_X[iX]^2 + fieldi.nsy_Y[iY]^2
			sz1 = √(1 - nsr² / n1^2)
			sz2 = √(1 - nsr² / n2^2)
			r12[iX,iY] = reflectioncoefficientinterfaces(n1, sz1, n2, sz2)
			r21[iX,iY] = -r12[iX,iY]
			t21[iX,iY] = transmissioncoefficientinterfaces(n2, sz2, n1, sz1)
			t12[iX,iY] = transmissioncoefficientinterfaces(n1, sz1, n2, sz2)

			ir12[iX,iY] = -im * k / 2 * (n2^2 - n1^2)  * (1 + r12[iX,iY])
			sr12[iX,iY] = (1 + r12[iX,iY]) / n1 / sz1
			Δ[iX,iY] = rmls.Δz(x_X[iX], y_Y[iY]) * fftconst

			it12[iX,iY] = -im * k / 2 * (n2^2 - n1^2)  * (1 + r12[iX,iY])
			st12[iX,iY] = t21[iX,iY] / n2 / sz2

			ir21[iX,iY] = im * k / 2 * (n1^2 - n2^2)  * (1 + r21[iX,iY])
			sr21[iX,iY] = (1 + r21[iX,iY]) / n2 / sz2

			it21[iX,iY] = im * k / 2 * (n1^2 - n2^2)  * (1 + r21[iX,iY])
			st21[iX,iY] = t12[iX,iY] / n1 / sz1
		end
	end

	if !checkposition(fieldi.ref, rmls.ref)
		propM = propagationmatrix(fieldi, rmls.ref)
		if fieldi.dir > 0
			vec(r12) .*= propM.diag .* propM.diag
			vec(ir12) .*= propM.diag
			vec(sr12) .*= propM.diag
			vec(t12) .*= propM.diag
			vec(it12) .*= propM.diag
			vec(t21) .*= propM.diag
			vec(st21) .*= propM.diag
		else
			conj!(propM.diag)
			vec(r21) .*= propM.diag .* propM.diag
			vec(ir21) .*= propM.diag
			vec(sr21) .*= propM.diag
			vec(t12) .*= propM.diag
			vec(st12) .*= propM.diag
			vec(t21) .*= propM.diag
			vec(it21) .*= propM.diag
		end
	end
	if fieldi.dir > 0
		fieldl = FieldAngularSpectrum{T,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, -1, fieldi.ref)
		fieldr = FieldAngularSpectrum{T,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n2(fieldi.λ), 1, rmls.ref)
	else
		fieldl = FieldAngularSpectrum{T,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n1(fieldi.λ), -1, rmls.ref)
		fieldr = FieldAngularSpectrum{T,X}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, 1, fieldi.ref)
	end
	tmp = Matrix{Complex{T}}(undef, sizeX, sizeY) # preallocate matrix to avoid allocation after
	planfft = plan_fft(tmp)  # precalculates the fft plan
	inv(planfft) # precalculates the inverse fft plan

	return RoughInterfaceConvolutionCoefficient{T, FieldAngularSpectrum{T,X}, FieldAngularSpectrum{T,X}, Matrix{Complex{T}}, typeof(planfft)}(r12, t12, r21, t21, ir12, sr12, it12, st12, ir21, sr21, it21, st21, Δ, fieldl, fieldr, planfft, tmp)
end

function roughfft(rmls::RoughInterface{T}, nsx::AbstractRange, nsy::AbstractRange, λ) where T
	x = fftshift(FFTW.fftfreq(length(nsx), 1 / (nsx[2] - nsx[1]))) * λ
	y = fftshift(FFTW.fftfreq(length(nsy), 1 / (nsy[2] - nsy[1]))) * λ
	z = rmls.Δz.(x, y')
	(nsxfft, nsyfft, fftz) = fourriertransformfft(x, y, z, λ)
	return (x,y)
end

function coefficient_general(rmls::RoughInterface{T}, fieldi::FieldAngularSpectrum{T,X}) where {T, X<:AbstractRange}
	isapprox(fieldi.n, fieldi.dir > 0 ? rmls.n1(fieldi.λ) : rmls.n2(fieldi.λ), atol = @tol) || error("Field medium and rmls incident medium are different")
	checkorientation(fieldi.ref, rmls.ref) || errorToDo()

	sizeX = length(fieldi.nsx_X)
	sizeY = length(fieldi.nsy_Y)
	k = 2π / fieldi.λ
	n1 = rmls.n1(fieldi.λ)
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
		fieldl = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n1(fieldi.λ), -1, fieldi.ref)
		fieldr = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n2(fieldi.λ), 1, rmls.ref)
	else
		fieldl = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n1(fieldi.λ), -1, rmls.ref)
		fieldr = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, rmls.n1(fieldi.λ), 1, fieldi.ref)
	end
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), Matrix{Complex{T}}, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
end

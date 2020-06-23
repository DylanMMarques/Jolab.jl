mutable struct RoughInterface{T<:Real} <: AbstractOpticalComponent{T}
	n1::JolabFunction1D{T,Complex{T}}
	n2::JolabFunction1D{T,Complex{T}}
	Δz::JolabFunction2D{T,T}
	ref::ReferenceFrame{T}
	RoughInterface{T}(n1, n2, Δz, ref) where T = new{T}(n1, n2, Δz, ref)
end
RoughInterface(n1, n2, Δz, ref) = RoughInterface{Float64}(n1, n2, Δz, ref)

function ScatteringConvolutionCoefficientScalar_maxtrixform(rmls::RoughInterface{T}, nsx_X::AbstractRange{<:Real}, nsy_Y::AbstractRange{<:Real}, λ::Real) where T
	k = convert(T, 2π / λ)
	n1 = rmls.n1(λ)
	n2 = rmls.n2(λ)
	(sizeX, sizeY) = (length(nsx_X), length(nsy_Y))

	ir12 = Array{Complex{T}, 2}(undef, sizeX, sizeY)
	sr12 = Array{Complex{T}, 2}(undef, sizeX, sizeY)
	Δr12 = Array{Complex{T}, 2}(undef, sizeX, sizeY)
	it12 = Array{Complex{T}, 2}(undef, sizeX, sizeY)
	st12 = Array{Complex{T}, 2}(undef, sizeX, sizeY)
	ir21 = Array{Complex{T}, 2}(undef, sizeX, sizeY)
	sr21 = Array{Complex{T}, 2}(undef, sizeX, sizeY)
	it21 = Array{Complex{T}, 2}(undef, sizeX, sizeY)
	st21 = Array{Complex{T}, 2}(undef, sizeX, sizeY)

	x_X = FFTW.fftfreq(sizeX, 1 / (nsx_X[2] - nsx_X[1])) * λ
	y_Y = FFTW.fftfreq(sizeY, 1 / (nsy_Y[2] - nsy_Y[1])) * λ

	#Constant to convert fft values to fourrier transform to apply the convolution theorem
	fftconst = (x_X[2] - x_X[1]) * (y_Y[2] - y_Y[1]) / 4π^2 * sizeX * sizeY * (nsx_X[2] - nsx_X[1]) * (nsy_Y[2] - nsy_Y[1]) * k^2
	@inbounds Threads.@threads for iY in eachindex(nsy_Y)
		@simd for iX in eachindex(nsx_X)
			nsr² = nsx_X[iX]^2 + nsy_Y[iY]^2
			sz1 = √(1 - nsr² / n1^2)
			sz2 = √(1 - nsr² / n2^2)
			r12 = reflectioncoefficientinterfaces(n1, sz1, n2, sz2)
			t21 = transmissioncoefficientinterfaces(n2, sz2, n1, sz1)
			t12 = transmissioncoefficientinterfaces(n1, sz1, n2, sz2)

			ir12[iX,iY] = -im * k / 2 * (n2^2 - n1^2)  * (1 + r12)
			sr12[iX,iY] = (1 + r12) / n1 / sz1
			Δr12[iX,iY] = rmls.Δz(x_X[iX], y_Y[iY]) * fftconst

			it12[iX,iY] = -im * k / 2 * (n2^2 - n1^2)  * (1 + r12)
			st12[iX,iY] = t21 / n2 / sz2

			ir21[iX,iY] = im * k / 2 * (n1^2 - n2^2)  * (1 - r12)
			sr21[iX,iY] = (1 - r12) / n2 / sz2

			it21[iX,iY] = im * k / 2 * (n1^2 - n2^2)  * (1 - r12)
			st21[iX,iY] = t12 / n1 / sz1
		end
	end
	return ScatteringConvolutionCoefficientScalar{T, Array{Complex{T},2}, Array{Complex{T},2}}(ir12, sr12, Δr12, it12, st12, Δr12, ir21, sr21, Δr12, ir21, st21, Δr12, λ, n1, rmls.ref, n2, rmls.ref);
end

function roughfft(rmls::RoughInterface{T}, nsx::AbstractRange, nsy::AbstractRange, λ) where T
	x = fftshift(FFTW.fftfreq(length(nsx), 1 / (nsx[2] - nsx[1]))) * λ
	y = fftshift(FFTW.fftfreq(length(nsy), 1 / (nsy[2] - nsy[1]))) * λ
	z = rmls.Δz.(x, y')
	(nsxfft, nsyfft, fftz) = fourriertransformfft(x, y, z, λ)
	return fftz
end

function coefficientscallar(rmls::RoughInterface{T}, nsx::AbstractRange{<:Real}, nsy::AbstractRange{<:Real}, λ::Real) where T
	sizeX = length(nsx)
	sizeY = length(nsy)
	k = convert(T, 2π / λ)
	n1 = rmls.n1(λ)
	n2 = rmls.n2(λ)
	x = fftshift(FFTW.fftfreq(sizeX, 1 / (nsx[2] - nsx[1]))) * λ
	y = fftshift(FFTW.fftfreq(sizeY, 1 / (nsy[2] - nsy[1]))) * λ
	z = rmls.Δz.(x, y')
	(nsxfft, nsyfft, fftz) = fourriertransformfft(x, y, z, λ)

	fftz_itp = LinearInterpolation((nsxfft, nsyfft), fftz)

	r12 = zeros(Complex{T}, sizeX, sizeY, sizeX, sizeY)
	t12 = zeros(Complex{T}, sizeX, sizeY, sizeX, sizeY)
	r21 = zeros(Complex{T}, sizeX, sizeY, sizeX, sizeY)
	t21 = zeros(Complex{T}, sizeX, sizeY, sizeX, sizeY)
	cons = (nsx[2]- nsx[1]) * (nsy[2] - nsy[1]) * k^2
	@inbounds Threads.@threads for iX1 in 1:sizeX
		@simd for iY1 in 1:sizeY
			sz1_i = √(1 - (nsx[iX1]^2 + nsy[iY1]^2) / n1^2)
			sz2_i = √(1 - (nsx[iX1]^2 + nsy[iY1]^2) / n2^2)
			ri12 = reflectioncoefficientinterfaces(n1, sz1_i, n2, sz2_i)
			for iX2 in 1:sizeX
				for iY2 in 1:sizeY
					if (nsxfft[1] < nsx[iX2] - nsx[iX1] < nsxfft[sizeX] && nsyfft[1] < nsy[iY2] - nsy[iY1] < nsyfft[sizeY])
						aux = fftz_itp(nsx[iX2] - nsx[iX1], nsy[iY2] - nsy[iY1])
						sz1_s = √(1 - (nsx[iX2]^2 + nsy[iY2]^2) / n1^2)
						sz2_s = √(1 - (nsx[iX2]^2 + nsy[iY2]^2) / n2^2)
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
	return PropagationCoefficientScalar{T,Array{Complex{T},2}}(r12, t12, r21, t21, λ, n1, rmls.ref, n2, rmls.ref)
end

function PropagationScatteringConvolutionCoefficientScalar(rmls::RoughInterface{T}, λ::Real) where T
	propCoef = PropagationCoefficientScalar([rmls.n1(λ), rmls.n2(λ)], zeros(T,0), λ, rmls.ref)
	scatConvCoef = ScatteringConvolutionCoefficientScalar(rmls, λ)
	return PropagationScatteringConvolutionCoefficientScalar{T}(propCoef, scatConvCoef);
end

function PropagationScatteringConvolutionCoefficientScalar_matrixform(rmls::RoughInterface{T}, nsx_X::AbstractRange{<:Real}, nsy_Y::AbstractRange{<:Real}, λ::Real) where T
	propCoef = PropagationCoefficientScalar_matrixform([rmls.n1(λ), rmls.n2(λ)], zeros(T,0), λ, rmls.ref)
	scatConvCoef = ScatteringConvolutionCoefficientScalar_matrixform(rmls, λ)
	return PropagationScatteringConvolutionCoefficientScalar{T}(propCoef, scatConvCoef);
end

coefficientscallar(rmls::RoughInterface, λ) = PropagationScatteringConvolutionCoefficientScalar(rmls, λ)

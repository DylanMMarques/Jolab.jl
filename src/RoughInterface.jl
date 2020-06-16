mutable struct RoughInterface{T<:Real} <: AbstractOpticalComponent{T}
	n1::JolabFunction1D{T,Complex{T}}
	n2::JolabFunction1D{T,Complex{T}}
	Δz::JolabFunction2D{T,T}
	ref::ReferenceFrame{T}
	RoughInterface{T}(n1, n2, Δz, ref) where T = new{T}(n1, n2, Δz, ref)
end
RoughInterface(n1, n2, Δz, ref) = RoughInterface{Float64}(n1, n2, Δz, ref)

function ScatteringConvolutionCoefficientScalar(rmls::RoughInterface{T}, λ::Real) where T
	k = convert(T, 2π / λ)
	n1 = rmls.n1(λ)
	n2 = rmls.n2(λ)
	#MISS Multiply by dX dY / 4 π^2
	ir12(x) = im * k^2 / 2 * (n2^2 - n1^2)
	sr12(x) = 1 / √(n1^2 - x^2)
	Δr12 = rmls.Δz
	it12(x) = im * k^2 / 2 * (n2^2 - n1^2)
	st12(x) = 1 / √(n2^2 - x^2)
	Δt12 = rmls.Δz
	ir21(x) = im * k^2 / 2 * (n1^2 - n2^2)
	sr21(x) = 1 / √(n2^2 - x^2)
	Δr21 = rmls.Δz
	it21(x) = im * k^2 / 2 * (n1^2 - n2^2)
	st21(x) = 1 / √(n1^2 - x^2)
	Δt21 = rmls.Δz
	ScatteringConvolutionCoefficientScalar{T}(ir12, sr12, Δr12, it12, st12, Δt12, ir21, sr21, Δr21, ir21, st21, Δt21, λ, n1, rmls.ref, n2, rmls.ref);
end

function matrix(rmls::RoughInterface{T}, nsx::AbstractRange{<:Real}, nsy::AbstractRange{<:Real}, λ::Real) where T
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

	r12 = zeros(ComplexF64, sizeX, sizeY, sizeX, sizeY)
	t12 = zeros(ComplexF64, sizeX, sizeY, sizeX, sizeY)
	r21 = zeros(ComplexF64, sizeX, sizeY, sizeX, sizeY)
	t21 = zeros(ComplexF64, sizeX, sizeY, sizeX, sizeY)
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
						r21[iX1, iY1, iX2, iY2] -= ri12
						t21[iX1, iY1, iX2, iY2] += ts21
					end
				end
			end
		end
	end
	r12 = reshape(r12, sizeX * sizeY, sizeX * sizeY)
	t12 = reshape(t12, sizeX * sizeY, sizeX * sizeY)
	r21 = reshape(r21, sizeX * sizeY, sizeX * sizeY)
	t21 = reshape(t21, sizeX * sizeY, sizeX * sizeY)
	return (r12, t12, r21, t21)
	return (r12, t12)
end

function PropagationScatteringConvolutionCoefficientScalar(rmls::RoughInterface{T}, λ::Real) where T
	propCoef = PropagationCoefficientScalar([rmls.n1(λ), rmls.n2(λ)], zeros(T,0), λ, rmls.ref)
	scatConvCoef = ScatteringConvolutionCoefficientScalar(rmls, λ)
	return PropagationScatteringConvolutionCoefficientScalar{T}(propCoef, scatConvCoef);
end

coefficientscallar(rmls::RoughInterface, λ) = PropagationScatteringConvolutionCoefficientScalar(rmls, λ)

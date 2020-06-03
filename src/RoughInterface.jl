mutable struct RoughInterface{T<:Real} <: AbstractOpticalComponent{T}
	n1::JolabFunction1D{T,Complex{T}}
	n2::JolabFunction1D{T,Complex{T}}
	Δz::JolabFunction2D{T,T}
	ref::ReferenceFrame{T}
	RoughInterface{T}(n1, n2, Δz, ref) where T = new{T}(n1, n2, Δz, ref)
end
RoughInterface(n1, n2, Δz, ref) = RoughInterface{Float64}(n1, n2, Δz, ref)
rotatestructure!(::PointDetector, ref_ref::ReferenceFrame, θ::Real, ϕ::Real) = rotatereferential!(ref_ref, pointdet.ref, θ, ϕ)

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
	ir12(x) = im * k^2 / 2 * (n2^2 - n1^2)
	sr12(x) = 1 / √(n1^2 - x^2)
	Δr12 = rmls.Δz
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
	fftz_itp = interpolate(fftz, BSpline(Linear()))
 	fftz_itp = scale(fftz_itp, nsxfft, nsyfft)
    fftz_itp = extrapolate(fftz_itp, zero(Complex{T}))

	r12 = Array{Complex{T}, 4}(undef, sizeX, sizeY, sizeX, sizeY)
	t12 = Array{Complex{T}, 4}(undef, sizeX, sizeY, sizeX, sizeY)
	r21 = Array{Complex{T}, 4}(undef, sizeX, sizeY, sizeX, sizeY)
	t21 = Array{Complex{T}, 4}(undef, sizeX, sizeY, sizeX, sizeY)
	@time for iX1 in 1:sizeX
		for iY1 in 1:sizeY
			for iX2 in 1:sizeX
				for iY2 in 1:sizeY
					aux = fftz_itp(nsx[iX2] - nsx[iX1], nsy[iY2] - nsy[iY1])
					r12[iX1, iY1, iX2, iY2] = im * k^2 / 2 * (n2^2 - n1^2) / √(n1^2 - nsx[iX2]^2 - nsy[iY2]^2) * aux
					t12[iX1, iY1, iX2, iY2] = im * k^2 / 2 * (n2^2 - n1^2) / √(n2^2 - nsx[iX2]^2 - nsy[iY2]^2) * aux
					r21[iX1, iY1, iX2, iY2] = im * k^2 / 2 * (n1^2 - n2^2) / √(n2^2 - nsx[iX2]^2 - nsy[iY2]^2) * aux
					t21[iX1, iY1, iX2, iY2] = im * k^2 / 2 * (n1^2 - n2^2) / √(n1^2 - nsx[iX2]^2 - nsy[iY2]^2) * aux
					if (iX1 == iX2 && iY1 == iY2)
						sz1 = √(n1^2 - nsx[iX1]^2 - nsy[iY1]^2)
						sz2 = √(n2^2 - nsx[iX1]^2 - nsy[iY1]^2)
						r12[iX1, iY1, iX2, iY2] += reflectioncoefficientinterfaces(n1, sz1, n2, sz2)
						t12[iX1, iY1, iX2, iY2] += transmissioncoefficientinterfaces(n1, sz1, n2, sz2)
						r21[iX1, iY1, iX2, iY2] += reflectioncoefficientinterfaces(n2, sz2, n1, sz1)
						t21[iX1, iY1, iX2, iY2] += transmissioncoefficientinterfaces(n2, sz2, n1, sz1)
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
end

function PropagationScatteringConvolutionCoefficientScalar(rmls::RoughInterface{T}, λ::Real) where T
	propCoef = PropagationCoefficientScalar([rmls.n1(λ), rmls.n2(λ)], zeros(T,0), λ, rmls.ref)
	scatConvCoef = ScatteringConvolutionCoefficientScalar(rmls, λ)
	return PropagationScatteringConvolutionCoefficientScalar{T}(propCoef, scatConvCoef);
end

coefficientscallar(rmls::RoughInterface, λ) = PropagationScatteringConvolutionCoefficientScalar(rmls, λ)

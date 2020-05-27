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

function PropagationScatteringConvolutionCoefficientScalar(rmls::RoughInterface{T}, λ::Real) where T
	propCoef = PropagationCoefficientScalar([rmls.n1(λ), rmls.n2(λ)], zeros(T,0), λ, rmls.ref)
	scatConvCoef = ScatteringConvolutionCoefficientScalar(rmls, λ)
	return PropagationScatteringConvolutionCoefficientScalar{T}(propCoef, scatConvCoef);
end

coefficientscallar(rmls::RoughInterface, λ) = PropagationScatteringConvolutionCoefficientScalar(rmls, λ)

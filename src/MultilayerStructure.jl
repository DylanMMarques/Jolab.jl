mutable struct MultilayerStructure{T} <: AbstractPropagationComponent{T}
	n_A::Vector{JolabFunction1D{T,Complex{T}}}
	h_A::Vector{T}
	ref::ReferenceFrame{T}
	function MultilayerStructure{T}(n_A, h_A, ref) where T
		length(n_A) != length(h_A) + 2 && error("The length of n_A must be equal to the length of h_A + 2")
		new{T}(n_A, h_A, ref)
	end
end
MultilayerStructure(n_A, h_A, ref) = MultilayerStructure{Float64}(n_A, h_A, ref)

@inline reflectioncoefficientinterfacep(n1::Number, sz1::Number, n2::Number, sz2::Number) = (n2 * sz1 - n1 * sz2) / (n2 * sz1 + n1 * sz2);

@inline reflectioncoefficientinterfaces(n1::Number, sz1::Number, n2::Number, sz2::Number) = (n1 * sz1 - n2 * sz2) / (n1 * sz1 + n2 * sz2);

@inline transmissioncoefficientinterfaces(n1::Number, sz1::Number, n2::Number, sz2::Number) = (2 * n1 * sz1) / (n1 * sz1 + n2 * sz2);

@inline transmissioncoefficientinterfacep(n1::Number, sz1::Number, n2::Number, sz2::Number) = (2 * n1 * sz1) / (n2 * sz1 + n1 * sz2);

function tss₁₂(nsr::Real, n_A::AbstractVector{<:Number}, h_A::AbstractVector{<:Real}, λ::Real)
	sizeA = length(n_A)
	sizeA == length(h_A) + 2 || error("Length n_A must be length h_A + 2")
	sz₂ = √(complex(1 - (nsr / n_A[sizeA])^2))
	sz₁ = √(complex(1 - (nsr / n_A[sizeA-1])^2))
	ri = reflectioncoefficientinterfaces(n_A[sizeA-1], sz₁, n_A[sizeA], sz₂)
	ti = transmissioncoefficientinterfaces(n_A[sizeA-1], sz₁, n_A[sizeA], sz₂)
	@inbounds for iA in (sizeA-2):-1:1
		sz₂ = sz₁
		sz₁ = √(complex(1 - (nsr / n_A[iA])^2))
		propagationTerm = exp(im * 2 * π * n_A[iA+1] / λ * sz₂ * h_A[iA])
		rinterface = reflectioncoefficientinterfaces(n_A[iA], sz₁, n_A[iA+1], sz₂)
		tinterface = transmissioncoefficientinterfaces(n_A[iA], sz₁, n_A[iA+1], sz₂)
		ti = tinterface * ti * propagationTerm / (1 + rinterface * ri * propagationTerm^2)
		ri = (rinterface + ri * propagationTerm^2) / (1 + rinterface * ri * propagationTerm^2)
	end
	return ti
end

function rss₁₂(nsr::Real, n_A::AbstractVector{<:Number}, h_A::AbstractVector{<:Real}, λ::Real)
	sizeA = length(n_A)
	sizeA == length(h_A) + 2 || error("Length n_A must be length h_A + 2")
	sz₂ = √(complex(1 - (nsr / n_A[sizeA])^2))
	sz₁ = √(complex(1 - (nsr / n_A[sizeA-1])^2))
	ri = reflectioncoefficientinterfaces(n_A[sizeA-1], sz₁, n_A[sizeA], sz₂)
	@inbounds for iA in (sizeA-2):-1:1
		sz₂ = sz₁;
		sz₁ = √(complex(1 - (nsr / n_A[iA])^2))
		propagationTerm = exp(im * 2 * π * n_A[iA+1] / λ * sz₂ * h_A[iA])
		rinterface = reflectioncoefficientinterfaces(n_A[iA], sz₁, n_A[iA+1], sz₂)
		ri = (rinterface + ri * propagationTerm^2) / (1 + rinterface * ri * propagationTerm^2)
	end
	return ri
end

function tpp₁₂(nsr::Real, n_A::AbstractVector{<:Number}, h_A::AbstractVector{<:Real}, λ::Real)
	sizeA = length(n_A)
	sizeA == length(h_A) + 2 || error("Length n_A must be length h_A + 2")
	sz₂ = √(complex(1 - (nsr / n_A[sizeA])^2))
	sz₁ = √(complex(1 - (nsr / n_A[sizeA-1])^2))
	ri = reflectioncoefficientinterfacep(n_A[sizeA-1], sz₁, n_A[sizeA], sz₂)
	ti = transmissioncoefficientinterfacep(n_A[sizeA-1], sz₁, n_A[sizeA], sz₂)
	@inbounds for iA in (sizeA-2):-1:1
		sz₂ = sz₁;
		sz₁ = √(complex(1 - (nsr / n_A[iA])^2))
		propagationTerm = exp(im * 2 * π * n_A[iA+1] / λ * sz₂ * h_A[iA])
		rinterface = reflectioncoefficientinterfacep(n_A[iA], sz₁, n_A[iA+1], sz₂)
		tinterface = transmissioncoefficientinterfacep(n_A[iA], sz₁, n_A[iA+1], sz₂)
		ti = tinterface * ti * propagationTerm / (1 + rinterface * ri * propagationTerm^2)
		ri = (rinterface + ri * propagationTerm^2) / (1 + rinterface * ri * propagationTerm^2)
	end
	return ti
end

function rpp₁₂(nsr::Real, n_A::AbstractVector{<:Number}, h_A::AbstractVector{<:Real}, λ::Real)
	sizeA = length(n_A)
	sizeA == length(h_A) + 2 || error("Length n_A must be length h_A + 2")
	sz₂ = √(complex(1 - (nsr / n_A[sizeA])^2))
	sz₁ = √(complex(1 - (nsr / n_A[sizeA-1])^2))
	ri = reflectioncoefficientinterfacep(n_A[sizeA-1], sz₁, n_A[sizeA], sz₂)
	@inbounds for iA in (sizeA-2):-1:1
		sz₂ = sz₁;
		sz₁ = √(complex(1 - (nsr / n_A[iA])^2))
		propagationTerm = exp(im * 2 * π * n_A[iA+1] / λ * sz₂ * h_A[iA])
		rinterface = reflectioncoefficientinterfacep(n_A[iA], sz₁, n_A[iA+1], sz₂)
		ri = (rinterface + ri * propagationTerm^2) / (1 + rinterface * ri * propagationTerm^2)
	end
	return ri
end

function PropagationCoefficientScalar(n_A::AbstractVector{<:Number}, h_A::AbstractVector{T}, λ::Real, ref::ReferenceFrame) where {T<:Real}
	inv_n_A = @view n_A[end:-1:1];
	inv_h_A = @view h_A[end:-1:1];

	r12(nsr) = rss₁₂(nsr, n_A, h_A, λ)
	t12(nsr) = tss₁₂(nsr, n_A, h_A, λ)
	r21(nsr) = rss₁₂(nsr, inv_n_A, inv_h_A, λ)
	t21(nsr) = tss₁₂(nsr, inv_n_A, inv_h_A, λ)

	totalh = sum(h_A);
	ref2 = ref + ReferenceFrame(sin(ref.θ) * cos(ref.ϕ) * totalh, sin(ref.θ) * sin(ref.ϕ) * totalh, cos(ref.θ) * totalh, ref.θ, ref.ϕ);
	return PropagationCoefficientScalar{T}(r12, t12, r21, t21, λ, n_A[1], ref, n_A[end], ref2);
end

function PropagationCoefficientVectorial(n_A::AbstractVector{<:Number}, h_A::AbstractVector{T}, λ::Real, ref::ReferenceFrame) where {T<:Real}
	inv_n_A = @view n_A[end:-1:1];
	inv_h_A = @view h_A[end:-1:1];

	rss12(nsr) = rss₁₂(nsr, n_A, h_A, λ)
	tss12(nsr) = tss₁₂(nsr, n_A, h_A, λ)
	rss21(nsr) = rss₁₂(nsr, inv_n_A, inv_h_A, λ)
	tss21(nsr) = tss₁₂(nsr, inv_n_A, inv_h_A, λ)
	rpp12(nsr) = rpp₁₂(nsr, n_A, h_A, λ)
	tpp12(nsr) = tpp₁₂(nsr, n_A, h_A, λ)
	rpp21(nsr) = rpp₁₂(nsr, inv_n_A, inv_h_A, λ)
	tpp21(nsr) = tpp₁₂(nsr, inv_n_A, inv_h_A, λ)

	totalh = sum(h_A);
	ref2 = ref + ReferenceFrame(sin(ref.θ) * cos(ref.ϕ) * totalh, sin(ref.θ) * sin(ref.ϕ) * totalh, cos(ref.θ) * totalh, ref.θ, ref.ϕ);
	return PropagationCoefficientVectorial{T}(rss12, tss12, rss21, tss21, rpp12, tpp12, rpp21, tpp21, λ, n_A[1], ref, n_A[end], ref2);
end

@inline function PropagationCoefficientScalar(mls::MultilayerStructure, λ::Real)
	n_A = [mls.n_A[i](λ) for i in eachindex(mls.n_A)];
	return PropagationCoefficientScalar(n_A, mls.h_A, λ, mls.ref);
end

@inline function PropagationCoefficientVectorial(mls::MultilayerStructure, λ::Real)
	n_A = [mls.n_A[i](λ) for i in eachindex(mls.n_A)];
	return PropagationCoefficientVectorial(n_A, mls.h_A, λ, mls.ref);
end

function lightinteraction(comps::AbstractVector{<:AbstractPropagationComponent{T}}, angspe::AbstractFieldAngularSpectrum) where T
	vectorial = (size(angspe.e_SXY, 1) == 3)
	if vectorial
		error("not done yet")
	else
		coefs = [coefficientscallar(compi, angspe.λ) for compi in comps]
	end
	coef = mergeorientated_propagationcoefficient(coefs)
	return lightinteraction(coef, angspe)
end

function lightinteraction(mls::MultilayerStructure, angspe::AbstractFieldAngularSpectrum)
	vectorial = size(angspe.e_SXY, 1) == 3;
	if vectorial
		error("not done yet")
	else
		coef = coefficientscallar(mls, angspe.λ)
	end
	return lightinteraction(coef, angspe)
end

coefficientscallar(mls::MultilayerStructure{T}, λ) where T = PropagationCoefficientScalar(mls, λ)
function coefficientscallar(comps::AbstractVector{<:AbstractPropagationComponent{T}}, λ) where T
	coefs = [coefficientscallar(compi, λ) for compi in comps]
	return mergeorientated_propagationcoefficient(coefs)
end

function lightinteraction_recursive(components::AbstractVector{T}, angspe::AbstractFieldAngularSpectrum; nsxOut = angspe.nsx_X::AbstractVector{<:Real}, nsyOut = angspe.nsy_Y::AbstractVector{<:Real}, thresold = 1E-4, sizeM=100::Integer) where T
	vectorial = (size(angspe.e_SXY, 1) == 3)
	vectorial && tobedone()

	coefs = [coefficientscallar(compi, angspe.λ) for compi in components]
	return lightinteraction_recursive(coefs, angspe, nsxOut, nsyOut, thresold = thresold, sizeM = sizeM);
end

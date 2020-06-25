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

function rtss₁₂(nsr::Real, n_A::AbstractVector{<:Number}, h_A::AbstractVector{<:Real}, λ::Real)
	sizeA = length(n_A)
	sizeA == length(h_A) + 2 || error("Length n_A must be length h_A + 2")
	sz₂ = √(complex(1 - (nsr / n_A[sizeA])^2))
	sz₁ = √(complex(1 - (nsr / n_A[sizeA-1])^2))
	ri = reflectioncoefficientinterfaces(n_A[sizeA-1], sz₁, n_A[sizeA], sz₂)
	ti = transmissioncoefficientinterfaces(n_A[sizeA-1], sz₁, n_A[sizeA], sz₂)
	aux = im * 2π / λ
	@inbounds for iA in (sizeA-2):-1:1
		sz₂ = sz₁
		sz₁ = √(complex(1 - (nsr / n_A[iA])^2))
		propagationTerm = exp(aux * n_A[iA+1] * sz₂ * h_A[iA])
		rinterface = reflectioncoefficientinterfaces(n_A[iA], sz₁, n_A[iA+1], sz₂)
		tinterface = transmissioncoefficientinterfaces(n_A[iA], sz₁, n_A[iA+1], sz₂)
		ti = tinterface * ti * propagationTerm / (1 + rinterface * ri * propagationTerm^2)
		ri = (rinterface + ri * propagationTerm^2) / (1 + rinterface * ri * propagationTerm^2)
	end
	return (ri, ti)
end

function rtpp₁₂(nsr::Real, n_A::AbstractVector{<:Number}, h_A::AbstractVector{<:Real}, λ::Real)
	sizeA = length(n_A)
	sizeA == length(h_A) + 2 || error("Length n_A must be length h_A + 2")
	sz₂ = √(complex(1 - (nsr / n_A[sizeA])^2))
	sz₁ = √(complex(1 - (nsr / n_A[sizeA-1])^2))
	ri = reflectioncoefficientinterfacep(n_A[sizeA-1], sz₁, n_A[sizeA], sz₂)
	ti = transmissioncoefficientinterfacep(n_A[sizeA-1], sz₁, n_A[sizeA], sz₂)
	aux = im * 2π / λ
	@inbounds for iA in (sizeA-2):-1:1
		sz₂ = sz₁;
		sz₁ = √(complex(1 - (nsr / n_A[iA])^2))
		propagationTerm = exp(aux * n_A[iA+1] * sz₂ * h_A[iA])
		rinterface = reflectioncoefficientinterfacep(n_A[iA], sz₁, n_A[iA+1], sz₂)
		tinterface = transmissioncoefficientinterfacep(n_A[iA], sz₁, n_A[iA+1], sz₂)
		ti = tinterface * ti * propagationTerm / (1 + rinterface * ri * propagationTerm^2)
		ri = (rinterface + ri * propagationTerm^2) / (1 + rinterface * ri * propagationTerm^2)
	end
	return (ri, ti)
end

function coefficient_geral(mls::MultilayerStructure{T}, fieldi::FieldAngularSpectrum) where {T<:Real}
	n_A = [mls.n_A[i](λ) for i in eachindex(mls.n_A)]

	inv_n_A = n_A[end:-1:1]
	inv_h_A = mls.h_A[end:-1:1]

	sizeXY = length(fieldi.nsx_X) * length(fieldi.nsy_Y)
	r12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	i = 1
	@inbounds for iY in eachindex(fieldi.nsy_Y)
		for iX in eachindex(fieldi.nsx_X)
			nsr = √(fieldi.nsx_X[iX]^2 + fieldi.nsy_Y[iY]^2)
			(r12.diag[i], t12.diag[i]) = rtss₁₂(nsr, n_A, mls.h_A, fieldi.λ)
			(r21.diag[i], t21.diag[i]) = rtss₁₂(nsr, inv_n_A, inv_h_A, fieldi.λ)
			i += 1
		end
	end
	fieldl = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, first(n_A), -1, mls.ref)
	fieldr = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, last(n_A), 1, ref2(mls))
	return ScaterringMatrix{T, Diagonal{Complex{T},Vector{Complex{T}}}, typeof(fieldl), typeof(fieldr)}(r12, t12, r21, t21, fieldl, fieldr)
end
@inline coefficient_specific(mls::MultilayerStructure, nsx_X, nsy_Y, λ) = coefficient_geral(mls, nsx_X, nsy_Y, λ)

@inline function ref2(mls::MultilayerStructure)
	totalh = sum(mls.h_A);
	return mls.ref + ReferenceFrame(sin(mls.ref.θ) * cos(mls.ref.ϕ) * totalh, sin(mls.ref.θ) * sin(mls.ref.ϕ) * totalh, cos(mls.ref.θ) * totalh, mls.ref.θ, mls.ref.ϕ);
end

@inline ref1(mls::MultilayerStructure) = mls.ref

function lightinteraction(comps::AbstractVector{<:AbstractPropagationComponent{T}}, angspe::AbstractFieldAngularSpectrum) where T
	vectorial = (size(angspe.e_SXY, 1) == 3)
	if vectorial
		error("not done yet")
	else
		coefs = [coefficient_geral(compi, angspe) for compi in comps]
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

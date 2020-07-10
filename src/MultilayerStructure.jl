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

function coefficient_general(mls::MultilayerStructure{T}, fieldi::FieldAngularSpectrum) where {T<:Real}
	isapprox(fieldi.n, fieldi.dir > 0 ? mls.n_A[1](fieldi.λ) : mls.n_A[end](fieldi.λ), atol = @tol) || error("Field medium and multilayer structure medium are different")
	checkorientation(fieldi.ref, mls.ref) || errorToDo()

	n_A = [mls.n_A[i](fieldi.λ) for i in eachindex(mls.n_A)]
	inv_n_A = n_A[end:-1:1]
	inv_h_A = mls.h_A[end:-1:1]

	sizeXY = length(fieldi.nsx_X) * length(fieldi.nsy_Y)
	r12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	sizeX = length(fieldi.nsx_X)
	sizeY = length(fieldi.nsy_Y)
	cart = LinearIndices((sizeX, sizeY))
	Threads.@threads for iY in eachindex(fieldi.nsy_Y)
		for iX in eachindex(fieldi.nsx_X)
			nsr = √(fieldi.nsx_X[iX]^2 + fieldi.nsy_Y[iY]^2)
			i = cart[iX, iY]
			(r12.diag[i], t12.diag[i]) = rtss₁₂(nsr, n_A, mls.h_A, fieldi.λ)
			(r21.diag[i], t21.diag[i]) = rtss₁₂(nsr, inv_n_A, inv_h_A, fieldi.λ)
		end
	end

	ref = (fieldi.dir > 0 ? ref1(mls) : ref2(mls))
	if !checkposition(fieldi.ref, ref)
		propM = propagationmatrix(fieldi, ref)
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
		fieldl = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, first(n_A), -1, fieldi.ref)
		fieldr = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, last(n_A), 1, ref2(mls))
	else
		fieldl = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, first(n_A), -1, mls.ref)
		fieldr = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, last(n_A), 1, fieldi.ref)
	end
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), Diagonal{Complex{T},Vector{Complex{T}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

@inline function ref2(mls::MultilayerStructure)
	totalh = sum(mls.h_A);
	return mls.ref + ReferenceFrame(sin(mls.ref.θ) * cos(mls.ref.ϕ) * totalh, sin(mls.ref.θ) * sin(mls.ref.ϕ) * totalh, cos(mls.ref.θ) * totalh, mls.ref.θ, mls.ref.ϕ);
end

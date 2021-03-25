mutable struct MultilayerStructure{T} <: AbstractPropagationComponent{T}
	n_A::Vector{Material{T}}
	h_A::Vector{T}
	ref::ReferenceFrame{T}
	function MultilayerStructure{T}(n_A, h_A, ref) where T
		length(n_A) != length(h_A) + 2 && error("The length of n_A must be equal to the length of h_A + 2")
		new{T}(n_A, h_A, ref)
	end
end

n1(mls::MultilayerStructure, λ) = first(mls.n_A)(λ)
n2(mls::MultilayerStructure, λ) = last(mls.n_A)(λ)

n(mls::MultilayerStructure, λ) = [ni(λ) for ni in mls.n_A]

"""
	MultilayerStructure(T, n, h, ref)

Initializes a multilayer structure.
- `T` number type specifing data precision. Ex: `Float32`, `BigFloat`, ... The default is Float64;
- `n` - vector specifying the refractive index of each layer; See how to define a refractive index.
- `h` - vector specifying the thickness of each layer. The layer with refractive index `n[i]` has a thickness of `h[i-1]` meters;
- `ref` - `ReferenceFrame` specifying the multilayer position and orientation. The first interface of the multilayer structure is intersect the reference frame.

**Examples:**

```julia
mls = MultilayerStructure([1, 2, 3], [100E-9], ReferenceFrame(0,0,0.,0,0))
```

```julia
n_1(λ) = √(complex(8.393 + .14383 / ((λ*1E6)^2 - 0.2421^2) + 4430.99 / ((λ*1E6)^2 - 36.71^2)))
mls = MultilayerStructure([1, n_1], [], ReferenceFrame(0,0,0.,0,0))
```
"""
MultilayerStructure(::Type{T}, n_A, h_A, ref) where T = MultilayerStructure{T}(n_A, h_A, ref)
MultilayerStructure(n_A, h_A, ref) = MultilayerStructure{Float64}(n_A, h_A, ref)

reflectioncoefficientinterfacep(n1, sz1, n2, sz2) = (n2 * sz1 - n1 * sz2) / (n2 * sz1 + n1 * sz2);

reflectioncoefficientinterfaces(n1, sz1, n2, sz2) = (n1 * sz1 - n2 * sz2) / (n1 * sz1 + n2 * sz2);

transmissioncoefficientinterfaces(n1, sz1, n2, sz2) = (2 * n1 * sz1) / (n1 * sz1 + n2 * sz2);

transmissioncoefficientinterfacep(n1, sz1, n2, sz2) = (2 * n1 * sz1) / (n2 * sz1 + n1 * sz2);

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

"""
     coefficient_general(::MultilayerStructure, ::FieldAngularSpectrumScalar)

Calculates the scattering matrix of a multilayer structure for an incident angular spectrum
- **Type:** Transmission and reflection matrices are diagonal
- **Time:** very short; scales with Nx Ny
- **RAM:** very small; scales with Nx Ny
- **Convergence** sampling of nsx and nsy
"""
function coefficient_general(mls::MultilayerStructure, fieldi::FieldAngularSpectrumScalar)
	checkapplicability(mls, fieldi)

	n_A = n(mls, fieldi.λ)
	inv_n_A = n_A[end:-1:1]
	inv_h_A = mls.h_A[end:-1:1]

	scat = get_scatteringmatrixtype(mls, fieldi)
	cart = CartesianIndices(fieldi)
	@inbounds Threads.@threads for i in iterator_index(fieldi)
		nsr = √(fieldi.nsx_X[cart[i][2]]^2 + fieldi.nsy_Y[cart[i][3]]^2)
		(scat.r₁₂.diag[i], scat.t₁₂.diag[i]) = rtss₁₂(nsr, n_A, mls.h_A, fieldi.λ)
		(scat.r₂₁.diag[i], scat.t₂₁.diag[i]) = rtss₁₂(nsr, inv_n_A, inv_h_A, fieldi.λ)
	end
	correctscatteringmatrix_referenceframes!(scat, mls, fieldi)
	return scat
end

"""
    lightinteraction(::MultilayerStructure, ::FieldAngularSpectrumScalar)

Calculates the reflected and transmitted fields from a multilayer structure for an incident angular spectrum
- **Time:** very short; scales with Nx Ny
- **RAM:** None
- **Convergence** sampling of nsx and nsy
"""
function lightinteraction(mls::MultilayerStructure, fieldi::FieldAngularSpectrumScalar)
	checkapplicability(mls, fieldi)
	(fieldl, fieldr) = getfields_lr(mls, fieldi)
	fieldi_newref = changereferenceframe(fieldi, dir(fieldi) > 0 ? ref1(mls) : ref2(mls))

	n_A = n(mls, fieldi.λ)
	inv_n_A = n_A[end:-1:1]
	inv_h_A = mls.h_A[end:-1:1]

	cart = CartesianIndices(fieldi)

	@inbounds Threads.@threads for i in iterator_index(fieldi)
		nsr = √(fieldi_newref.nsx_X[cart[i][2]]^2 + fieldi_newref.nsy_Y[cart[i][3]]^2)
		if dir(fieldi_newref) > 0
			(r, t) = rtss₁₂(nsr, n_A, mls.h_A, fieldi.λ)
			fieldl.e_SXY[i] = fieldi_newref.e_SXY[i] * r
			fieldr.e_SXY[i] = fieldi_newref.e_SXY[i] * t
		else
			(r, t) = rtss₁₂(nsr, inv_n_A, inv_h_A, fieldi.λ)
			fieldl.e_SXY[i] = fieldi_newref.e_SXY[i] * t
			fieldr.e_SXY[i] = fieldi_newref.e_SXY[i] * r
		end
	end
	return (fieldl, fieldr)
end

function get_scatteringmatrixtype(mls::MultilayerStructure{T}, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,D,X,B}
	sizeXY = length(fieldi.e_SXY)
	r12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))

	fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, n1(mls, fieldi.λ), dir(fieldi) > 0 ? fieldi.ref : ref1(mls))
	fieldr = FieldAngularSpectrumScalar{T,1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, n2(mls, fieldi.λ), dir(fieldi) > 0 ? ref2(mls) : fieldi.ref)

	return ScatteringMatrix{T, FieldAngularSpectrumScalar{T,-1,X,B}, FieldAngularSpectrumScalar{T,1,X,B}, Diagonal{Complex{T},Vector{Complex{T}}}, Diagonal{Complex{T},Vector{Complex{T}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

function getfields_lr(mls::MultilayerStructure{T}, fieldi::FieldAngularSpectrumScalar{T,A,X,B}) where {T,X,A,B}
	fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, n1(mls, fieldi.λ), ref1(mls))
	fieldr = FieldAngularSpectrumScalar{T,1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, n2(mls, fieldi.λ), ref2(mls))
	return (fieldl, fieldr)
end

@inline function ref2(mls::MultilayerStructure)
	totalh = sum(mls.h_A);
	return mls.ref + ReferenceFrame(sin(mls.ref.θ) * cos(mls.ref.ϕ) * totalh, sin(mls.ref.θ) * sin(mls.ref.ϕ) * totalh, cos(mls.ref.θ) * totalh, mls.ref.θ, mls.ref.ϕ);
end

function checkapplicability(mls::MultilayerStructure, fieldi::AbstractFieldAngularSpectrum)
	isapprox(fieldi.n, dir(fieldi) > 0 ? n1(mls, fieldi.λ) : n2(mls, fieldi.λ), atol = @tol) || error("Field medium and multilayer structure medium are different")
	checkorientation(fieldi.ref, mls.ref) || errorToDo()
end

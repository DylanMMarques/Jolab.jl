mutable struct FieldAngularSpectrumScalarRadialSymmetric{T,D,X<:AbstractVector{T},Y<:AbstractVector{Complex{T}}} <: AbstractFieldAngularSpectrum{T,D}
	nsr_R::X
	e_SXY::Y
	λ::T
	n::Complex{T}
	ref::ReferenceFrame{T}
	function FieldAngularSpectrumScalarRadialSymmetric{T,D,X,Y}(nsr_R, e_SXY, λ, n, ref) where {T,D,X,Y}
		length(nsr_R) == length(e_SXY) || error("The length of nsr_R must be the same as the length of e_SXY.")
		(abs(ref.x) > 1E-10 || abs(ref.y) > 1E-10 || abs(ref.θ) > 1E-5 || abs(ref.ϕ) > 1E-5) && error("The referenceframe of a symetric field must have x, y, θ and ϕ equal to 0")
		return new{T,D,X,Y}(nsr_R, e_SXY, λ, n, ref);
	end
	function FieldAngularSpectrumScalarRadialSymmetric{T,D}(nsr_R::X, e_SXY::Y, λ, n, ref) where {T,D,X,Y}
		return FieldAngularSpectrumScalarRadialSymmetric{T,D,X,Y}(nsr_R, e_SXY, λ, n, ref);
	end
end

function FieldAngularSpectrumScalarRadialSymmetric_gaussian(::Type{T}, nsr_R, ω, λ, n, dir, ref) where T
	k = 2π / T(λ);
	norm = 	T(ω) * √(T(1 / 32 / π^3));
	e_SXY = Vector{Complex{T}}(undef, length(nsr_R))
	angspe = FieldAngularSpectrumScalarRadialSymmetric{T,dir}(nsr_R, e_SXY, λ, n, ref);
	@inbounds @simd for i in iterator_index(angspe)
		angspe.e_SXY[i] = norm * exp((-T(ω)^2 * k^2 / 16) * nsr_R[i]^2)
	end
	return angspe
end
FieldAngularSpectrumScalarRadialSymmetric_gaussian(arg...) =  FieldAngularSpectrumScalarRadialSymmetric_gaussian(Float64, arg...)

function FieldAngularSpectrumScalarRadialSymmetric_uniform(::Type{T}, nsr_R, λ, n, dir, ref) where T
	e_SXY = ones(Complex{T}, length(nsr_R))
	angspe = FieldAngularSpectrumScalarRadialSymmetric{T,dir}(nsr_R, e_SXY, λ, n, ref);
	angspe.e_SXY ./= √intensity(angspe)
	return angspe
end
FieldAngularSpectrumScalarRadialSymmetric_uniform(arg...) =  FieldAngularSpectrumScalarRadialSymmetric_uniform(Float64, arg...)

function changereferenceframe!(angspe::FieldAngularSpectrumScalarRadialSymmetric, refnew::ReferenceFrame)
	(abs(angspe.ref.x - refnew.x) > 1E-10 || abs(angspe.ref.y - refnew.y) > 1E-10 || abs(angspe.ref.θ - refnew.θ) > 1E-5 || abs(angspe.ref.ϕ - refnew.ϕ) > 1E-5) && error("FieldAngularSpectrumScalarRadialSymmetric cannot be translate in x or y or rotated")
	checkposition(angspe.ref, refnew) || translatereferenceframe!(angspe, refnew);
end

function add!(fielda::FieldAngularSpectrumScalarRadialSymmetric, fieldb::FieldAngularSpectrumScalarRadialSymmetric)
	samedefinitions(fielda, fieldb) || error("Cannot sum the fields. Different definitions")
	@inbounds @simd for i in iterator_index(fielda)
		fielda.e_SXY[i] += fieldb.e_SXY[i]
	end
end

function translatereferenceframe!(angspe::FieldAngularSpectrumScalarRadialSymmetric, refnew::ReferenceFrame)
	refΔz = refnew.z - angspe.ref.z;
	kim = im * 2 * π / angspe.λ;
	@inbounds @simd for i in iterator_index(angspe)
		nsz_a = dir(angspe) * nsz(angspe.n, angspe.nsr_R[i], 0)
		angspe.e_SXY[i] *= exp(kim * nsz_a * refΔz)
	end
	angspe.ref.z = refnew.z;
end

function propagationmatrix(fieldl::FieldAngularSpectrumScalarRadialSymmetric{T}, ref::ReferenceFrame) where T
	checkorientation(fieldl.ref, ref) || error("Cannot calculate propagation matrix as the referenceframe are not oriented")
	(abs(ref.x - 0) < 1E-12 && abs(ref.y - 0) < 1E-12) || error("Cannot propagate a radial symmetric angular spectrum to out of optical axis")
	refΔz = ref.z - fieldl.ref.z;

	imk = im * 2π / fieldl.λ;
	propMatrix = Diagonal(Vector{Complex{T}}(undef, length(fieldl.e_SXY)))
	@inbounds for i in iterator_index(fieldl)
		nsz_a = dir(fieldl) * nsz(fieldl.n, fieldl.nsr_R[i], 0)
		propMatrix.diag[i] = exp(imk * nsz_a * refΔz)
	end
	return propMatrix
end

function intensity(angspe::FieldAngularSpectrumScalarRadialSymmetric{T}) where T
	int = zero(T)
	@inbounds @simd for i in iterator_index(angspe)
		Δnsr = Δvector(angspe.nsr_R, i)
		int += abs2(angspe.e_SXY[i]) * Δnsr * angspe.nsr_R[i]
	end
	return 32π^5 / angspe.λ^2 * real(angspe.n) * int
end

CartesianIndices(field::FieldAngularSpectrumScalarRadialSymmetric) = Base.CartesianIndices((1, length(field.nsr_R), 1))
LinearIndices(field::FieldAngularSpectrumScalarRadialSymmetric) = Base.LinearIndices((1, length(field.nsr_R), 1))

function Base.:copy(field::FieldAngularSpectrumScalarRadialSymmetric{T,D,X,Y}) where {T,D,X,Y}
	return FieldAngularSpectrumScalarRadialSymmetric{T,D,X,Y}(deepcopy(field.nsr_R), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

function copy_differentD(field::FieldAngularSpectrumScalarRadialSymmetric{T,D,X,Y}) where {T,D,X,Y}
	return FieldAngularSpectrumScalarRadialSymmetric{T,-D,X,Y}(deepcopy(field.nsr_R), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

function samedefinitions(fieldl::L, fieldr::R) where {L <: FieldAngularSpectrumScalarRadialSymmetric, R <: FieldAngularSpectrumScalarRadialSymmetric}
	isapprox(fieldl.nsr_R, fieldr.nsr_R, atol = @tol) || error("nsr_R are different")
	isapprox(fieldl.n, fieldr.n, atol = @tol) || error("n are different")
	isapprox(fieldl.λ, fieldr.λ, atol = @tol) || error("λ is different")
	fieldl.ref == fieldr.ref || error("reference frames are different")
	return true
end

mutable struct FieldAngularSpectrumScalar{T,D, X<:AbstractVector{T}, Y<:AbstractVector{Complex{T}}} <: AbstractFieldAngularSpectrum{T,D}
	nsx_X::X
	nsy_Y::X
	e_SXY::Y
	λ::T
	n::Complex{T}
	ref::ReferenceFrame{T}
	function FieldAngularSpectrumScalar{T,D,X,Y}(nsx_X, nsy_Y, e_SXY, λ, n, ref) where {T,D,X,Y}
		length(nsx_X) * length(nsy_Y) == length(e_SXY) || error("length(nsx_X) * length(nsy_Y) must be equal to length(e_SXY)")
		return new{T,D,X,Y}(nsx_X, nsy_Y, e_SXY, λ, n, ref);
	end
	function FieldAngularSpectrumScalar{T,D}(nsx_X::X1, nsy_Y::X2, e_SXY::Y, λ, n, ref) where {T,D,X1,X2,Y}
		length(nsx_X) * length(nsy_Y) == length(e_SXY) || error("length(nsx_X) * length(nsy_Y) must be equal to length(e_SXY)")
		X = promote_type(X1,X2)
		return new{T,D,X,Y}(nsx_X, nsy_Y, e_SXY, λ, n, ref);
	end
	function FieldAngularSpectrumScalar{T,D,X}(nsx_X, nsy_Y, e_SXY::Y, λ, n, ref) where {T,D,X,Y}
		length(nsx_X) * length(nsy_Y) == length(e_SXY) || error("length(nsx_X) * length(nsy_Y) must be equal to length(e_SXY)")
		return new{T,D,X,Y}(nsx_X, nsy_Y, e_SXY, λ, n, ref);
	end
end

function FieldAngularSpectrumScalar_uniform(nsx_X, nsy_Y, λ, n, dir, ref; T = Float64)
	e_SXY = ones(Complex{T}, length(nsx_X) * length(nsy_Y))
	angspe = FieldAngularSpectrumScalar{T,dir}(nsx_X, nsy_Y, e_SXY, λ, n, ref);
	int = √intensity(angspe)
	@inbounds @simd for i in iterator_index(angspe)
		angspe.e_SXY[i] /= int
	end
	return angspe
end

function FieldAngularSpectrumScalar_gaussian(nsx_X, nsy_Y, ω, λ, n, dir, ref; T = Float64)
	k = 2π / λ;
	norm = 	ω * √(1 / 32 / π^3);
	e_SXY = Vector{Complex{T}}(undef, length(nsx_X) * length(nsy_Y))
	angspe = FieldAngularSpectrumScalar{T,dir}(nsx_X, nsy_Y, e_SXY, λ, n, ref);
	cart = CartesianIndices(angspe)
	@inbounds @simd for i in iterator_index(angspe)
		angspe.e_SXY[i] = norm * exp((-ω^2 * k^2 / 16) * (nsx_X[cart[i][2]]^2 + nsy_Y[cart[i][3]]^2))
	end
	return angspe
end

function changereferenceframe!(angspe::AbstractFieldAngularSpectrum, refnew::ReferenceFrame)
	#Needs checking after doing the 2D interpolation
	checkposition(angspe.ref, refnew) || translatereferenceframe!(angspe, refnew);
	checkorientation(angspe.ref, refnew) || errorToDo() #rotatereferenceframe!(angspe, refnew);
end

function translatereferenceframe!(angspe::FieldAngularSpectrumScalar, refnew::ReferenceFrame)
	(refΔx, refΔy, refΔz) = (refnew.x - angspe.ref.x, refnew.y - angspe.ref.y, refnew.z - angspe.ref.z)
	# Can't use based on referenceframe as the rotation is inverted
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, angspe.ref.θ, angspe.ref.ϕ);
	kim = im * 2π / angspe.λ;
	cart = CartesianIndices(angspe)
	@inbounds @simd for i in iterator_index(angspe)
		nsz = dir(angspe) * √(angspe.n^2 - angspe.nsx_X[cart[i][2]]^2 - angspe.nsy_Y[cart[i][3]]^2)
		angspe.e_SXY[i] *= exp(kim * (angspe.nsx_X[cart[i][2]] * refΔx + angspe.nsy_Y[cart[i][3]] * refΔy + nsz * refΔz))
	end
	(angspe.ref.x, angspe.ref.y, angspe.ref.z) = (refnew.x, refnew.y, refnew.z)
end

function intensity(angspe::FieldAngularSpectrumScalar{T}) where T
	int = zero(T)
	cart = CartesianIndices(angspe)
	@inbounds @simd for i in iterator_index(angspe)
		Δx = Δvector(angspe.nsx_X, cart[i][2])
		Δy = Δvector(angspe.nsy_Y, cart[i][3])
		int += abs2(angspe.e_SXY[i]) * Δx * Δy
	end
	return 16π^4 / angspe.λ^2 * real(angspe.n) * int
end

function samedefinitions(fieldl::L, fieldr::R) where {L <: FieldAngularSpectrumScalar, R <: FieldAngularSpectrumScalar}
	isapprox(fieldl.nsx_X, fieldr.nsx_X, atol = @tol) || error("nsx_X are different")
	isapprox(fieldl.nsy_Y, fieldr.nsy_Y, atol = @tol) || error("nsy_Y are different")
	isapprox(fieldl.n, fieldr.n, atol = @tol) || error("n are different")
	isapprox(fieldl.λ, fieldr.λ, atol = @tol) || error("λ is different")
	fieldl.ref == fieldr.ref || error("reference frames are different")
	return true
end

function add!(fielda::FieldAngularSpectrumScalar{T}, fieldb::FieldAngularSpectrumScalar{T}) where T
	samedefinitions(fielda, fieldb) || error("Cannot sum the fields. Different definitions")
	@inbounds @simd for i in iterator_index(fielda)
		fielda.e_SXY[i] += fieldb.e_SXY[i]
	end
end

function Base.:copy(field::FieldAngularSpectrumScalar{T,D,X}) where {T,D,X}
	return FieldAngularSpectrumScalar{T,D,X}(deepcopy(field.nsx_X), deepcopy(field.nsy_Y), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

function copy_differentD(field::FieldAngularSpectrumScalar{T,1,X}) where {T,X}
	return FieldAngularSpectrumScalar{T,-1,X}(deepcopy(field.nsx_X), deepcopy(field.nsy_Y), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

function copy_differentD(field::FieldAngularSpectrumScalar{T,-1,X}) where {T,X}
	return FieldAngularSpectrumScalar{T,1,X}(deepcopy(field.nsx_X), deepcopy(field.nsy_Y), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

CartesianIndices(field::FieldAngularSpectrumScalar) = Base.CartesianIndices((1, length(field.nsx_X), length(field.nsy_Y)))

function propagationmatrix(fieldl::AbstractFieldAngularSpectrum{T}, ref::ReferenceFrame) where T
	checkorientation(fieldl.ref, ref) || error("Cannot calculate propagation matrix as the referenceframe are not oriented")
	refΔx = ref.x - fieldl.ref.x;
	refΔy = ref.y - fieldl.ref.y;
	refΔz = ref.z - fieldl.ref.z;

	# Can't use based on referenceframe as the rotation is inverted
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, fieldl.ref.θ, fieldl.ref.ϕ);
	imk = im * 2π / fieldl.λ;
	propMatrix = Diagonal(Vector{Complex{T}}(undef, length(fieldl.e_SXY)))
	cart = CartesianIndices(fieldl)
	@inbounds for i in iterator_index(fieldl)
		nsz = dir(fieldl) * √(complex(fieldl.n^2 - fieldl.nsx_X[cart[i][2]]^2 - fieldl.nsy_Y[cart[i][3]]^2))
		propMatrix.diag[i] = exp(imk * (fieldl.nsx_X[cart[i][2]] * refΔx + fieldl.nsy_Y[cart[i][3]] * refΔy + nsz * refΔz))
	end
	return propMatrix
end

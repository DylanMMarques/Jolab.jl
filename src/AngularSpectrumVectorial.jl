mutable struct FieldAngularSpectrumVectorial{T,D, X<:AbstractVector{T}, Y<:AbstractVector{Complex{T}}} <: AbstractFieldAngularSpectrum{T,D}
	nsx_X::X
	nsy_Y::X
	e_SXY::Y
	λ::T
	n::Complex{T}
	ref::ReferenceFrame{T}
	function FieldAngularSpectrumVectorial{T,D,X,Y}(nsx_X, nsy_Y, e_SXY, λ, n, ref) where {T,D,X,Y}
		length(nsx_X) * length(nsy_Y) * 2 == length(e_SXY) || error("length(nsx_X) * length(nsy_Y) must be equal to length(e_SXY)")
		return new{T,D,X,Y}(nsx_X, nsy_Y, e_SXY, λ, n, ref);
	end
	function FieldAngularSpectrumVectorial{T,D}(nsx_X::X1, nsy_Y::X2, e_SXY::Y, λ, n, ref) where {T,D,X1,X2,Y}
		length(nsx_X) * length(nsy_Y) * 2 == length(e_SXY) || error("length(nsx_X) * length(nsy_Y) must be equal to length(e_SXY)")
		X = promote_type(X1,X2)
		return new{T,D,X,Y}(nsx_X, nsy_Y, e_SXY, λ, n, ref);
	end
	function FieldAngularSpectrumVectorial{T,D,X}(nsx_X, nsy_Y, e_SXY::Y, λ, n, ref) where {T,D,X,Y}
		length(nsx_X) * length(nsy_Y) * 2 == length(e_SXY) || error("length(nsx_X) * length(nsy_Y) must be equal to length(e_SXY)")
		return new{T,D,X,Y}(nsx_X, nsy_Y, e_SXY, λ, n, ref);
	end
end

function FieldAngularSpectrumVectorial_uniform(::Type{T}, nsx_X, nsy_Y, λ, n, dir, ref) where T
	e_SXY = ones(Complex{T}, 2 * length(nsx_X) * length(nsy_Y))
	angspe = FieldAngularSpectrumScalar{T,dir}(nsx_X, nsy_Y, e_SXY, λ, n, ref);
	angspe.e_SXY ./= √intensity(angspe)
	return angspe
end
FieldAngularSpectrumVectorial_uniform(arg...) = FieldAngularSpectrumVectorial_uniform(Float64, arg...)

function FieldAngularSpectrumVectorial_gaussian(::Type{T}, nsx_X, nsy_Y, ω, λ, n, dir, ref) where T
	k = 2π / T(λ);
	norm = 	T(ω) * √(T(1 / 32 / π^3));
	e_SXY = Vector{Complex{T}}(undef, 2 * length(nsx_X) * length(nsy_Y))
	angspe = FieldAngularSpectrumScalar{T,dir}(nsx_X, nsy_Y, e_SXY, λ, n, ref);
	cart = CartesianIndices(angspe)
	@inbounds @simd for i in iterator_index(angspe)
		angspe.e_SXY[i] = norm * exp((-T(ω)^2 * k^2 / 16) * (nsx_X[cart[i][2]]^2 + nsy_Y[cart[i][3]]^2))
	end
	return angspe
end
FieldAngularSpectrumVectorial_gaussian(arg...) = FieldAngularSpectrumVectorial_gaussian(Float64, arg...)

function translatereferenceframe!(angspe::FieldAngularSpectrumVectorial, refnew::ReferenceFrame)
	(refΔx, refΔy, refΔz) = (refnew.x - angspe.ref.x, refnew.y - angspe.ref.y, refnew.z - angspe.ref.z)
	# Can't use based on referenceframe as the rotation is inverted
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, angspe.ref.θ, angspe.ref.ϕ);
	kim = im * 2π / angspe.λ;
	cart = CartesianIndices(angspe)
	@inbounds @simd for i in iterator_index(angspe)
		nsz_a = dir(angspe) * nsz(angspe.n, angspe.nsx_X[cart[i][2]], angspe.nsy_Y[cart[i][3]])
		angspe.e_SXY[i] *= exp(kim * (angspe.nsx_X[cart[i][2]] * refΔx + angspe.nsy_Y[cart[i][3]] * refΔy + nsz_a * refΔz))
	end
	(angspe.ref.x, angspe.ref.y, angspe.ref.z) = (refnew.x, refnew.y, refnew.z)
end

function rotatereferenceframe!(angspe::FieldAngularSpectrumVectorial{T}, refnew::ReferenceFrame) where T
	errorToDo()
	sizeX, sizeY = length(angspe.nsx_X), length(angspe.nsy_Y)
	nsz_1 = dir(angspe) * T(nsz(angspe.n, first(angspe.nsx_X), first(angspe.nsy_Y)))
	nsz_2 = dir(angspe) * T(nsz(angspe.n, first(angspe.nsx_X), last(angspe.nsy_Y)))
	nsz_3 = dir(angspe) * T(nsz(angspe.n, last(angspe.nsx_X), first(angspe.nsy_Y)))
	nsz_4 = dir(angspe) * T(nsz(angspe.n, last(angspe.nsx_X), last(angspe.nsy_Y)))
	(nsx_1, nsy_1, ~) = rotatecoordinatesfromto(first(angspe.nsx_X), first(angspe.nsy_Y), nsz_1, angspe.ref, refnew)
	(nsx_2, nsy_2, ~) = rotatecoordinatesfromto(first(angspe.nsx_X), last(angspe.nsy_Y), nsz_2, angspe.ref, refnew)
	(nsx_3, nsy_3, ~) = rotatecoordinatesfromto(last(angspe.nsx_X), first(angspe.nsy_Y), nsz_3, angspe.ref, refnew)
	(nsx_4, nsy_4, ~) = rotatecoordinatesfromto(last(angspe.nsx_X), last(angspe.nsy_Y), nsz_4, angspe.ref, refnew)
	nsx_r_X = range(min(nsx_1, nsx_2, nsx_3, nsx_4), max(nsx_1, nsx_2, nsx_3, nsx_4), length = sizeX)
	nsy_r_Y = range(min(nsy_1, nsy_2, nsy_3, nsy_4), max(nsy_1, nsy_2, nsy_3, nsy_4), length = sizeY)

	cart = CartesianIndices(angspe)
	e_itp = LinearInterpolation((angspe.nsx_X, angspe.nsy_Y), reshape(angspe.e_SXY, sizeX, sizeY), extrapolation_bc = zero(Complex{T}))
	@inbounds @simd for i in iterator_index(angspe)
		nsx, nsy = nsx_r_X[cart[i][2]], nsy_r_Y[cart[i][3]]
		(nsx_old, nsy_old, ~) = rotatecoordinatesfromto(nsx, nsy, dir(angspe) * T(nsz(angspe.n, nsx, nsy)), refnew, angspe.ref)
		angspe.e_SXY[i] = e_itp(nsx_old, nsy_old)
	end
	angspe.nsx_X, angspe.nsy_Y = nsx_r_X, nsy_r_Y
	angspe.ref = refnew
end

function intensity(angspe::FieldAngularSpectrumVectorial{T}) where T
	int = zero(T)
	cart = CartesianIndices(angspe)
	@inbounds @simd for i in iterator_index(angspe)
		Δnsx = Δvector(angspe.nsx_X, cart[i][2])
		Δnsy = Δvector(angspe.nsy_Y, cart[i][3])
		int += abs2(angspe.e_SXY[i]) * Δnsx * Δnsy
	end
	return 16π^4 / angspe.λ^2 * real(angspe.n) * int
end

function samedefinitions(fieldl::L, fieldr::R) where {L <: FieldAngularSpectrumVectorial, R <: FieldAngularSpectrumVectorial}
	isapprox(fieldl.nsx_X, fieldr.nsx_X, atol = @tol) || error("nsx_X are different")
	isapprox(fieldl.nsy_Y, fieldr.nsy_Y, atol = @tol) || error("nsy_Y are different")
	isapprox(fieldl.n, fieldr.n, atol = @tol) || error("n are different")
	isapprox(fieldl.λ, fieldr.λ, atol = @tol) || error("λ is different")
	fieldl.ref == fieldr.ref || error("reference frames are different")
	return true
end

function add!(fielda::FieldAngularSpectrumVectorial{T}, fieldb::FieldAngularSpectrumVectorial{T}) where T
	samedefinitions(fielda, fieldb) || error("Cannot sum the fields. Different definitions")
	@inbounds @simd for i in iterator_index(fielda)
		fielda.e_SXY[i] += fieldb.e_SXY[i]
	end
end

function Base.:copy(field::FieldAngularSpectrumVectorial{T,D,X,Y}) where {T,D,X,Y}
	return FieldAngularSpectrumVectorial{T,D,X}(deepcopy(field.nsx_X), deepcopy(field.nsy_Y), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

function copy_differentD(field::FieldAngularSpectrumVectorial{T,D,X,Y}) where {T,D,X,Y}
	return FieldAngularSpectrumVectorial{T,-D,X,Y}(deepcopy(field.nsx_X), deepcopy(field.nsy_Y), deepcopy(field.e_SXY), deepcopy(field.λ), deepcopy(field.n), deepcopy(field.ref))
end

CartesianIndices(field::FieldAngularSpectrumVectorial) = Base.CartesianIndices((2, length(field.nsx_X), length(field.nsy_Y)))
LinearIndices(field::FieldAngularSpectrumVectorial) = Base.LinearIndices((2, length(field.nsx_X), length(field.nsy_Y)))

function propagationmatrix(fieldl::FieldAngularSpectrumVectorial{T}, ref::ReferenceFrame) where T
	checkorientation(fieldl.ref, ref) || error("Cannot calculate propagation matrix as the referenceframe are not oriented")
	refΔx, refΔy, refΔz = ref.x - fieldl.ref.x, ref.y - fieldl.ref.y, ref.z - fieldl.ref.z;

	# Can't use based on referenceframe as the rotation is inverted
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, fieldl.ref.θ, fieldl.ref.ϕ);
	imk = im * 2π / fieldl.λ;
	propMatrix = Diagonal(Vector{Complex{T}}(undef, length(fieldl.e_SXY)))
	cart = CartesianIndices(fieldl)
	lin = LinearIndices(fieldl)
	@inbounds @simd for i in lin[1,:,:]
		nsz_a = dir(fieldl) * nsz(fieldl.n, fieldl.nsx_X[cart[i][2]], fieldl.nsy_Y[cart[i][3]])
		propMatrix.diag[i] = exp(imk * (fieldl.nsx_X[cart[i][2]] * refΔx + fieldl.nsy_Y[cart[i][3]] * refΔy + nsz_a * refΔz))
		propMatrix.diag[i+1] = propMatrix[i]
	end
	return propMatrix
end

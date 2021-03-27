struct Sphere{T,A<:Union{T, Complex{T}, Function, AbstractExtrapolation{T}, AbstractInterpolation{T}}} <: AbstractOpticalComponent{T}
	n::JolabFunction{T,A}
	radius::T
	ref::ReferenceFrame{T}
	function Sphere{T}(n::JolabFunction{T1,A}, r, ref) where {T,T1,A}
		return new{T,A}(n, r, ref)
	end
end

Sphere(::Type{T}, n, radius, ref) where T = Sphere{T}(convert(JolabFunction{T},n), radius, ref)
Sphere(n, radius, ref) = Sphere(Float64, n, radius, ref)

n(sphere::Sphere, λ) = sphere.n(λ)

an(n, m, x) = (m^2 * sphericalbesselj(n, m * x) * d_besselj(n, x) - sphericalbesselj(n,x) * d_besselj(n, m * x)) /
    (m^2 * sphericalbesselj(n, m * x) * d_hankelh1(n, x) - sphericalhankelh1(n,x) * d_besselj(n, m * x))
bn(n, m, x) = (sphericalbesselj(n, m * x) * d_besselj(n, x) - sphericalbesselj(n, x) * d_besselj(n, m * x)) /
    (sphericalbesselj(n, m * x) * d_hankelh1(n, x) - sphericalhankelh1(n, x) * d_besselj(n, m * x))

d_besselj(n,x) = x * sphericalbesselj(n - 1, x) - n * sphericalbesselj(n, x)
d_hankelh1(n,x) = x * sphericalhankelh1(n - 1, x) - n * sphericalhankelh1(n, x)
sphericalhankelh1(n,x) = sphericalbesselj(n,x) + im * sphericalbessely(n,x)

checkapplicability(sphere::Sphere, field::FieldAngularSpectrumScalar) = true

function get_scatteringmatrixtype(sphere::Sphere, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,D,X,B}
	sizeI = length(fieldi.e_SXY)
	r12 = Matrix{Complex{T}}(undef, sizeI, sizeI)
	t12 = Matrix{Complex{T}}(undef, sizeI, sizeI)
	r21 = Matrix{Complex{T}}(undef, sizeI, sizeI)
	t21 = Matrix{Complex{T}}(undef, sizeI, sizeI)

	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrumScalar{T,1,X,B}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, sphere.ref)
	else
		fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, sphere.ref)
		fieldr = FieldAngularSpectrumScalar{T,1,X,B}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, fieldi.ref)
	end
	return ScatteringMatrix{T, FieldAngularSpectrumScalar{T,-1,X,B}, FieldAngularSpectrumScalar{T,1,X,B}, Matrix{Complex{T}}, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
end

@inline function S1S2(cosθ::Union{T,Complex{T}}, an_N, bn_N) where T
	Πn2, Πn1, Πn = zero(Complex{T}), one(Complex{T}), one(Complex{T})
	τn = cosθ
	s1 = (2 + 1) / (1 * (1 + 1)) * (an_N[1] * Πn + bn_N[1] * τn)
	s2 = (2 + 1) / (1 * (1 + 1)) * (an_N[1] * τn + bn_N[1] * Πn)
	@inbounds for n in 2:length(an_N)
    	Πn = (2n - 1) / (n - 1) * cosθ * Πn1 - n / (n-1) * Πn2
    	τn = n * cosθ * Πn - (n + 1) * Πn1
		s1 += (2n + 1) / (n * (n + 1)) * (an_N[n] * Πn + bn_N[n] * τn)
		s2 += (2n + 1) / (n * (n + 1)) * (an_N[n] * τn + bn_N[n] * Πn)
    	Πn2 = Πn1
    	Πn1 = Πn
	end
	return (s1, s2)
end

function coefficient_general(sphere::Sphere{T}, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,D,X,B}
	checkapplicability(sphere, fieldi)

	sizeI = length(fieldi.e_SXY)
	m = n(sphere, fieldi.λ) / fieldi.n
	x = sphere.radius / fieldi.λ

	k = 2π / fieldi.λ
	scat = get_scatteringmatrixtype(sphere, fieldi)

    n_max = round(Int, (2 + x + 4 * x^(1/3)))
	an_N = an.(1:n_max, m, x)
	bn_N = bn.(1:n_max, m, x)
	scat.t₁₂ .= 0
	scat.r₁₂ .= 0
	cons = 1 / (2π * k^2 * fieldi.n^2) # To correct for the MIE theory which gives the far field
	cart_i = CartesianIndices(fieldi)
	@inbounds Threads.@threads for i_i in iterator_index(fieldi)
		iX1, iY1 = cart_i[i_i][2], cart_i[i_i][3]
		nsx_i = fieldi.nsx_X[iX1]
		nsy_i = fieldi.nsy_Y[iY1]
		dkxdky = k^2 * Δvector(fieldi.nsx_X, iX1) * Δvector(fieldi.nsy_Y, iY1)
		nsz_i = nsz(fieldi.n, nsx_i, nsy_i)
		if imag(nsz_i) > @tol
			view(scat.t₁₂,:,i_i) .= zero(Complex{T})
			view(scat.r₁₂,:,i_i) .= zero(Complex{T})
			continue
		end
		real_nsz_i = real(nsz_i)
		for i_s in iterator_index(fieldi)
			iX2, iY2 = cart_i[i_s][2], cart_i[i_s][3]
			nsx_s = fieldi.nsx_X[iX2]
			nsy_s = fieldi.nsy_Y[iY2]
			nsz_s = -nsz(fieldi.n, nsx_s, nsy_s)
			if imag(nsz_s) < @tol
				real_nsz_s = real(nsz_s)
				cosθr = (nsx_s * nsx_i + nsy_i * nsy_s + real_nsz_s * real_nsz_i) / real(fieldi.n)^2
				(scat.r₁₂[i_s, i_i], ~) = S1S2(cosθr, an_N, bn_N)
				scat.r₁₂[i_s, i_i] *= dkxdky * cons / -real_nsz_s

				real_nsz_s = -real_nsz_s
				cosθr = (nsx_s * nsx_i + nsy_i * nsy_s + real_nsz_s * real_nsz_i) / real(fieldi.n)^2
				(scat.t₁₂[i_s, i_i], ~) = S1S2(cosθr, an_N, bn_N)
				scat.t₁₂[i_s, i_i] *= dkxdky * cons / real_nsz_s
			else
				scat.t₁₂[i_s, i_i] = zero(Complex{T})
				scat.r₁₂[i_s, i_i] = zero(Complex{T})
			end
		end
    end
	for i in 1:size(scat.t₁₂, 1)
		scat.t₁₂[i,i] += 1
	end
	scat.t₂₁ .= scat.t₁₂
	scat.r₂₁ .= scat.r₁₂

	correctscatteringmatrix_referenceframes!(scat, sphere, fieldi)
	return scat
end

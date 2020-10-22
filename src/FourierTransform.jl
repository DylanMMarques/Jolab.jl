struct FourierTransform{T, X <:Union{AbstractVector{T}, Nothing}, Y<:Union{ReferenceFrame{T}, Nothing}} <: AbstractOpticalComponent{T}
	x_X::X
	y_Y::X
	nsx_X::X
	nsy_Y::X
	ref::Y
end

function FourierTransform(x_X::X1, y_Y::X2, nsx_X::X3, nsy_Y::X4) where {X1 <: AbstractVector{T}, X2 <: AbstractVector{T}, X3 <: AbstractVector{T}, X4 <: AbstractVector{T}} where T <: Real
	X = promote_type(X1, X2, X3, X4)
	return FourierTransform{T, X, Nothing}(x_X, y_Y, nsx_X, nsy_Y, nothing)
end

FourierTransform{T}() where T = FourierTransform{T, Nothing}(nothing, nothing, nothing, nothing)
FourierTransform() = FourierTransform{Float64}()

@inline checkapplicability(fourier::FourierTransform, field::Union{FieldAngularSpectrum{T}, FieldSpace{T}}) where {T} = true

@inline function t(fourier::FourierTransform{T,X,Y}, field::FieldSpace, iX, iY, iA, iB) where {T, X, Y<:Nothing}
	(xmin, xmax) = integralExtremes(field.x_X, iX)
	(ymin, ymax) = integralExtremes(field.y_Y, iY)

	k = 2π / field.λ

	return integrate_exp_x_y(-k * fourier.nsx_X[iA], -k * fourier.nsy_Y[iB], zero(T), xmin, xmax, ymin, ymax) / 4π^2
end

@inline function tinv(fourier::FourierTransform{T,X,Y}, field::FieldSpace, iX, iY, iA, iB) where {T, X, Y<:Nothing}
	(nsxmin, nsxmax) = integralExtremes(fourier.nsx_X, iA)
	(nsymin, nsymax) = integralExtremes(fourier.nsy_Y, iB)

	k = 2π / field.λ

	return integrate_exp_x_y(k * field.x_X[iX], k * field.y_Y[iY], zero(T), nsxmin, nsxmax, nsymin, nsymax) * k^2
end

@inline function t(fourier::FourierTransform{T,X,Y}, field::FieldAngularSpectrum, iX, iY, iA, iB) where {T, X, Y<:Nothing}
	(nsxmin, nsxmax) = integralExtremes(field.nsx_X, iX)
	(nsymin, nsymax) = integralExtremes(field.nsy_Y, iY)

	k = 2π / field.λ

	return k^2 * integrate_exp_x_y(k * fourier.x_X[iA], k*fourier.y_Y[iB], zero(T), nsxmin, nsxmax, nsymin, nsymax)
end

@inline function tinv(fourier::FourierTransform{T,X,Y}, field::FieldAngularSpectrum, iX, iY, iA, iB) where {T, X, Y<:Nothing}
	(xmin, xmax) = integralExtremes(fourier.x_X, iA)
	(ymin, ymax) = integralExtremes(fourier.y_Y, iB)

	k = 2π / field.λ

	return integrate_exp_x_y(-k * field.nsx_X[iX], -k * field.nsy_Y[iY], zero(T), xmin, xmax, ymin, ymax) / 4π^2
end

@inline function t(fourier::FourierTransform{T,X,Y}, field::FieldAngularSpectrum{T}, iX, iY, iA, iB) where {T, X, Y<:ReferenceFrame}
	(nsxmin, nsxmax) = integralExtremes(field.nsx_X, iX)
	(nsymin, nsymax) = integralExtremes(field.nsy_Y, iY)

	k = 2π / field.λ

	refΔx, refΔy, refΔz = fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z
	(refΔx2, refΔy2, refΔz2) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ)

	imagina_waves = (imag(field.n) > @tol) || nsxmin^2 + nsymin^2 > real(field.n)^2 || nsxmin^2 + nsymax^2 > real(field.n)^2 || nsxmax^2 + nsymin^2 > real(field.n)^2 || nsxmax^2 + nsymax^2 > real(field.n)^2

	if imagina_waves
		f(nsr) = k * ((fourier.x_X[iA] + refΔx2) * nsr[1] + (fourier.y_Y[iB] + refΔy2) * nsr[2] + dir(field) * √(complex(field.n^2 - nsr[1]^2 - nsr[2]^2)) * refΔz2)
		t12 = k^2 * hcubature(f, SVector(nsxmin, nsymin), SVector(nsxmax, nsymax))[1]
	else
		@inline g(nsx, nsy) = k * ((fourier.x_X[iA] + refΔx2) * nsx + (fourier.y_Y[iB] + refΔy2) * nsy + dir(field) * √(real(field.n)^2 - nsx^2 - nsy^2) * refΔz2)
		t12 = k^2 * integrate_exp_xy_x_y(g, nsxmin, nsxmax, nsymin, nsymax)
	end
	return t12
end

@inline function tinv(fourier::FourierTransform{T,X,Y}, field::FieldAngularSpectrum{T}, iX, iY, iA, iB) where {T, X, Y<:ReferenceFrame}
	(xmin, xmax) = integralExtremes(fourier.x_X, iA)
	(ymin, ymax) = integralExtremes(fourier.y_Y, iB)

	k = 2π / field.λ

	refΔx, refΔy, refΔz = fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z
	(refΔx2, refΔy2, refΔz2) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ)

	aux = exp(-im * k * refΔz2 * dir(field) * √(field.n^2 - fourier.nsx_X[iX]^2 - fourier.nsy_Y[iY]^2)) # miss ref.x factors

	return aux * integrate_exp_x_y(-k * field.nsx_X[iX], -k * field.nsy_Y[iY], zero(T), xmin, xmax, ymin, ymax) / 4π^2
end

@inline function t(fourier::FourierTransform{T,X,Y}, field::FieldSpace, iX, iY, iA, iB) where {T, X, Y<:ReferenceFrame}
	(xmin, xmax) = integralExtremes(field.x_X, iX)
	(ymin, ymax) = integralExtremes(field.y_Y, iY)

	k = 2π / field.λ

	refΔx, refΔy, refΔz = fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z
	(refΔx2, refΔy2, refΔz2) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ)

	aux = exp(im * k * refΔz2 * dir(field) * √(field.n^2 - fourier.nsx_X[iA]^2 - fourier.nsy_Y[iB]^2)) #MISS ref factor

	return aux * integrate_exp_x_y(-k * fourier.nsx_X[iA], -k * fourier.nsy_Y[iB], zero(T), xmin, xmax, ymin, ymax) / 4π^2
end

@inline function tinv(fourier::FourierTransform{T,X,Y}, field::FieldSpace, iX, iY, iA, iB) where {T, X, Y<:ReferenceFrame}
	(nsxmin, nsxmax) = integralExtremes(fourier.nsx_X, iA)
	(nsymin, nsymax) = integralExtremes(fourier.nsy_Y, iB)

	k = 2π / field.λ

	refΔx, refΔy, refΔz = fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z
	(refΔx2, refΔy2, refΔz2) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ)

	imagina_waves = (imag(field.n) > @tol) || nsxmin^2 + nsymin^2 > real(field.n)^2 || nsxmin^2 + nsymax^2 > real(field.n)^2 || nsxmax^2 + nsymin^2 > real(field.n)^2 || nsxmax^2 + nsymax^2 > real(field.n)^2
	if imagina_waves
		@inline f(nsr) = k * ((fourier.x_X[iX] + refΔx2) * nsr[1] + (fourier.y_Y[iY] + refΔy2) * nsr[2] - dir(field) * √(complex(field.n^2 - nsr[1]^2 - nsr[2]^2)) * refΔz2)
		t = k^2 * hcubature(f, SVector(nsxmin, nsymin), SVector(nsxmax, nsymax))
	else
		@inline g(nsx, nsy) = k * ((fourier.x_X[iX] + refΔx2) * nsx + (fourier.y_Y[iY] + refΔy2) * nsy + dir(field) * √(real(field.n)^2 - nsx^2 - nsy^2) * refΔz2)
		t = k^2 * integrate_exp_xy_x_y(g, nsxmin, nsxmax, nsymin, nsymax)
	end
	return t
end

@inline function et(fourier::FourierTransform{T,X,Y}, field::FieldAngularSpectrum, angleE, iX,iY,iA,iB) where {T, X, Y<:ReferenceFrame}
	nsxmin, nsxmax = field.nsx_X[iX], field.nsx_X[iX+1]
	nsymin, nsymax = field.nsy_Y[iY], field.nsy_Y[iY+1]

	k = 2π / field.λ

	refΔx, refΔy, refΔz = fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z
	(refΔx2, refΔy2, refΔz2) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ)

	imagina_waves = (imag(field.n) > @tol) || nsxmin^2 + nsymin^2 > real(field.n)^2 || nsxmin^2 + nsymax^2 > real(field.n)^2 || nsxmax^2 + nsymin^2 > real(field.n)^2 || nsxmax^2 + nsymax^2 > real(field.n)^2

	# if imagina_waves
	# 	f(nsr) = k * ((fourier.x_X[iA] + refΔx2) * nsr[1] + (fourier.y_Y[iB] + refΔy2) * nsr[2] + dir(field) * √(complex(field.n^2 - nsr[1]^2 - nsr[2]^2)) * refΔz2)
	# 	et = space.e_SXY[1,iX,iY] * k^2 * hcubature(f, SVector(nsxmin, nsymin), SVector(nsxmax, nsymax))[1]
	# 	@show "this is wrong"
	# else
		@inline g(nsx, nsy) = k * ((fourier.x_X[iA] + refΔx2) * nsx + (fourier.y_Y[iB] + refΔy2) * nsy + dir(field) * √(real(field.n)^2 - nsx^2 - nsy^2) * refΔz2)
		(α, β, γ, δ) = bilinearinterpolation(angleE[iX,iY] + g(nsxmin,nsymin), angleE[iX,iY+1] + g(nsxmin, nsymax), angleE[iX+1,iY] + g(nsxmax, nsymin), angleE[iX+1,iY+1] + g(nsxmax,nsymax), nsxmin, nsxmax, nsymin, nsymax)
		(a,b,c,d) = bilinearinterpolation(abs(field.e_SXY[1,iX,iY]), abs(field.e_SXY[1,iX,iY+1]), abs(field.e_SXY[1,iX+1,iY]), abs(field.e_SXY[1,iX+1,iY+1]), nsxmin, nsxmax, nsymin, nsymax)

		et = k^2 * integrate_xy_x_y_d_exp_xy_x_y(a, b, c, d, α, β, γ, δ, nsxmin, nsxmax, nsymin, nsymax)
	# end
	return et

end


getsizes(fourier::FourierTransform{T,X}, field::FieldSpace) where {T, X<:AbstractVector} = (length(fourier.nsx_X), length(fourier.nsy_Y))
getsizes(fourier::FourierTransform{T,X}, field::FieldAngularSpectrum) where {T, X<:AbstractVector} = (length(fourier.x_X), length(fourier.y_Y))

function get_scatteringmatrixtype(fourier::FourierTransform{T,Y,R}, fieldi::FieldSpace{T,D,X}) where {T,D,Y, X<:AbstractVector, R<:ReferenceFrame}
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(fourier, fieldi)
	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)
	m = Array{Complex{T},3}(undef,1,sizeA, sizeB)
	if dir(fieldi) > 0
		fieldl = FieldSpace{T,-1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrum{T,1,Y}(deepcopy(fourier.nsx_X), deepcopy(fourier.nsy_Y), m, fieldi.λ, fieldi.n, fourier.ref)
		return ScatteringMatrix{T,FieldSpace{T,-1,X}, FieldAngularSpectrum{T,1,Y},Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldAngularSpectrum{T,-1,Y}(deepcopy(fourier.nsx_X), deepcopy(fourier.nsy_Y), m, fieldi.λ, fieldi.n, fourier.ref)
		fieldr = FieldSpace{T,1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldAngularSpectrum{T,-1,Y}, FieldSpace{T,1,X}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	end
end

function get_scatteringmatrixtype(fourier::FourierTransform{T,Y,R}, fieldi::FieldSpace{T,D,X}) where {T, D,X<:AbstractVector, R<:Nothing, Y}
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(fourier, fieldi)
	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)
	m = Array{Complex{T},3}(undef,1,sizeA, sizeB)
	if dir(fieldi) > 0
		fieldl = FieldSpace{T,-1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrum{T,1,Y}(deepcopy(fourier.nsx_X), deepcopy(fourier.nsy_Y), m, fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T,FieldSpace{T,-1,X}, FieldAngularSpectrum{T,1,Y},Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldAngularSpectrum{T,-1,Y}(deepcopy(fourier.nsx_X), deepcopy(fourier.nsy_Y), m, fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpace{T,1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldAngularSpectrum{T,-1,Y}, FieldSpace{T,1,X}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	end
end

function get_scatteringmatrixtype(fourier::FourierTransform{T,Y,R}, fieldi::FieldAngularSpectrum{T,D,X}) where {T,D,X,Y,R<:Nothing}
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(fourier, fieldi)

	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)

	m = Array{Complex{T},3}(undef,1,length(fourier.x_X), length(fourier.y_Y))
	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrum{T,-1,Y}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpace{T,1,X}(deepcopy(fourier.x_X), deepcopy(fourier.y_Y), m, fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldAngularSpectrum{T,-1,Y}, FieldSpace{T,1,X}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldSpace{T,-1,X}(deepcopy(fourier.x_X), deepcopy(fourier.y_Y), m, fieldi.λ, fieldi.n, field.ref)
		fieldr = FieldAngularSpectrum{T,1,Y}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldSpace{T,-1,X}, FieldAngularSpectrum{T,1,Y}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	end
end

function get_scatteringmatrixtype(fourier::FourierTransform{T,Y,R}, fieldi::FieldAngularSpectrum{T,D,X}) where {T,D,X,Y,R<:ReferenceFrame}
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(fourier, fieldi)

	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeA * sizeB, sizeX * sizeY)
	t21 = Matrix{Complex{T}}(undef, sizeX * sizeY, sizeA * sizeB)

	m = Array{Complex{T},3}(undef,1,length(fourier.x_X), length(fourier.y_Y))
	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrum{T,-1,Y}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpace{T,1,X}(deepcopy(fourier.x_X), deepcopy(fourier.y_Y), m, fieldi.λ, fieldi.n, fourier.ref)
		return ScatteringMatrix{T, FieldAngularSpectrum{T,-1,Y}, FieldSpace{T,1,X}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldSpace{T,-1,X}(deepcopy(fourier.x_X), deepcopy(fourier.y_Y), m, fieldi.λ, fieldi.n, fourier.ref)
		fieldr = FieldAngularSpectrum{T,1,Y}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldSpace{T,-1,X}, FieldAngularSpectrum{T,1,Y}, Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	end
end

function getfields_lr(fourier::FourierTransform{T,Y,R}, fieldi::FieldSpace{T,D,X}) where {T,D,Y, X<:AbstractVector, R<:ReferenceFrame}
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(fourier, field)

	m = Array{Complex{T},3}(undef,1,sizeA, sizeB)
	if dir(fieldi) > 0
		fieldl = FieldSpace{T,-1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrum{T,1,Y}(deepcopy(fourier.nsx_X), deepcopy(fourier.nsy_Y), m, fieldi.λ, fieldi.n, fourier.ref)
	else
		fieldl = FieldAngularSpectrum{T,-1,Y}(deepcopy(fourier.nsx_X), deepcopy(fourier.nsy_Y), m, fieldi.λ, fieldi.n, fourier.ref)
		fieldr = FieldSpace{T,1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
	end
	return (fieldl, fieldr)
end

function getfields_lr(fourier::FourierTransform{T,Y,R}, fieldi::FieldSpace{T,D,X}) where {T,Y,D,X<:AbstractVector, R<:Nothing}
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(fourier, fieldi)
	m = Array{Complex{T},3}(undef,1,sizeA, sizeB)
	if dir(fieldi) > 0
		fieldl = FieldSpace{T,-1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrum{T,1,Y}(deepcopy(fourier.nsx_X), deepcopy(fourier.nsy_Y), m, fieldi.λ, fieldi.n, fieldi.ref)
	else
		fieldl = FieldAngularSpectrum{T,-1,Y}(deepcopy(fourier.nsx_X), deepcopy(fourier.nsy_Y), m, fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpace{T,1,X}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
	end
	return (fieldl, fieldr)
end

function getfields_lr(fourier::FourierTransform{T,Y,R}, fieldi::FieldAngularSpectrum{T,D,X}) where {T,D,X,Y,R<:Nothing}
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(fourier, fieldi)

	m = Array{Complex{T},3}(undef,1,length(fourier.x_X), length(fourier.y_Y))
	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrum{T,-1,Y}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpace{T,1,X}(deepcopy(fourier.x_X), deepcopy(fourier.y_Y), m, fieldi.λ, fieldi.n, fieldi.ref)
	else
		fieldl = FieldSpace{T,-1,X}(deepcopy(fourier.x_X), deepcopy(fourier.y_Y), m, fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrum{T,1,Y}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
	end
	return (fieldl, fieldr)
end

function getfields_lr(fourier::FourierTransform{T,Y,R}, fieldi::FieldAngularSpectrum{T,D,X}) where {T,D,X,Y,R<:ReferenceFrame}
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(fourier, fieldi)

	m = Array{Complex{T},3}(undef,1,length(fourier.x_X), length(fourier.y_Y))
	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrum{T,-1,Y}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpace{T,1,X}(deepcopy(fourier.x_X), deepcopy(fourier.y_Y), m, fieldi.λ, fieldi.n, fourier.ref)
	else
		fieldl = FieldSpace{T,-1,X}(deepcopy(fourier.x_X), deepcopy(fourier.y_Y), m, fieldi.λ, fieldi.n, fourier.ref)
		fieldr = FieldAngularSpectrum{T,1,Y}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
	end
	return (fieldl, fieldr)
end

function coefficient_general(fourier::FourierTransform{T,X}, field::AbstractFieldMonochromatic) where {T, X <: AbstractVector}
	checkapplicability(fourier, field) || tobedone()

	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(fourier, field)

	scat = get_scatteringmatrixtype(fourier, field)

	coordAB = LinearIndices((sizeA, sizeB))
	coordXY = LinearIndices((sizeX, sizeY))

	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let scat=scat, fourier = fourier, field = field
		Threads.@threads for iX1 in 1:sizeX
	 	 	for iY1 in 1:sizeY
				i1 = coordXY[iX1, iY1]
				for iA2 in 1:sizeA
					for iB2 in 1:sizeB
						i2 = coordAB[iA2, iB2]
						if dir(field) > 0
							scat.t₁₂[i2, i1] = t(fourier, field, iX1, iY1, iA2, iB2)
							scat.t₂₁[i1, i2] = tinv(fourier, field, iX1, iY1, iA2, iB2)
						else
							scat.t₂₁[i2, i1] = t(fourier, field, iX1, iY1, iA2, iB2)
							scat.t₁₂[i1, i2] = tinv(fourier, field, iX1, iY1, iA2, iB2)
						end
					end
				end
			end
		end
	end
	return scat
end

function lightinteraction(fourier::FourierTransform{T,X}, field::AbstractFieldMonochromatic) where {T, X <: AbstractVector}
	checkapplicability(fourier, field) || tobedone()

	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	(sizeA, sizeB) = getsizes(fourier, field)

	(fieldl, fieldr) = getfields_lr(fourier, field)

	coordAB = LinearIndices((sizeA, sizeB))
	coordXY = LinearIndices((sizeX, sizeY))
	fieldl.e_SXY .= zero(Complex{T})
	fieldr.e_SXY .= zero(Complex{T})

	# angleE = unwrap!(angle.(view(field.e_SXY,1,:,:)), dims = 1:2)

	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let fieldl = fieldl, fieldr = fieldr, fourier = fourier, field = field
		 @inbounds Threads.@threads for iA2 in 1:sizeA
			for iB2 in 1:sizeB
				i2 = coordAB[iA2, iB2]
				for iX1 in 1:sizeX
	 	 			for iY1 in 1:sizeY
						i1 = coordXY[iX1, iY1]
						if dir(field) > 0
							fieldr.e_SXY[i2] += t(fourier, field, iX1, iY1, iA2, iB2) * field.e_SXY[i1]
						else
							fieldl.e_SXY[i2] += t(fourier, field, iX1, iY1, iA2, iB2) * field.e_SXY[i1]
						end
					end
				end
			end
		end
	end
	return (fieldl, fieldr)
end

@inline function coefficient_general(fourier::FourierTransform{T,X}, field::FieldSpace{T,D,Y}) where {T, X<:Nothing, Y <: AbstractRange, D}
	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	Δx = Δvector(field.x_X, 1)
	Δy = Δvector(field.y_Y, 1)
	nsx_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δx)) * field.λ;
	nsy_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δy)) * field.λ;
	fourieraux = FourierTransform{T,Y}(field.x_X, field.y_Y, nsx_X, nsy_Y)
	return coefficient_general(fourieraux, field)
end

@inline function coefficient_general(fourier::FourierTransform{T,X}, field::FieldAngularSpectrum{T,D,Y}) where {T,D, X <: Nothing, Y <: AbstractRange}
	(sizeX, sizeY) = size(field.e_SXY)[2:3]
	Δnsx = Δvector(field.nsx_X, 1)
	Δnsy = Δvector(field.nsy_Y, 1)
	x_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δnsx)) * field.λ;
	y_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δnsy)) * field.λ;
	fourieraux = FourierTransform{T,Y}(x_X, y_Y, field.nsx_X, field.nsy_Y)
	return coefficient_general(fourieraux, field)
end

function coefficient_specific(fourier::FourierTransform{T,X}, field::F) where {T, X <: Nothing, F <: FieldSpace{T, <:AbstractRange}}
	planfft = FFTW.plan_fft(view(field.e_SXY,1,:,:))
	(sizeX, sizeY) = (length(field.x_X), length(field.y_Y))

	(Δx, Δy) = (Δvector(field.x_X, 1), Δvector(field.y_Y, 1))

	nsx_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δx)) * field.λ;
	nsy_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δy)) * field.λ;
	scalar = Δx * Δy / 4 / π^2

	(Δnsx, Δnsy) = (Δvector(nsx_X, 1), Δvector(nsy_Y, 1))
	i_scalar = sizeX * sizeY * Δnsx * Δnsy * 4 * π^2 / field.λ^2

	tmp = Array{Complex{T},2}(undef, sizeX, sizeY)
	if dir(field) > 0
		fieldr = FieldAngularSpectrum{T, 1}(nsx_X, nsy_Y, field.e_SXY, field.λ, field.n, field.ref)
		return FFTCoefficient{T,F,typeof(fieldr),typeof(planfft)}(planfft, scalar, i_scalar, tmp, field, fieldr)
	else
		fieldl = FieldAngularSpectrum{T,-1}(nsx_X, nsy_Y, field.e_SXY, field.λ, field.n, field.ref)
		return FFTCoefficient{T,typeof(fieldl),F,typeof(planfft)}(planfft, scalar, i_scalar, tmp, fieldl, field)
	end
end

function coefficient_specific(fourier::FourierTransform{T,X}, field::F) where {T, X <: Nothing, F <: FieldAngularSpectrum{T, <:AbstractRange}}
	planfft = FFTW.plan_fft(view(field.e_SXY,1,:,:))
	(sizeX, sizeY) = (length(field.nsx_X), length(field.nsy_Y))

	(Δnsx, Δnsy) = (Δvector(field.nsx_X, 1), Δvector(field.nsy_Y, 1))

	x_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δnsx)) * field.λ;
	y_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δnsy)) * field.λ;
	(Δx, Δy) = (Δvector(x_X, 1), Δvector(y_Y, 1))
	scalar = Δx * Δy / 4 / π^2

	i_scalar = sizeX * sizeY * Δnsx * Δnsy * 4 * π^2 / field.λ^2

	tmp = Array{Complex{T},2}(undef, sizeX, sizeY)
	if dir(field) > 0
		fieldr = FieldSpace{T,1}(x_X, y_Y, field.e_SXY, field.λ, field.n, field.ref)
		return FFTCoefficient{T,F,typeof(fieldr),typeof(planfft)}(planfft, scalar, i_scalar, tmp, field, fieldr)
	else
		fieldl = FieldSpace{T,-1}(x_X, y_Y, field.e_SXY, field.λ, field.n, field.ref)
		return FFTCoefficient{T,typeof(fieldl),F,typeof(planfft)}(planfft, scalar, i_scalar, tmp, fieldl, field)
	end
end

struct FFTCoefficient{T,L,R, F <: AbstractFFTs.Plan} <: AbstractCoefficient{T,L,R}
	planfft::F
	scalar::T
	i_scalar::T
	tmp::Array{Complex{T}, 2}
	fieldl::L
	fieldr::R
end

function lightinteraction!(fieldl::L, fieldr::R, coef::FFTCoefficient{T,L,R}, fieldi::Union{L,R}) where {T,L,R}
	if dir(fieldi) > 0
		samedefinitions(fieldl, fieldi) || tobedone()
		if fieldi isa FieldSpace
			mul!(coef.tmp, coef.planfft, fftshift(view(fieldi.e_SXY,1,:,:)))
			vec(fieldr.e_SXY) .= vec(fftshift(coef.tmp)) * coef.scalar
		else
			ldiv!(coef.tmp, coef.planfft, fftshift(view(fieldi.e_SXY,1,:,:)))
			vec(fieldr.e_SXY) .= vec(fftshift(coef.tmp)) * coef.i_scalar
		end
		fieldl.e_SXY .= zero(Complex{T})
	else
		samedefinitions(fieldr, fieldi) || tobedone()
		if fieldi isa FieldSpace
			mul!(coef.tmp, coef.planfft, fftshift(view(fieldi.e_SXY,1,:,:)))
			vec(fieldl.e_SXY) .= vec(fftshift(coef.tmp)) * coef.scalar
		else
			ldiv!(coef.tmp, coef.planfft, fftshift(view(fieldi.e_SXY,1,:,:)))
			vec(fieldl.e_SXY) .= vec(fftshift(coef.tmp)) * coef.i_scalar
		end
		fieldr.e_SXY .= zero(Complex{T})
	end
end

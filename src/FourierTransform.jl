struct FourierTransform{T, X <:Union{AbstractVector{T}, Nothing}, Y<:Union{ReferenceFrame{T}, Nothing}} <: AbstractOpticalComponent{T}
	x_X::X
	y_Y::X
	nsx_X::X
	nsy_Y::X
	ref::Y
end

function FourierTransform(x_X::X1, y_Y::X2, nsx_X::X3, nsy_Y::X4, ref::R) where {X1 <: AbstractVector{T}, X2 <: AbstractVector{T}, X3 <: AbstractVector{T}, X4 <: AbstractVector{T}, R <: ReferenceFrame{T}} where T <: Real
	X = promote_type(X1, X2, X3, X4)
	return FourierTransform{T, X ,R}(x_X, y_Y, nsx_X, nsy_Y, ref)
end
function FourierTransform(x_X::X1, y_Y::X2, nsx_X::X3, nsy_Y::X4) where {X1 <: AbstractVector{T}, X2 <: AbstractVector{T}, X3 <: AbstractVector{T}, X4 <: AbstractVector{T}} where T <: Real
	X = promote_type(X1, X2, X3, X4)
	return FourierTransform{T, X, Nothing}(x_X, y_Y, nsx_X, nsy_Y, nothing)
end

FourierTransform{T}() where T = FourierTransform{T, Nothing}(nothing, nothing, nothing, nothing)
FourierTransform() = FourierTransform{Float64}()

checkapplicability(fourier::FourierTransform, field::Union{AbstractFieldAngularSpectrum, AbstractFieldSpace}) = true

function t(fourier::FourierTransform{T,X,Y}, field::FieldSpaceScalar, iX, iY, iA, iB) where {T, X, Y<:Nothing}
	(xmin, xmax) = integralExtremes(field.x_X, iX)
	(ymin, ymax) = integralExtremes(field.y_Y, iY)

	k = 2π / field.λ

	return integrate_exp_x_y(-k * fourier.nsx_X[iA], -k * fourier.nsy_Y[iB], zero(T), xmin, xmax, ymin, ymax) / 4π^2
end

function tinv(fourier::FourierTransform{T,X,Y}, field::AbstractFieldSpace, iX, iY, iA, iB) where {T, X, Y<:Nothing}
	(nsxmin, nsxmax) = integralExtremes(fourier.nsx_X, iA)
	(nsymin, nsymax) = integralExtremes(fourier.nsy_Y, iB)

	k = 2π / field.λ

	return integrate_exp_x_y(k * field.x_X[iX], k * field.y_Y[iY], zero(T), nsxmin, nsxmax, nsymin, nsymax) * k^2
end

function t(fourier::FourierTransform{T,X,Y}, field::AbstractFieldAngularSpectrum, iX, iY, iA, iB) where {T, X, Y<:Nothing}
	(nsxmin, nsxmax) = integralExtremes(field.nsx_X, iX)
	(nsymin, nsymax) = integralExtremes(field.nsy_Y, iY)

	k = 2π / field.λ

	return k^2 * integrate_exp_x_y(k * fourier.x_X[iA], k*fourier.y_Y[iB], zero(T), nsxmin, nsxmax, nsymin, nsymax)
end

function tinv(fourier::FourierTransform{T,X,Y}, field::AbstractFieldAngularSpectrum, iX, iY, iA, iB) where {T, X, Y<:Nothing}
	(xmin, xmax) = integralExtremes(fourier.x_X, iA)
	(ymin, ymax) = integralExtremes(fourier.y_Y, iB)

	k = 2π / field.λ

	return integrate_exp_x_y(-k * field.nsx_X[iX], -k * field.nsy_Y[iY], zero(T), xmin, xmax, ymin, ymax) / 4π^2
end

function t(fourier::FourierTransform{T,X,Y}, field::AbstractFieldAngularSpectrum, iX, iY, iA, iB) where {T, X, Y<:ReferenceFrame}
	(nsxmin, nsxmax) = integralExtremes(field.nsx_X, iX)
	(nsymin, nsymax) = integralExtremes(field.nsy_Y, iY)

	k = 2π / field.λ

	refΔx, refΔy, refΔz = fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z
	(refΔx2, refΔy2, refΔz2) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ)

	imagina_waves = (imag(field.n) > @tol) || nsxmin^2 + nsymin^2 > real(field.n)^2 || nsxmin^2 + nsymax^2 > real(field.n)^2 || nsxmax^2 + nsymin^2 > real(field.n)^2 || nsxmax^2 + nsymax^2 > real(field.n)^2

	if imagina_waves
		f(nsx, nsy) = k * ((fourier.x_X[iA] + refΔx2) * nsx + (fourier.y_Y[iB] + refΔy2) * nsy)
		t12 = exp(im * k * dir(field) * √(complex(field.n^2 - field.nsx_X[iX]^2 - field.nsy_Y[iY]^2)) * refΔz2) * k^2 *  integrate_exp_xy_x_y(f, nsxmin, nsxmax, nsymin, nsymax)
	else
		g(nsx, nsy) = k * ((fourier.x_X[iA] + refΔx2) * nsx + (fourier.y_Y[iB] + refΔy2) * nsy + dir(field) * √(real(field.n)^2 - nsx^2 - nsy^2) * refΔz2)
		t12 = k^2 * integrate_exp_xy_x_y(g, nsxmin, nsxmax, nsymin, nsymax)
	end
	return t12
end

function tinv(fourier::FourierTransform{T,X,Y}, field::AbstractFieldAngularSpectrum, iX, iY, iA, iB) where {T, X, Y<:ReferenceFrame}
	(xmin, xmax) = integralExtremes(fourier.x_X, iA)
	(ymin, ymax) = integralExtremes(fourier.y_Y, iB)

	k = 2π / field.λ

	refΔx, refΔy, refΔz = fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z
	(refΔx2, refΔy2, refΔz2) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ)

	aux = exp(-im * k * refΔz2 * dir(field) * √(field.n^2 - fourier.nsx_X[iX]^2 - fourier.nsy_Y[iY]^2)) # miss ref.x factors

	return aux * integrate_exp_x_y(-k * field.nsx_X[iX], -k * field.nsy_Y[iY], zero(T), xmin, xmax, ymin, ymax) / 4π^2
end

function t(fourier::FourierTransform{T,X,Y}, field::AbstractFieldSpace, iX, iY, iA, iB) where {T, X, Y<:ReferenceFrame}
	(xmin, xmax) = integralExtremes(field.x_X, iX)
	(ymin, ymax) = integralExtremes(field.y_Y, iY)

	k = 2π / field.λ

	refΔx, refΔy, refΔz = fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z
	(refΔx2, refΔy2, refΔz2) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ)

	aux = exp(im * k * refΔz2 * dir(field) * √(field.n^2 - fourier.nsx_X[iA]^2 - fourier.nsy_Y[iB]^2)) #MISS ref factor

	return aux * integrate_exp_x_y(-k * fourier.nsx_X[iA], -k * fourier.nsy_Y[iB], zero(T), xmin, xmax, ymin, ymax) / 4π^2
end

function tinv(fourier::FourierTransform{T,X,Y}, field::AbstractFieldSpace, iX, iY, iA, iB) where {T, X, Y<:ReferenceFrame}
	(nsxmin, nsxmax) = integralExtremes(fourier.nsx_X, iA)
	(nsymin, nsymax) = integralExtremes(fourier.nsy_Y, iB)

	k = 2π / field.λ

	refΔx, refΔy, refΔz = fourier.ref.x - field.ref.x, fourier.ref.y - field.ref.y, fourier.ref.z - field.ref.z
	(refΔx2, refΔy2, refΔz2) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, field.ref.θ, field.ref.ϕ)

	imagina_waves = (imag(field.n) > @tol) || nsxmin^2 + nsymin^2 > real(field.n)^2 || nsxmin^2 + nsymax^2 > real(field.n)^2 || nsxmax^2 + nsymin^2 > real(field.n)^2 || nsxmax^2 + nsymax^2 > real(field.n)^2
	if imagina_waves
		f(nsx, nsy) = k * ((fourier.x_X[iX] + refΔx2) * nsx + (fourier.y_Y[iY] + refΔy2) * nsy)
		t = exp(-im * k * dir(field) * √(complex(field.n^2 - fourier.nsx_X[iA]^2 - fourier.nsy_Y[iB]^2)) * refΔz2) * k^2 *  integrate_exp_xy_x_y(f, nsxmin, nsxmax, nsymin, nsymax)
	else
		g(nsx, nsy) = k * ((fourier.x_X[iX] + refΔx2) * nsx + (fourier.y_Y[iY] + refΔy2) * nsy + dir(field) * √(real(field.n)^2 - nsx^2 - nsy^2) * refΔz2)
		t = k^2 * integrate_exp_xy_x_y(g, nsxmin, nsxmax, nsymin, nsymax)
	end
	return t
end

CartesianIndices(fourier::FourierTransform{T,X}, field::FieldSpaceScalar) where {T, X<:AbstractVector} = Base.CartesianIndices((1, length(fourier.nsx_X), length(fourier.nsy_Y)))
CartesianIndices(fourier::FourierTransform{T,X}, field::FieldAngularSpectrumScalar) where {T, X<:AbstractVector} = Base.CartesianIndices((1, length(fourier.x_X), length(fourier.y_Y)))

function get_scatteringmatrixtype(fourier::FourierTransform{T,Y,R}, fieldi::FieldSpaceScalar{T,D,X,B}) where {T,D,Y,X,R,B}
	sizeI = length(fieldi.e_SXY)
	sizeS = length(CartesianIndices(fourier, fieldi))
	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeS, sizeI)
	t21 = Matrix{Complex{T}}(undef, sizeI, sizeS)
	m = Vector{Complex{T}}(undef,sizeS)
	ref1 = R <: ReferenceFrame ? fourier.ref : fieldi.ref
	nsx_X = Y <: AbstractVector ? fourier.nsx_X : error("not done")
	nsy_Y = Y <: AbstractVector ? fourier.nsy_Y : error("not done")
	if dir(fieldi) > 0
		fieldl = FieldSpaceScalar{T,-1,X,B}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrumScalar{T,1,Y,B}(deepcopy(nsx_X), deepcopy(nsy_Y), m, fieldi.λ, fieldi.n, ref1)
		return ScatteringMatrix{T,FieldSpaceScalar{T,-1,X,B}, FieldAngularSpectrumScalar{T,1,Y,B},Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldAngularSpectrumScalar{T,-1,Y,B}(deepcopy(nsx_X), deepcopy(nsy_Y), m, fieldi.λ, fieldi.n, ref1)
		fieldr = FieldSpaceScalar{T,1,X,B}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldAngularSpectrumScalar{T,-1,Y,B}, FieldSpaceScalar{T,1,X,B}, Nothing, Matrix{Complex{T}}}(r12, t21, r21, t12, fieldl, fieldr)
	end
end

function get_scatteringmatrixtype(fourier::FourierTransform{T,Y,R}, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,D,Y,X,R,B}
	sizeI = length(fieldi.e_SXY)
	sizeS = length(CartesianIndices(fourier, fieldi))
	r12 = nothing
	r21 = nothing
	t12 = Matrix{Complex{T}}(undef, sizeS, sizeI)
	t21 = Matrix{Complex{T}}(undef, sizeI, sizeS)
	m = Vector{Complex{T}}(undef,sizeS)
	ref1 = R <: ReferenceFrame ? fourier.ref : fieldi.ref
	x_X = Y <: AbstractVector ? fourier.x_X : error("not done")
	y_Y = Y <: AbstractVector ? fourier.y_Y : error("not done")
	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpaceScalar{T,1,Y,B}(deepcopy(x_X), deepcopy(y_Y), m, fieldi.λ, fieldi.n, ref1)
		return ScatteringMatrix{T,FieldAngularSpectrumScalar{T,-1,X,B}, FieldSpaceScalar{T,1,Y,B},Nothing, Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
	else
		fieldl = FieldSpaceScalar{T,-1,Y,B}(deepcopy(x_X), deepcopy(y_Y), m, fieldi.λ, fieldi.n, ref1)
		fieldr = FieldAngularSpectrumScalar{T,1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		return ScatteringMatrix{T, FieldSpaceScalar{T,-1,Y,B}, FieldAngularSpectrumScalar{T,1,X,B}, Nothing, Matrix{Complex{T}}}(r12, t21, r21, t12, fieldl, fieldr)
	end
end

function getfields_lr(fourier::FourierTransform{T,Y,R}, fieldi::FieldSpaceScalar{T,D,X,B}) where {T,D,Y,X,R,B}
	sizeS = length(CartesianIndices(fourier, fieldi))

	m = Vector{Complex{T}}(undef,sizeS)
	ref1 = R <: ReferenceFrame ? fourier.ref : fieldi.ref
	nsx_X = Y <: AbstractVector ? fourier.nsx_X : error("not done")
	nsy_Y = Y <: AbstractVector ? fourier.nsy_Y : error("not done")
	if dir(fieldi) > 0
		fieldl = FieldSpaceScalar{T,-1,X,B}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldAngularSpectrumScalar{T,1,Y,B}(deepcopy(nsx_X), deepcopy(nsy_Y), m, fieldi.λ, fieldi.n, ref1)
	else
		fieldl = FieldAngularSpectrumScalar{T,-1,Y,B}(deepcopy(nsx_X), deepcopy(nsy_Y), m, fieldi.λ, fieldi.n, ref1)
		fieldr = FieldSpaceScalar{T,1,X,B}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
	end
	return (fieldl, fieldr)
end

function getfields_lr(fourier::FourierTransform{T,Y,R}, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,D,Y,X,R,B}
	sizeS = length(CartesianIndices(fourier, fieldi))

	m = Vector{Complex{T}}(undef,sizeS)
	ref1 = R <: ReferenceFrame ? fourier.ref : fieldi.ref
	x_X = Y <: AbstractVector ? fourier.x_X : error("not done")
	y_Y = Y <: AbstractVector ? fourier.y_Y : error("not done")
	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
		fieldr = FieldSpaceScalar{T,1,Y,B}(deepcopy(x_X), deepcopy(y_Y), deepcopy(m), fieldi.λ, fieldi.n, ref1)
	else
		fieldl = FieldSpaceScalar{T,-1,Y,B}(deepcopy(x_X), deepcopy(y_Y), m, fieldi.λ, fieldi.n, ref1)
		fieldr = FieldAngularSpectrumScalar{T,1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, fieldi.ref)
	end
	return (fieldl, fieldr)
end

function coefficient_general(fourier::FourierTransform{T,X}, field::AbstractFieldMonochromatic) where {T, X <: AbstractVector}
	checkapplicability(fourier, field) || tobedone()

	scat = get_scatteringmatrixtype(fourier, field)
	cart_i = CartesianIndices(field)
	fieldt = dir(field) > 0 ? scat.fieldr : scat.fieldl
	cart_s = CartesianIndices(fieldt)

	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let scat=scat, fourier = fourier, field = field
		@inbounds Threads.@threads for i_s in iterator_index(fieldt)
			for i_i in iterator_index(field)
				if dir(field) > 0
					scat.t₁₂[i_s, i_i] = t(fourier, field, cart_i[i_i][2], cart_i[i_i][3], cart_s[i_s][2], cart_s[i_s][3])
					scat.t₂₁[i_i, i_s] = tinv(fourier, field, cart_i[i_i][2], cart_i[i_i][3], cart_s[i_s][2], cart_s[i_s][3])
				else
					scat.t₁₂[i_i, i_s] = tinv(fourier, field, cart_i[i_i][2], cart_i[i_i][3], cart_s[i_s][2], cart_s[i_s][3])
					scat.t₂₁[i_s, i_i] = t(fourier, field, cart_i[i_i][2], cart_i[i_i][3], cart_s[i_s][2], cart_s[i_s][3])
				end
			end
		end
	end
	return scat
end

function lightinteraction(fourier::FourierTransform{T,X}, field::AbstractFieldMonochromatic) where {T, X <: AbstractVector}
	checkapplicability(fourier, field) || tobedone()

	(fieldl, fieldr) = getfields_lr(fourier, field)
	fieldl.e_SXY .= zero(Complex{T})
	fieldr.e_SXY .= zero(Complex{T})

	fieldt = dir(field) > 0 ? fieldr : fieldl
	cart_i = CartesianIndices(field)
	cart_s = CartesianIndices(fieldt)

	#Work around to make good compiled code - https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373
	let fieldl = fieldl, fieldr = fieldr, fourier = fourier, field = field
		@inbounds Threads.@threads for i_s in iterator_index(fieldt)
			for i_i in iterator_index(field)
				fieldt.e_SXY[i_s] += t(fourier, field, cart_i[i_i][2], cart_i[i_i][3], cart_s[i_s][2], cart_s[i_s][3]) * field.e_SXY[i_i]
			end
		end
	end
	return (fieldl, fieldr)
end

@inline function coefficient_general(fourier::FourierTransform{T,X}, field::FieldSpaceScalar{T,D,Y}) where {T, X<:Nothing, Y <: AbstractRange, D}
	(sizeX, sizeY) = size(CartesianIndices(field))[2:3]
	Δx = Δvector(field.x_X, 1)
	Δy = Δvector(field.y_Y, 1)
	nsx_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δx)) * field.λ;
	nsy_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δy)) * field.λ;
	fourieraux = FourierTransform{T,Y}(field.x_X, field.y_Y, nsx_X, nsy_Y)
	return coefficient_general(fourieraux, field)
end

@inline function coefficient_general(fourier::FourierTransform{T,X}, field::FieldAngularSpectrumScalar{T,D,Y}) where {T,D, X <: Nothing, Y <: AbstractRange}
	(sizeX, sizeY) = size(CartesianIndices(field))[2:3]
	Δnsx = Δvector(field.nsx_X, 1)
	Δnsy = Δvector(field.nsy_Y, 1)
	x_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δnsx)) * field.λ;
	y_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δnsy)) * field.λ;
	fourieraux = FourierTransform{T,Y}(x_X, y_Y, field.nsx_X, field.nsy_Y)
	return coefficient_general(fourieraux, field)
end

function FieldSpaceScalar_fromAngularSpectrum(angspe::FieldAngularSpectrumScalar, x, y)
	fourier = FourierTransform(x,y,angspe.nsx_X, angspe.nsy_Y)
	return lightinteraction(fourier, angspe)
end

function FieldSpaceVectorial_fromAngularSpectrum(angspe::FieldAngularSpectrumVectorial, x, y)
	fourier = FourierTransform(x,y,angspe.nsx_X, angspe.nsy_Y)
	return lightinteraction(fourier, angspe)
end

function FieldAngularSpectrumScalar_fromSpace(space::FieldSpaceScalar, nsx, nsy)
	fourier = FourierTransform(space.x, space.y, nsx, nsy)
	return lightinteraction(fourier, space)
end

function FieldAngularSpectrumVectorial_fromSpace(space::FieldSpaceVectorial, nsx, nsy)
	fourier = FourierTransform(space.x, space.y, nsx, nsy)
	return lightinteraction(fourier, space)
end

#=
function coefficient_specific(fourier::FourierTransform{T,X}, field::F) where {T, X <: Nothing, F <: FieldSpaceScalar{T, <:AbstractRange}}
	(sizeX, sizeY) = size(CartesianIndices(field))[2:3]
	planfft = FFTW.plan_fft(reshape(field.e_SXY, sizeX, sizeY))

	(Δx, Δy) = (Δvector(field.x_X, 1), Δvector(field.y_Y, 1))

	nsx_X = fftshift(FFTW.fftfreq(sizeX, 1 / Δx)) * field.λ;
	nsy_Y = fftshift(FFTW.fftfreq(sizeY, 1 / Δy)) * field.λ;
	scalar = Δx * Δy / 4 / π^2

	(Δnsx, Δnsy) = (Δvector(nsx_X, 1), Δvector(nsy_Y, 1))
	i_scalar = sizeX * sizeY * Δnsx * Δnsy * 4 * π^2 / field.λ^2

	tmp = Array{Complex{T},2}(undef, sizeX, sizeY)
	if dir(field) > 0
		fieldr = FieldAngularSpectrumScalar{T, 1}(nsx_X, nsy_Y, field.e_SXY, field.λ, field.n, field.ref)
		return FFTCoefficient{T,F,typeof(fieldr),typeof(planfft)}(planfft, scalar, i_scalar, tmp, field, fieldr)
	else
		fieldl = FieldAngularSpectrumScalar{T,-1}(nsx_X, nsy_Y, field.e_SXY, field.λ, field.n, field.ref)
		return FFTCoefficient{T,typeof(fieldl),F,typeof(planfft)}(planfft, scalar, i_scalar, tmp, fieldl, field)
	end
end

function coefficient_specific(fourier::FourierTransform{T,X}, field::F) where {T, X <: Nothing, F <: FieldAngularSpectrumScalar{T, <:AbstractRange}}
	(sizeX, sizeY) = size(CartesianIndices(field))[2:3]
	planfft = FFTW.plan_fft(reshape(field.e_SXY, sizeX, sizeY))
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
=#

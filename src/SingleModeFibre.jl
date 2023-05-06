mutable struct SingleModeFibre{T} <: AbstractOpticalComponent{T}
	mfd::T
	n1::Material{T}
	dir1::Int8
	ref1::ReferenceFrame{T}
	length::T
	n2::Material{T}
	dir2::Int8
	ref2::ReferenceFrame{T}
	function SingleModeFibre{T}(mfd, n1, dir1 = 1, ref1 = ReferenceFrame(), length=0, n2=n1, dir2=dir1, ref2=ref1) where T
		return new{T}(mfd, n1, dir1 >= 0 ? 1 : -1, ref1, length, n2, dir2 >= 0 ? 1 : -1, ref2);
	end
end
SingleModeFibre(mfd, n1, dir1 = 1, ref1 = ReferenceFrame(), length=0, n2=n1, dir2=dir1, ref2=ref1) = SingleModeFibre{Float64}(mfd, n1, dir1, ref1, length, n2, dir2, ref2)

function FieldSpaceScalar_fromfibre(fibre::SingleModeFibre, x_X, y_Y, λ, tip=1, intensity=1.)
	tip == 2 ? (n = fibre.n2; dir2 = fibre.dir2; ref = fibre.ref2) : (n = fibre.n1; dir2 = fibre.dir1; ref = fibre.ref1);
	fieldspace = FieldSpaceScalar_gaussian(x_X, y_Y, fibre.mfd, λ, n(λ), Int64(dir), ref);
	fieldspace.e_SXY .*= √intensity
	return fieldspace;
end

# function FieldSpaceSymmetric_fromfibre(fibre::SingleModeFibre, x_X, λ, tip=1, intensity=1.)
# 	tip == 2 ? (n = fibre.n2; dir2 = fibre.dir2; ref = fibre.ref2) : (n = fibre.n1; dir2 = fibre.dir1; ref = fibre.ref1);
# 	fieldspace = FieldSpaceSymmetric_gaussian(x_X, fibre.mfd, λ, n(λ), Int64(dir2), ref);
# 	fieldspace.e_SXY = fieldspace.e_SXY .* sqrt(intensity);
# 	return fieldspace;
# end

function FieldAngularSpectrumScalar_fromfibre(fibre::SingleModeFibre{T}, nsx_X, nsy_Y, λ, tip=1, intensity=1.) where T
	tip == 2 ? (n = fibre.n2; dir2 = fibre.dir2; ref = fibre.ref2) : (n = fibre.n1; dir2 = fibre.dir1; ref = fibre.ref1);

	angspe = FieldAngularSpectrumScalar_gaussian(nsx_X, nsy_Y, fibre.mfd, λ, n(λ), Int64(dir2), ref);
	angspe.e_SXY .*= √intensity
	return angspe
end

function FieldAngularSpectrumScalarRadialSymmetric_fromfibre(fibre::SingleModeFibre{T}, nsr_R, λ, tip=1, intensity=1.) where T
	tip == 2 ? (n = fibre.n2; dir2 = fibre.dir2; ref = fibre.ref2) : (n = fibre.n1; dir2 = fibre.dir1; ref = fibre.ref1);

	angspe = FieldAngularSpectrumScalarRadialSymmetric_gaussian(nsr_R, fibre.mfd, λ, n(λ), Int64(dir2), ref);
	angspe.e_SXY .*= √intensity
	return angspe
end

# function FieldAngularSpectrumSymmetric_fromfibre(fibre::SingleModeFibre, sx_X::AbstractVector{<:Number}, λ::Real, tip::Integer=1, intensity::Real=1.)::FieldAngularSpectrumSymmetric
# 	tip == 2 ? (n = fibre.n2; dir2 = fibre.dir2; ref = fibre.ref2) : (n = fibre.n1; dir2 = fibre.dir1; ref = fibre.ref1);
#
# 	angspe = FieldAngularSpectrumSymmetric_gaussian(sx_X, fibre.mfd, λ, n(λ), Int64(dir), ref);
# 	angspe.e_SXY *= √(intensity);
# 	return angspe
# end

function signal_complex(fibre::SingleModeFibre, angspe::FieldAngularSpectrumScalar{T}, tip=1) where T
	tip == 2 ? (ref = fibre.ref2; dir2 = fibre.dir2;) : (ref = fibre.ref1; dir2 = fibre.dir1;)
	dir2 == dir(angspe) ? error("The direction that the fibre tip is pointing at must be the oposite from the field propagation (fibre.dir must be equal to -angspe.dir)") : nothing
	angsperef = changereferenceframe(angspe, ref);
	k = 2 * π / angsperef.λ;
	cons = fibre.mfd * √(1 / 32 / π^3) * k^2;
	int = zero(Complex{T})
	cart = CartesianIndices(angsperef)
	@inbounds @simd for i in iterator_index(angsperef)
		(xmin, xmax) = integralExtremes(angsperef.nsx_X, cart[i][2])
		(ymin, ymax) = integralExtremes(angsperef.nsy_Y, cart[i][3])
		sens = exp(-fibre.mfd^2 / 16 * k^2 * (angsperef.nsx_X[cart[i][2]]^2 + angsperef.nsy_Y[cart[i][3]]^2))
		int += angsperef.e_SXY[i] * sens * (xmax - xmin) * (ymax - ymin)
	end
	return √(16π^4 * real(angsperef.n)) * cons * int
end

function signal_complex(fibre::SingleModeFibre, angspe::FieldAngularSpectrumScalarRadialSymmetric{T}, tip=1) where T
	tip == 2 ? (ref = fibre.ref2; dir2 = fibre.dir2;) : (ref = fibre.ref1; dir2 = fibre.dir1;)
	dir2 == dir(angspe) ? error("The direction that the fibre tip is pointing at must be the oposite from the field propagation (fibre.dir must be equal to -angspe.dir)") : nothing
	angsperef = changereferenceframe(angspe, ref);
	k = 2 * π / angsperef.λ;
	cons = fibre.mfd * √(1 / 32 / π^3) * k^2;
	int = zero(Complex{T})
	cart = CartesianIndices(angsperef)
	@inbounds @simd for i in iterator_index(angsperef)
		(rmin, rmax) = integralExtremes(angsperef.nsr_R, cart[i][2])
		sens = exp(-fibre.mfd^2 / 16 * k^2 * angsperef.nsr_R[cart[i][2]]^2)
		int += angsperef.e_SXY[i] * sens * (rmax - rmin) * angsperef.nsr_R[cart[i][2]]
	end
	return √(64π^6 * real(angsperef.n)) * cons * int
end

function signal_complex(fibre::SingleModeFibre, fieldspace::FieldSpaceScalar{T}, tip::Integer=1) where T
	tip == 2 ? (ref = fibre.ref2; dir2 = fibre.dir2;) : (ref = fibre.ref1; dir2 = fibre.dir1);
	dir2 == dir(fieldspace) && error("The direction that the fibre tip is pointing at must be the oposite from the field propagation (fibre.dir must be equal to -fieldspace.dir)")

	fieldspaceref = changereferenceframe(fieldspace, ref);

	int = zero(Complex{T})
	cart = CartesianIndices(fieldspaceref)
	@inbounds @simd for i in iterator_index(fieldspaceref)
		(xmin, xmax) = integralExtremes(fieldspaceref.x_X, cart[i][2])
		(ymin, ymax) = integralExtremes(fieldspaceref.y_Y, cart[i][3])
		sens = exp(- 4 / fibre.mfd^2 * (fieldspaceref.x_X[cart[i][2]]^2 + fieldspaceref.y_Y[cart[i][3]]^2))
		int += fieldspaceref.e_SXY[i] * sens * (xmin - xmax) * (ymax - ymin)
	end
	return √(8 / π / fibre.mfd^2) * int
end

function signal(fibre::SingleModeFibre, field, tip::Integer = 1)
	return abs2(signal_complex(fibre, field, tip))
end
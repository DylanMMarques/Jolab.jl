mutable struct SingleModeFibre{T} <: AbstractOpticalComponent{T}
	mfd::T
	n1::JolabFunction1D{T,Complex{T}}
	dir1::Int8
	ref1::ReferenceFrame{T}
	length::T
	n2::JolabFunction1D{T,Complex{T}}
	dir2::Int8
	ref2::ReferenceFrame{T}
	function SingleModeFibre{T}(mfd, n1, dir1 = 1, ref1 = ReferenceFrame(), length=0, n2=n1, dir2=dir1, ref2=ref1) where T
		return new{T}(mfd, n1, dir1 >= 0 ? 1 : -1, ref1, length, n2, dir2 >= 0 ? 1 : -1, ref2);
	end
end
SingleModeFibre(mfd, n1, dir1 = 1, ref1 = ReferenceFrame(), length=0, n2=n1, dir2=dir1, ref2=ref1) = SingleModeFibre{Float64}(mfd, n1, dir1, ref1, length, n2, dir2, ref2)

function FieldSpace_fromfibre(fibre::SingleModeFibre, x_X::AbstractVector{<:Real}, y_Y::AbstractVector{<:Real}, λ::Real, tip::Integer=1, intensity::Real=1.)::FieldSpace
	tip == 2 ? (n = fibre.n2; dir = fibre.dir2; ref = fibre.ref2) : (n = fibre.n1; dir = fibre.dir1; ref = fibre.ref1);
	fieldspace = FieldSpace_gaussian(x_X, y_Y, fibre.mfd, λ, n(λ), dir, ref);
	fieldspace.e_SXY = fieldspace.e_SXY .* sqrt(intensity);
	return fieldspace;
end

function FieldSpaceSymmetric_fromfibre(fibre::SingleModeFibre, x_X::AbstractVector{<:Real}, λ::Real, tip::Integer=1, intensity::Real=1.)::FieldSpaceSymmetric
	tip == 2 ? (n = fibre.n2; dir = fibre.dir2; ref = fibre.ref2) : (n = fibre.n1; dir = fibre.dir1; ref = fibre.ref1);
	fieldspace = FieldSpaceSymmetric_gaussian(x_X, fibre.mfd, λ, n(λ), dir, ref);
	fieldspace.e_SXY = fieldspace.e_SXY .* sqrt(intensity);
	return fieldspace;
end

function FieldAngularSpectrum_fromfibre(fibre::SingleModeFibre, sx_X::AbstractVector{<:Number}, sy_Y::AbstractVector{<:Number}, λ::Real, tip::Integer=1, intensity::Real=1.)::FieldAngularSpectrum
	tip == 2 ? (n = fibre.n2; dir = fibre.dir2; ref = fibre.ref2) : (n = fibre.n1; dir = fibre.dir1; ref = fibre.ref1);

	angspe = FieldAngularSpectrum_gaussian(sx_X, sy_Y, fibre.mfd, λ, n(λ), dir, ref);
	angspe.e_SXY *= sqrt(intensity);
	return angspe
end

function FieldAngularSpectrumSymmetric_fromfibre(fibre::SingleModeFibre, sx_X::AbstractVector{<:Number}, λ::Real, tip::Integer=1, intensity::Real=1.)::FieldAngularSpectrumSymmetric
	tip == 2 ? (n = fibre.n2; dir = fibre.dir2; ref = fibre.ref2) : (n = fibre.n1; dir = fibre.dir1; ref = fibre.ref1);

	angspe = FieldAngularSpectrumSymmetric_gaussian(sx_X, fibre.mfd, λ, n(λ), dir, ref);
	angspe.e_SXY *= √(intensity);
	return angspe
end

function signal(fibre::SingleModeFibre, angspe::FieldAngularSpectrum, tip::Integer=1)::Float64

	tip == 2 ? (ref = fibre.ref2; dir = fibre.dir2;) : (ref = fibre.ref1; dir = fibre.dir1;)

	dir == angspe.dir ? error("The direction that the fibre tip is pointing at must be the oposite from the field propagation (fibre.dir must be equal to -angspe.dir)") : nothing
	angsperef = changereferenceframe(angspe, ref);

	k = 2 * π / angsperef.λ;
	kx_X = reshape(k * angsperef.nsx_X, 1, :);
	ky_Y = reshape(k * angsperef.nsy_Y, 1, 1, :);

	cons = fibre.mfd * √(1 / 32 / π^3);
	sens_XY = cons * exp.(.- fibre.mfd.^2 ./ 16 .* (kx_X.^2 .+ ky_Y.^2));

	if size(angsperef.e_SXY, 1) == 3
		e_SXY = angsperef.e_SXY[1:2,:,:] .* sens_XY;
		e_S = [∫∫(e_SXY[i,:,:], vec(kx_X), vec(ky_Y)) for i in 1:2]
		return 16 * π^4 * dotdim(e_S, e_S, 1);
	else
		e_SXY = angsperef.e_SXY .* sens_XY;
		e_S = ∫∫(e_SXY[1,:,:], vec(kx_X), vec(ky_Y));
		return 16 * π^4 * abs2(e_S);
	end
end

function signal(fibre::SingleModeFibre, angspe::FieldAngularSpectrumSymmetric, tip::Integer=1)::Float64

	tip == 2 ? (ref = fibre.ref2; dir = fibre.dir2;) : (ref = fibre.ref1; dir = fibre.dir1;)

	dir == angspe.dir ? error("The direction that the fibre tip is pointing at must be the oposite from the field propagation (fibre.dir must be equal to -angspe.dir)") : nothing
	angsperef = angspe#changereferenceframe(angspe, ref);

	k = 2 * π * angsperef.n / angsperef.λ;
	kx_X = adddims(k * angsperef.sx_X, (1,));

	cons = fibre.mfd * √(1 / 32 / π^3);
	sens_XY = cons .* exp.(.- fibre.mfd.^2 ./ 16 .* kx_X.^2);

	e_SXY = angsperef.e_SXY .* sens_XY .* kx_X;
	e_S = ∫(vec(e_SXY), vec(kx_X));
	return 64 .* π.^6 .* abs2(e_S);
end

function signal(fibre::SingleModeFibre, fieldspace::FieldSpace, tip::Integer=1)::Float64
	tip == 2 ? (ref = fibre.ref2; dir = fibre.dir2;) : (ref = fibre.ref1; dir = fibre.dir1);
	dir == fieldspace.dir ? error("The direction that the fibre tip is pointing at must be the oposite from the field propagation (fibre.dir must be equal to -fieldspace.dir)") : nothing

	fieldspaceref = changereferenceframe(fieldspace, ref);

	x_X = adddims(fieldspaceref.x_X, (1,));
	y_Y = adddims(fieldspaceref.y_Y, (1,2));

	cons = √(8 / π / fibre.mfd.^2);
	sens_XY = cons .* exp.(.- 4 ./ fibre.mfd.^2 .* (x_X.^2 .+ y_Y.^2));

	if size(fieldspaceref.e_SXY, 1) == 3
		e_SXY = fieldspaceref.e_SXY[1:2,:,:] .* sens_XY;
		e_S = [∫∫(e_SXY[i,:,:], vec(x_X), vec(y_Y)) for i in 1:2]
		return dotdim(e_S, e_S, 1);
	else
		e_SXY = fieldspaceref.e_SXY .* sens_XY;
		e_S = ∫∫(e_SXY[1,:,:], vec(x_X), vec(y_Y));
		return abs2(e_S);
	end
end

function signal(fibre::SingleModeFibre, fieldspace::FieldSpaceSymmetric, tip::Integer=1)::Float64
	tip == 2 ? (ref = fibre.ref2; dir = fibre.dir2;) : (ref = fibre.ref1; dir = fibre.dir1);
	dir == fieldspace.dir ? error("The direction that the fibre tip is pointing at must be the oposite from the field propagation (fibre.dir must be equal to -fieldspace.dir)") : nothing

	fieldspaceref = changereferenceframe(fieldspace, ref);

	x_X = adddims(fieldspaceref.x_X, (1,));

	cons = √(8 / π / fibre.mfd.^2);
	sens_XY = cons .* exp.(.- 4 ./ fibre.mfd.^2 .* x_X.^2);

	e_SXY = fieldspaceref.e_SXY .* sens_XY .* x_X;
	e_S = ∫(vec(e_SXY), vec(x_X));
	return abs2(e_S) * 4 * π^2;
end

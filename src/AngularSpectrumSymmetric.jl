mutable struct FieldAngularSpectrumSymmetric{T,X<:AbstractArray{T}} <: AbstractFieldAngularSpectrum{T}
	sx_X::X
	sy_Y::X
	e_SXY::Array{Complex{T}, 3}
	λ::T
	n::Complex{T}
	dir::Int8
	ref::ReferenceFrame{T}
	function FieldAngularSpectrumSymmetric{T}(sx_X::X, e_SXY, λ, n, dir, ref) where {T,X}
		length(sx_X) != size(e_SXY, 2) ? error("The length of sx_X must be the same as the size of e_SXY in the second dimension.") : nothing;
		size(e_SXY, 3) != 1 ? error("The length of e_SXY must be 1 in the third dimension.") : nothing;
		abs(ref.x) > 1E-10 || abs(ref.y) > 1E-10 || abs(ref.θ) > 1E-5 || abs(ref.ϕ) > 1E-5 ? error("The referenceframe of a symetric field must be with x, y, θ and ϕ equal to 0") : nothing;
		size(e_SXY, 1) != 1 ? error("The size of e_SXY of a symmetric field must be 1 (scallar field).") : nothing;
		return new{T,X}(sx_X, [0.], e_SXY, λ, n, dir >= 0 ? 1 : -1, ref);
	end
end

function FieldAngularSpectrumSymmetric_gaussian(sx_X::AbstractVector{<:Number}, ω::Real, λ::Real, n::Number=1. + 0 .* im, dir::Int = 1, ref::ReferenceFrame=ReferenceFrame())::FieldAngularSpectrumSymmetric
	k = 2 * π * n / λ;
	norm = 	ω * √(1 / 32 / π^3);
	sx_X = reshape(sx_X, 1, :, 1);
	e_SXY = norm .* exp.(- ω.^2 .* k.^2 ./ 16 .* sx_X.^2);
	sx_X = dropdims(sx_X, dims=(1, 3));
	return FieldAngularSpectrumSymmetric(sx_X, e_SXY, λ, n, dir, ref);
end

function changereferenceframe!(angspe::FieldAngularSpectrumSymmetric, refnew::ReferenceFrame)
	#Needs checking after doing the 2D interpolation on C
	abs(angspe.ref.x - refnew.x) > 1E-10 || abs(angspe.ref.y - refnew.y) > 1E-10 || abs(angspe.ref.θ - refnew.θ) > 1E-5 || abs(angspe.ref.ϕ - refnew.ϕ) > 1E-5 ? error("FieldAngularSpectrumSymmetric cannot be translate in x or y or rotated") : nothing

	if !checkposition(angspe.ref, refnew)
	 	translatereferenceframe!(angspe, refnew);
	end
end

function translatereferenceframe!(angspe::FieldAngularSpectrumSymmetric, refnew::ReferenceFrame)
	refΔz = refnew.z - angspe.ref.z;

	k = 2 * π / angspe.λ;
	kx_X = k .* reshape(angspe.nsx_X, 1, :, 1);
	kz_Z = angspe.dir .* .√( k.^2 .- kx_X.^2);

	angspe.e_SXY = angspe.e_SXY .* exp.(im .* (kz_Z .* refΔz));
	angspe.ref.z = refnew.z;
end

function intensity(angspe::FieldAngularSpectrumSymmetric)::Float64
	k = 2 .* π / angspe.λ;
	kx_X = real.(k .* angspe.nsx_X);

	e_SXY = abs2.(angspe.e_SXY) .* adddims(kx_X, (1,));
	return 4 * π^2 * 2 * π * angspe.n * ∫(e_SXY, kx_X);
end

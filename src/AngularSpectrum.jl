mutable struct FieldAngularSpectrum{T,X<:AbstractVector{T}} <: AbstractFieldAngularSpectrum{T}
	nsx_X::X
	nsy_Y::X
	e_SXY::Array{Complex{T}, 3}
	λ::T
	n::Complex{T}
	dir::Int8
	ref::ReferenceFrame{T}
	function FieldAngularSpectrum{T,X}(nsx_X, nsy_Y, e_SXY, λ, n, dir, ref) where {T,X}
		length(nsx_X) != size(e_SXY, 2) && error("The length of sx_X must be the same as the size of e_SXY in the second dimension.")
		length(nsy_Y) != size(e_SXY, 3) && error("The length of sy_Y must be the same as the size of e_SXY in the third dimension.")
		(size(e_SXY, 1) != 1 && size(e_SXY, 1) != 3) && error("The size of e_SXY must be 1 (scallar field) or 3 (vectorial field).")
		return new{T,X}(nsx_X, nsy_Y, e_SXY, λ, n, dir >= 0 ? 1 : -1, ref);
	end
	function FieldAngularSpectrum{T}(nsx_X::X, nsy_Y::Y, e_SXY, λ, n, dir, ref) where {T,X,Y}
		M = promote_type(X,Y)
		return FieldAngularSpectrum{T,M}(nsx_X, nsy_Y, e_SXY, λ, n, dir, ref)
	end
end

function FieldAngularSpectrum_uniform(nsx_X, nsy_Y, λ, n, dir, ref)
	T = Float64
	e_SXY = ones(Complex{T}, 1, length(nsx_X), length(nsy_Y))
	angspe = FieldAngularSpectrum{T}(nsx_X, nsy_Y, e_SXY, λ, n, dir, ref);
	int = √(intensity(angspe))
	vec(angspe.e_SXY) ./= int
	return angspe
end

# FieldAngularSpectrum(nsx_X, nsy_Y, e_SXY, λ, n, dir, ref) = FieldAngularSpectrum{Float64}(nsx_X, nsy_Y, e_SXY, λ, n, dir, ref)
function FieldAngularSpectrum_gaussian(nsx_X, nsy_Y, ω, λ, n, dir, ref)
	T = Float64
	k = 2 * π / λ;
	norm = 	ω * √(1 / 32 / π^3);
	e_SXY = Array{Complex{T},3}(undef,1,length(nsx_X), length(nsy_Y))
	@inbounds @simd for iY in eachindex(nsy_Y)
		for iX in eachindex(nsx_X)
			e_SXY[1,iX,iY] = norm * exp((-ω^2 * k^2 / 16) * (nsx_X[iX]^2 + nsy_Y[iY]^2))
		end
	end
	return FieldAngularSpectrum{T}(nsx_X, nsy_Y, e_SXY, λ, n, dir, ref);
end

function FieldSpace_fromangspe(angspe::FieldAngularSpectrum{T}, x_X, y_Y) where {T}
	e_SXY = inversefourriertransform(angspe.nsx_X, angspe.nsy_Y, angspe.e_SXY, angspe.λ, angspe.n, x_X, y_Y);
	return FieldSpace{T}(x_X, y_Y, e_SXY, angspe.λ, angspe.n, angspe.dir, angspe.ref);
end

function FieldSpace_fromangspefft(angspe::FieldAngularSpectrum{T,X}; padding=0::Integer) where {T, X<:AbstractRange{T}}
	(x_X, y_Y, e_SXY) = inversefourriertransformfft(angspe.nsx_X, angspe.nsy_Y, angspe.e_SXY, angspe.λ, padding=padding);
	return FieldSpace{T}(x_X, y_Y, e_SXY, angspe.λ, angspe.n, angspe.dir, angspe.ref);
end

function angspeto3Dspace(angspe::FieldAngularSpectrum, x_X::AbstractVector{<:Real}, y_Y::AbstractVector{<:Real}, z_Z::AbstractVector{<:Real})
	e_SXYZ = zeros(eltype(angspe.e_SXY), size(angspe.e_SXY, 1), length(x_X), length(y_Y), length(z_Z));
 	angspeWithZ = Array{eltype(angspe.e_SXY), 3}(undef, size(angspe.e_SXY));
	kz_XY = (2π / angspe.λ) .* .√(angspe.n^2 .- reshape(angspe.nsx_X.^2, 1,:,1) .- reshape(angspe.nsy_Y.^2, 1, 1, :))
	@inbounds @simd for iZ in 1:length(z_Z)
		angspeWithZ .= angspe.e_SXY .* exp.(im .* kz_XY .* z_Z[iZ]);
		auxE = @view e_SXYZ[:,:,:,iZ];
		inversefourriertransform!(auxE, angspe.nsx_X, angspe.nsy_Y, angspeWithZ, angspe.λ, angspe.n, x_X, y_Y);
	end
	return e_SXYZ;
end

function changereferenceframe!(nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number,2}, nsz_XY::AbstractArray{<:Number}, e_SXY::AbstractArray{<:Number, 3}, λ::Real, refold::ReferenceFrame, refnew::ReferenceFrame)
	checkposition(refold, refnew) || translatereferenceframe!(e_SXY, nsx_XY, nsy_XY, nsz_XY, λ, refold, refnew)
	checkorientation(refold, refnew) || rotatereferenceframe!(nsx_XY, nsy_XY, nsz_XY, e_SXY, refold, refnew)
end

function propagationmatrix!(propMatrix::AbstractArray{<:Number, 2}, nsx_XY::AbstractArray{<:Number}, nsy_XY::AbstractArray{<:Number}, nsz_XY::AbstractArray{<:Number}, λ::Real, refold::ReferenceFrame, refnew::ReferenceFrame)
	checkorientation(refnew, refold) || error("Cannot calculate propagation matrix as the referenceframe are not oriented")
	length(nsz_XY) == length(nsy_XY) == length(nsx_XY) == size(propMatrix, 1) || error("wrong sizes")
	refΔx = refnew.x - refold.x;
	refΔy = refnew.y - refold.y;
	refΔz = refnew.z - refold.z;

	# Can't use based on referenceframe as the rotation is inverted
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, refold.θ, refold.ϕ);
	k = im * 2π / λ;
	@inbounds @simd for i in eachindex(nsx_XY)
		tmpphase = exp(k * (nsx_XY[i] * refΔx + nsy_XY[i] * refΔy + nsz_XY[i] * refΔz))
		propMatrix.diag[i] = tmpphase;
	end
end

function propagationmatrix(fieldl::AbstractFieldAngularSpectrum{T}, ref::ReferenceFrame) where T
	checkorientation(fieldl.ref, ref) || error("Cannot calculate propagation matrix as the referenceframe are not oriented")
	refΔx = ref.x - fieldl.ref.x;
	refΔy = ref.y - fieldl.ref.y;
	refΔz = ref.z - fieldl.ref.z;

	# Can't use based on referenceframe as the rotation is inverted
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, fieldl.ref.θ, fieldl.ref.ϕ);
	imk = im * 2π / fieldl.λ;
	propMatrix = Diagonal(Vector{Complex{T}}(undef, length(fieldl.nsx_X) * length(fieldl.nsy_Y)))
	i = 1
	@inbounds for iY in eachindex(fieldl.nsy_Y)
		for iX in eachindex(fieldl.nsx_X)
			nsz = √(complex(fieldl.n^2 - fieldl.nsx_X[iX]^2 - fieldl.nsy_Y[iY]^2))
			propMatrix.diag[i] = exp(imk * (fieldl.nsx_X[iX] * refΔx + fieldl.nsy_Y[iY] * refΔy + nsz * refΔz))
			i += 1
		end
	end
	return propMatrix
end

@inline propagationmatrix(fieldl::AbstractFieldAngularSpectrum, fieldr::AbstractFieldAngularSpectrum) = propagationmatrix(fieldl, fieldr.ref)

function propagationmatrix!(propMatrix::AbstractArray{<:Number, 2}, nsx_X::AbstractVector{<:Number}, nsy_Y::AbstractVector{<:Number}, nsz_XY::AbstractArray{<:Number}, λ::Real, refold::ReferenceFrame, refnew::ReferenceFrame)
	checkorientation(refnew, refold) || error("Cannot calculate propagation matrix as the referenceframe are not oriented")
	length(nsx_X) * length(nsy_Y) == size(propMatrix,1) == length(nsz_XY) || error("Wrong sizes")
	refΔx = refnew.x - refold.x;
	refΔy = refnew.y - refold.y;
	refΔz = refnew.z - refold.z;

	# Can't use based on referenceframe as the rotation is inverted
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, refold.θ, refold.ϕ);
	k = im * 2π / λ;
	i = 1
	@inbounds for iY in eachindex(nsy_Y)
		for iX in eachindex(nsx_X)
			propMatrix.diag[i] = exp(k * (nsx_X[iX] * refΔx + nsy_Y[iY] * refΔy + nsz_XY[iX,iY] * refΔz))
			i += 1
		end
	end
end

function translatereferenceframe!(e_SXY::AbstractArray{<:Number, 3}, nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number,2}, nsz_XY::AbstractArray{<:Number,2}, λ::Real, refold::ReferenceFrame, refnew::ReferenceFrame)
	refΔx = refnew.x - refold.x;
	refΔy = refnew.y - refold.y;
	refΔz = refnew.z - refold.z;

	# Can't use based on referenceframe as the rotation is inverted
	(refΔx, refΔy, refΔz) = rotatecoordinatesfrom(refΔx, refΔy, refΔz, refold.θ, refold.ϕ);
	k = im * 2π / λ;
	tmpe_SXY = reshape(e_SXY, size(e_SXY, 1), :);
	@inbounds @simd for i in eachindex(nsx_XY)
		tmpphase = exp(k * (nsx_XY[i] * refΔx + nsy_XY[i] * refΔy + nsz_XY[i] * refΔz))
		for iS in 1:size(e_SXY,1)
				tmpe_SXY[iS, i] *= tmpphase;
		end
	end
end

function rotatereferenceframe!(nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number,2}, nsz_XY::AbstractArray{<:Number,2}, e_SXY::AbstractArray{<:Number,3}, refold::ReferenceFrame, refnew::ReferenceFrame)
	# Rotate the coordinates to the referenceframe (0, 0)
	(sum(imag.(nsz_XY)) > @tol) && error("You can't rotate the field in a media that absorves light")
	rotatecoordinatesfromto!(nsx_XY, nsy_XY, nsz_XY, refold, refnew)
	# Rotate the direction respective to the referenceframe refnew
	if size(e_SXY) == 3 #if vectorial need to rotate the field as well
		ex = @view e_SXY[1,:,:]; ey = @view e_SXY[2,:,:]; ez = @view e_SXY[3,:,:];
		rotatecoordinatatesfromto!(ex, ey, ez, refold, refnew);
	end
end

function changereferenceframe!(angspe::FieldAngularSpectrum, refnew::ReferenceFrame)
	#Needs checking after doing the 2D interpolation
	checkposition(angspe.ref, refnew) || translatereferenceframe!(angspe, refnew);
	checkorientation(angspe.ref, refnew) || rotatereferenceframe!(angspe, refnew);
end

function translatereferenceframe!(angspe::FieldAngularSpectrum, refnew::ReferenceFrame)
	nsx_XY = repeat(angspe.nsx_X, 1, length(angspe.nsy_Y))
	nsy_XY = repeat(reshape(angspe.nsy_Y, 1, :), length(angspe.nsx_X), 1)
	nsz_XY = angspe.dir .* .√(angspe.n^2 .- angspe.nsx_X.^2 .- (angspe.nsy_Y').^2)
	translatereferenceframe!(angspe.e_SXY, nsx_XY, nsy_XY, nsz_XY, angspe.λ, angspe.ref, refnew)
	angspe.ref.x = refnew.x;
	angspe.ref.y = refnew.y;
	angspe.ref.z = refnew.z;
end

function rotatereferenceframe!(angspe::FieldAngularSpectrum, refnew::ReferenceFrame)
	(sizeX, sizeY) = size(angspe.e_SXY)[2:3]
	nsx_XY = repeat(angspe.nsx_X, 1, length(angspe.nsy_Y))
	nsy_XY = repeat(reshape(angspe.nsy_Y, 1, :), length(angspe.nsx_X), 1)
	nsz_XY = angspe.dir .* .√(angspe.n^2 .- angspe.nsx_X.^2 .- (angspe.nsy_Y').^2)
	rotatereferenceframe!(nsx_XY, nsy_XY, nsz_XY, angspe.e_SXY, angspe.ref, refnew)

	if (nsx_XY[2] - nsx_XY[1]) > 0
		angspe.nsx_X = range(minimum(view(nsx_XY,1,:)), stop = maximum(view(nsx_XY,sizeX,:)), length = sizeX)
	else
		angspe.nsx_X = range(maximum(view(nsx_XY,1,:)), stop = minimum(view(nsx_XY,sizeX,:)), length = sizeX)
	end
	if (nsy_XY[1,2] - nsy_XY[1,1]) > 0
		angspe.nsy_Y = range(minimum(view(nsy_XY,:,1)), stop = maximum(view(nsy_XY,:,sizeY)), length = sizeY)
	else
		angspe.nsy_Y = range(maximum(view(nsy_XY,:,1)), stop = minimum(view(nsy_XY,:,sizeY)), length = sizeY)
	end

	if size(angspe.e_SXY) == 3 #if vectorial need to rotate the field as well
		eaux_XY = @view angspe.e_SXY[1,:,:]
		angspe.e_SXY[1,:,:] = interpolatetogrid(nsx_XY, nsy_XY, eaux_XY, angspe.nsx_X, angspe.nsy_Y);# .* √(cos(angspe.ref.θ) / cos(refnew.θ));
		eaux_XY = @view angspe.e_SXY[2,:,:]
		angspe.e_SXY[2,:,:] = interpolatetogrid(nsx_XY, nsy_XY, eaux_XY, angspe.nsx_X, angspe.nsy_Y);# .* √(cos(angspe.ref.θ) / cos(refnew.θ));
		eaux_XY = @view angspe.e_SXY[3,:,:]
		angspe.e_SXY[3,:,:] = interpolatetogrid(nsx_XY, nsy_XY, eaux_XY, angspe.nsx_X, angspe.nsy_Y);# .* √(cos(angspe.ref.θ) / cos(refnew.θ));
	else
		angspe.e_SXY = interpolatetogrid(nsx_XY, nsy_XY, angspe.e_SXY, angspe.nsx_X, angspe.nsy_Y);# .* √(cos(angspe.ref.θ) / cos(refnew.θ));
	end
	angspe.ref.θ = refnew.θ;
	angspe.ref.ϕ = refnew.ϕ;
end

@inline function polarizedcomponents(E_S::AbstractVector{<:Number}, nsx::Number, nsy::Number, dir::Integer, n::Number)
	sr = √(nsx^2 + nsy^2) / n;
	if abs2(sr) < 1E-20
		@show "this has to be chekect(polarizedcomponents of Angularspectrum)"
		return (dir * E_S[1], - E_S[2]);
	else
		sz = √(complex(1 - sr^2));
		return (dir * (nsx / n * E_S[1] + nsy / n * E_S[2]) / sr * sz - sr * E_S[3], (nsy / n * E_S[1] - nsx / n * E_S[2]) / sr)
	end
end

function polarizedcomponents(angspe::AbstractFieldAngularSpectrum)::Tuple{Array{Complex{Float64}, 2}, Array{Complex{Float64}, 2}}
	sizeX = length(angspe.nsx_X);
	sizeY = length(angspe.nsy_Y);
	Eₚ_XY = Array{Complex{Float64}, 2}(undef, sizeX, sizeY)
	Eₛ_XY = Array{Complex{Float64}, 2}(undef, sizeX, sizeY)
	size(angspe.e_SXY, 1) == 1 && error("No polarized components can be calculated for scallar approximation");
	@inbounds @simd for iY in 1:sizeY
		for iX in 1:sizeX
			tmpe = @view angspe.e_SXY[:,iX,iY];
			(Eₚ_XY[iX,iY], Eₛ_XY[iX,iY]) = polarizedcomponents(tmpe, angspe.nsx_X[iX], angspe.nsy_Y[iY], angspe.dir, angspe.n);
		end
	end
	return (Eₚ_XY, Eₛ_XY)
end

@inline function polarizeddirections!(nₚ::AbstractVector{T}, nₛ::AbstractVector{T}, nsx::Real, nsy::Real, dir::Integer, n::Number) where {T<:Number}
	sr = √(nsx^2 + nsy^2) / n
	if abs2(sr) < 1E-20
		@show "this has to be chekect(polarizeddirections! of Angularspectrum)"
		nₚ[1], nₚ[2], nₚ[3] = dir * one(T), zero(T), zero(T)
		nₛ[1], nₛ[2], nₛ[3] = zero(T), -one(T), zero(T)
	else
		sz = √(complex(1 - sr^2));
		nₚ[1], nₚ[2], nₚ[3] = dir * nsx / n / sr * sz, dir * nsy / n / sr * sz, - sr
		nₛ[1], nₛ[2], nₛ[3] = nsy / n / sr, -nsx / n / sr, zero(T);
	end
end

function polarizeddirections(nsx_X::AbstractVector{<:Number}, nsy_Y::AbstractVector{<:Number}, dir::Integer)::Tuple{Array{Complex{Float64}, 3}, Array{Complex{Float64}, 3}}
	sizeX = length(nsx_X);
	sizeY = length(nsy_Y);
	nₛ_3XY = Array{Complex{Float64}, 3}(undef, 3, sizeX, sizeY)
	nₚ_3XY = Array{Complex{Float64}, 3}(undef, 3, sizeX, sizeY)

	@inbounds @simd for iY in 1:sizeY
		for iX in 1:sizeX
			auxnₚ = @view nₚ_3XY[:,iX,iY];
			auxnₛ = @view nₛ_3XY[:,iX,iY];
			polarizeddirections!(auxnₚ, auxnₛ, nsx_X[iX], nsy_Y[iY], dir, 1);
		end
	end
	return (nₚ_3XY, nₛ_3XY);
end

function polarizeddirections(angspe::AbstractFieldAngularSpectrum)::Tuple{Array{Complex{Float64}, 3}, Array{Complex{Float64}, 3}}
	size(angspe.e_SXY, 1) == 1 && error("No polarized directions can be calculated for scalar approximation");
	return polarizeddirections(angspe.nsx_X, angspe.nsy_Y, angspe.dir);
end

function scalartovectorial!(angspe::FieldAngularSpectrum; pol = :x )::FieldAngularSpectrum
	if (size(angspe.e_SXY, 1) == 1)
		oldE_SXY = copy(angspe.e_SXY);
		angspe.e_SXY = Array{Complex{Float64}, 3}(undef, 3, length(angspe.nsx_X), length(angspe.nsy_Y));

		if pol == :x
			ePol = [1; 0];
		elseif pol == :y
			ePol = [0; 1];
		elseif pol <: AbstractVector{<:Number}
			length(pol) != 2 && error("Polarization type not known. Use :x, :y or a vector specifiying the polarization as [.8; im]");
			ePol = pol ./ norm(pol);
		else
			error("Polarization type not known. Use :x, :y or a vector specifiying the polarization as [.8; im]");
		end
		nₚ = MVector{3, ComplexF64}(undef)
		nₛ = MVector{3, ComplexF64}(undef)
		for iX in eachindex(angspe.nsy_Y)
			for iY in eachindex(angspe.nsx_X)
				polarizeddirections!(nₚ, nₛ, angspe.nsx_X[iX], angspe.nsy_Y[iY], angspe.dir, angspe.n);
				nsr = √(angspe.nsx_X[iX]^2 + angspe.nsy_Y[iY]^2);
				Eₚ = oldE_SXY[1,iX,iY] * (ePol[1] * angspe.dir * angspe.nsx_X[iX] + ePol[2] * angspe.nsy_Y[iY] * angspe.dir) / nsr;
				Eₛ = oldE_SXY[1,iX,iY] * (ePol[1] * angspe.nsy_Y[iY] + ePol[2] * -angspe.nsx_X[iX]) / nsr;
				angspe.e_SXY[1,iX,iY] = nₚ[1] * Eₚ + nₛ[1] * Eₛ
				angspe.e_SXY[2,iX,iY] = nₚ[2] * Eₚ + nₛ[2] * Eₛ
				angspe.e_SXY[3,iX,iY] = nₚ[3] * Eₚ
			end
		end
	else
		println("Field is already vectorial. The function scallartovectorial did not apply any change to FieldAngularSpectrum");
	end
	return angspe;
end
scalartovectorial(angspe::FieldAngularSpectrum) = scalartovectorial!(deepcopy(angspe));
function vectorialtoscalar!(angspe::FieldAngularSpectrum)::FieldAngularSpectrum
	if (size(angspe.e_SXY, 1) == 3)
		(Eₚ_XY, Eₛ_XY) = polarizedcomponents(angspe);
		angspe.e_SXY = Eₚ_XY + Eₛ_XY;
		println("The function vectorial to scallar does not conserve energy as well as losing information on the phase of the plane waves. Avoid using the function when possible.");
	else
		println("Field is already vectorial. The function vectorialtoscallar did not apply any change to FieldAngularSpectrum.");
	end
	return angspe;
end
vectorialtoscalar(angspe::FieldAngularSpectrum) = vectorialtoscalar!(deepcopy(angspe));

function intensity(angspe::FieldAngularSpectrum)::Float64
	k = 2 * π / angspe.λ;
	kx_X = k .* angspe.nsx_X;
	ky_Y = k .* angspe.nsy_Y;

	size(angspe.e_SXY, 1) == 3 ? e_SXY = dotdim(angspe.e_SXY, angspe.e_SXY, 1) : e_SXY = abs2.(angspe.e_SXY[1,:,:]);
	return 4 * π^2 * angspe.n * ∫∫(e_SXY, kx_X, ky_Y);
end

function centerangspereferenceframe!(angspe::FieldAngularSpectrum)
	(tmp, arg) = findmax(abs.(vec(angspe.e_SXY)))
	arg = CartesianIndices(angspe.e_SXY)[arg]
	nsx = angspe.nsx_X[arg[2]];
	nsy = angspe.nsy_Y[arg[3]];
	nsz = √(angspe.n^2 - nsx^2 - nsy^2)
	(abs(imag(nsz)) < @tol) || error("Rotation can only be performed in non dispersive media")

	(Δx, Δy, Δz) = rotatecoordinatesto(nsx, nsy, nsz, angspe.ref.θ, angspe.ref.ϕ)
	θ = acos(real(Δz))
	ϕ = atan(real(Δy), real(Δx))
	# θ > π/2 && 	(θ = π - θ; ϕ = ϕ + π)
	ref = ReferenceFrame(angspe.ref.x, angspe.ref.y, angspe.ref.z, θ, ϕ)
	changereferenceframe!(angspe, ref)
end

function samedefinitions(fieldl::L, fieldr::L) where L <: FieldAngularSpectrum
	isapprox(fieldl.nsx_X, fieldr.nsx_X, atol = @tol) || return false
	isapprox(fieldl.nsy_Y, fieldr.nsx_X, atol = @tol) || return false
	isapprox(fieldl.n, fieldr.n, atol = @tol) || return false
	isapprox(fieldl.λ, fieldr.λ, atol = @tol) || return false
	checkorientation(fieldl.ref, fieldr.ref) || return false
	return true
end

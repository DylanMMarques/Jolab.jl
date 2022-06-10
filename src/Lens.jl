mutable struct Lens{T} <: AbstractOpticalComponent{T}
	f::JolabFunction{T}
	na::JolabFunction{T}
	ref::ReferenceFrame{T}
end

"""
	Lens(T, f, na, ref)

Initializes a multilayer structure.
- `T` number type specifing data precision. Ex: `Float32`, `BigFloat`, ... The default is Float64;
- `f` - defines the lens focal length in meters. A function or interpolation object can be used to describe a focal length wavelength dependent.
- `na` - defines the lens numerical aperture. A function or interpolation object can be used to describe a focal length wavelength dependent.
- `ref` - `ReferenceFrame` specifying the lens position and orientation. The focal planes are located at the plus and minus the focal length.

**Examples:**

```julia
lens = Lens(10E-3, 0.5, ReferenceFrame(0,0,0,0,0))
```

```julia
f(λ) = 10E-3 * λ
na(λ) = .1 * λ
lens = Lens(f, na, ReferenceFrame(0,0,0,0,0))
```
"""
Lens(::Type{T}, f, na, ref) where T = Lens{T}(f, na, ref)
Lens(f, na, ref) = Lens{Float64}(f, na, ref)
f(lens::Lens, λ) = real(lens.f(λ))
na(lens::Lens, λ) = real(lens.na(λ))

function ref1(lens::Lens{T}, λ::Real) where T
	x = lens.ref.x - f(lens, λ) * sin(lens.ref.θ) * cos(lens.ref.ϕ)
	y = lens.ref.y - f(lens, λ) * sin(lens.ref.θ) * sin(lens.ref.ϕ)
	z = lens.ref.z - f(lens, λ) * cos(lens.ref.θ)
	return ReferenceFrame{T}(x, y, z, lens.ref.θ, lens.ref.ϕ)
end

function ref2(lens::Lens{T}, λ::Real) where T
	x = lens.ref.x + f(lens, λ) * sin(lens.ref.θ) * cos(lens.ref.ϕ)
	y = lens.ref.y + f(lens, λ) * sin(lens.ref.θ) * sin(lens.ref.ϕ)
	z = lens.ref.z + f(lens, λ) * cos(lens.ref.θ)
	return ReferenceFrame{T}(x, y, z, lens.ref.θ, lens.ref.ϕ)
end

function t(lens::Lens, field::FieldAngularSpectrumScalarRadialSymmetric{T}, iR, iY) where T
	cosθ² = 1 - (field.nsr_R[iR] / real(field.n))^2
	k = 2π / field.λ
	if cosθ² < 0 # Evasnecent waves
		return (zero(Complex{T}), zero(Complex{T}))
	else # Plane waves
		if (1 - cosθ² > na(lens, field.λ)^2) # above the lens NA
			return (zero(Complex{T}), zero(Complex{T}))
		else
			aux = f(lens, field.λ) / im / k / 2π * cosθ²^(1/4) # might be wrong
			return (1 / aux, aux)
		end
	end
end

function t(lens::Lens, field::FieldAngularSpectrumScalar{T}, iX, iY) where T
	cosθ² = 1 - (field.nsx_X[iX] / real(field.n))^2 - (field.nsy_Y[iY] / real(field.n))^2
	k = 2π / field.λ
	if cosθ² < 0 # Evasnecent waves
		return (zero(Complex{T}), zero(Complex{T}))
	else # Plane waves
		if (1 - cosθ² > na(lens, field.λ)^2) # above the lens NA
			return (zero(Complex{T}), zero(Complex{T}))
		else
			aux = f(lens, field.λ) / im / k / 2π * cosθ²^(1/4) # might be wrong
			return (1 / aux, aux)
		end
	end
end

function t(lens::Lens, field::FieldSpaceScalar{T}, iX, iY) where T
	cosθ² = 1 - (field.x_X[iX] / f(lens, field.λ))^2 - (field.y_Y[iY] / f(lens, field.λ))^2
	k = 2π / field.λ
	if cosθ² < 0 # Evasnecent waves
		return (zero(Complex{T}), zero(Complex{T}))
	else # Plane waves
		if (1 - cosθ² > na(lens,field.λ)^2) # above the lens NA
			return (zero(Complex{T}), zero(Complex{T}))
		else
			aux = f(lens, field.λ) / im / k / 2π * cosθ²^(1/4) # might be wrong
			return (aux, 1 / aux)
		end
	end
end

function t(lens::Lens, field::FieldSpaceScalarRadialSymmetric{T}, iR, iY) where T
	cosθ² = 1 - (field.r_R[iR] / f(lens, field.λ))^2
	k = 2π / field.λ
	if cosθ² < 0 # Evasnecent waves
		return (zero(Complex{T}), zero(Complex{T}))
	else # Plane waves
		if (1 - cosθ² > na(lens,field.λ)^2) # above the lens NA
			return (zero(Complex{T}), zero(Complex{T}))
		else
			aux = f(lens, field.λ) / im / k / 2π * cosθ²^(1/4) # might be wrong
			return (aux, 1 / aux)
		end
	end
end

function checkapplicability(lens::Lens, field::AbstractFieldAngularSpectrum)
	(abs(imag(field.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
	checkorientation(field.ref, lens.ref) || error("need to be done")
	return true
end

function checkapplicability(lens::Lens, field::AbstractFieldSpace)
	(abs(imag(field.n)) < @tol) || error("To apply a lens the medium cannot absorb light.")
	#TO DO: update this to allow in plane. And making the correction on t(lens, field space for the reference frame change)
	field.ref == (dir(field) > 0 ? ref1(lens, field.λ) : ref2(lens, field.λ)) || error("The reference frame of the field must be a focal plane of the lens.")
	return true
end

function getfields_lr(lens::Lens, fieldi::FieldAngularSpectrumScalar{T,D,X,B}) where {T,D,X,B}
	f_val = f(lens, fieldi.λ)
	x_X = fieldi.nsx_X / real(fieldi.n) * f_val
	y_Y = fieldi.nsy_Y / real(fieldi.n) * f_val

	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref1(lens, fieldi.λ))
		fieldr = FieldSpaceScalar{T,1,X,B}(x_X, y_Y, 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref2(lens, fieldi.λ))
	else
		fieldl = FieldSpaceScalar{T,-1,X,B}(x_X, y_Y, 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref1(lens, fieldi.λ))
		fieldr = FieldAngularSpectrumScalar{T,1,X,B}(deepcopy(fieldi.nsx_X), deepcopy(fieldi.nsy_Y), 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref2(lens, fieldi.λ))
	end
	return (fieldl, fieldr)
end

function getfields_lr(lens::Lens, fieldi::FieldSpaceScalar{T,D,X,B}) where {T,D,X,B}
	f_val = f(lens, fieldi.λ)
	nsx_X = -fieldi.x_X / f_val * real(fieldi.n)
	nsy_Y = -fieldi.y_Y / f_val * real(fieldi.n)

	if dir(fieldi) > 0
		fieldl = FieldSpaceScalar{T,-1,X,B}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref1(lens, fieldi.λ))
		fieldr = FieldAngularSpectrumScalar{T,1,X,B}(nsx_X, nsy_Y, 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref2(lens, fieldi.λ))
	else
		fieldl = FieldAngularSpectrumScalar{T,-1,X,B}(nsx_X, nsy_Y, 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref1(lens, fieldi.λ))
		fieldr = FieldSpaceScalar{T,1,X,B}(deepcopy(fieldi.x_X), deepcopy(fieldi.y_Y), 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref2(lens, fieldi.λ))
	end
	return (fieldl, fieldr)
end

function getfields_lr(lens::Lens, fieldi::FieldSpaceScalarRadialSymmetric{T,D,X,B}) where {T,D,X,B}
	f_val = f(lens, fieldi.λ)
	nsr_R = -fieldi.r_R / f_val * real(fieldi.n)

	if dir(fieldi) > 0
		fieldl = FieldSpaceScalarRadialSymmetric{T,-1,X,B}(deepcopy(fieldi.r_R), 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref1(lens, fieldi.λ))
		fieldr = FieldAngularSpectrumScalarRadialSymmetric{T,1,X,B}(nsr_R, 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref2(lens, fieldi.λ))
	else
		fieldl = FieldAngularSpectrumScalarRadialSymmetric{T,-1,X,B}(nsr_R, 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref1(lens, fieldi.λ))
		fieldr = FieldSpaceScalarRadialSymmetric{T,1,X,B}(deepcopy(fieldi.r_R), 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref2(lens, fieldi.λ))
	end
	return (fieldl, fieldr)
end

function getfields_lr(lens::Lens, fieldi::FieldAngularSpectrumScalarRadialSymmetric{T,D,X,B}) where {T,D,X,B}
	f_val = f(lens, fieldi.λ)
	r_R = fieldi.nsr_R / real(fieldi.n) * f_val

	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrumScalarRadialSymmetric{T,-1,X,B}(deepcopy(fieldi.nsr_R), 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref1(lens, fieldi.λ))
		fieldr = FieldSpaceScalarRadialSymmetric{T,1,X,B}(r_R, 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref2(lens, fieldi.λ))
	else
		fieldl = FieldSpaceScalarRadialSymmetric{T,-1,X,B}(r_R, 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref1(lens, fieldi.λ))
		fieldr = FieldAngularSpectrumScalarRadialSymmetric{T,1,X,B}(deepcopy(fieldi.nsr_R), 0deepcopy(fieldi.e_SXY), fieldi.λ, fieldi.n, ref2(lens, fieldi.λ))
	end
	return (fieldl, fieldr)
end

function get_scatteringmatrixtype(lens::Lens, field::FieldSpaceScalar{T,D,Y,B}) where {T,D,Y,B}
	(fieldl, fieldr) = getfields_lr(lens, field)
	r12 = nothing
	r21 = nothing
	t12 = Diagonal(Vector{Complex{T}}(undef, length(field.e_SXY)))
	t21 = Diagonal(Vector{Complex{T}}(undef, length(field.e_SXY)))
	if dir(field) > 0
    	return ScatteringMatrix{T,FieldSpaceScalar{T,-1,Y,B}, FieldAngularSpectrumScalar{T,1,Y,B}, Nothing, Diagonal{Complex{T},Vector{Complex{T}}}}(r12,t12,r21,t21,fieldl,fieldr)
	else
    	return ScatteringMatrix{T,FieldAngularSpectrumScalar{T,-1,Y,B}, FieldSpaceScalar{T,1,Y,B}, Nothing, Diagonal{Complex{T},Vector{Complex{T}}}}(r12,t12,r21,t21,fieldl,fieldr)
	end
end

function get_scatteringmatrixtype(lens::Lens, field::FieldAngularSpectrumScalar{T,D,Y,B}) where {T,D,Y,B}
	(fieldl, fieldr) = getfields_lr(lens, field)
	r12 = nothing
	r21 = nothing
	t12 = Diagonal(Vector{Complex{T}}(undef, length(field.e_SXY)))
	t21 = Diagonal(Vector{Complex{T}}(undef, length(field.e_SXY)))
	if dir(field) > 0
    	return ScatteringMatrix{T,FieldAngularSpectrumScalar{T,-1,Y,B}, FieldSpaceScalar{T,1,Y,B}, Nothing, Diagonal{Complex{T},Vector{Complex{T}}}}(r12,t12,r21,t21,fieldl,fieldr)
	else
    	return ScatteringMatrix{T,FieldSpaceScalar{T,-1,Y,B}, FieldAngularSpectrumScalar{T,1,Y,B}, Nothing, Diagonal{Complex{T},Vector{Complex{T}}}}(r12,t12,r21,t21,fieldl,fieldr)
	end
end

function get_scatteringmatrixtype(lens::Lens, field::FieldAngularSpectrumScalarRadialSymmetric{T,D,Y,B}) where {T,D,Y,B}
	(fieldl, fieldr) = getfields_lr(lens, field)
	r12 = nothing
	r21 = nothing
	t12 = Diagonal(Vector{Complex{T}}(undef, length(field.e_SXY)))
	t21 = Diagonal(Vector{Complex{T}}(undef, length(field.e_SXY)))
	if dir(field) > 0
    	return ScatteringMatrix{T,FieldAngularSpectrumScalarRadialSymmetric{T,-1,Y,B}, FieldSpaceScalarRadialSymmetric{T,1,Y,B}, Nothing, Diagonal{Complex{T},Vector{Complex{T}}}}(r12,t12,r21,t21,fieldl,fieldr)
	else
    	return ScatteringMatrix{T,FieldSpaceScalarRadialSymmetric{T,-1,Y,B}, FieldAngularSpectrumScalarRadialSymmetric{T,1,Y,B}, Nothing, Diagonal{Complex{T},Vector{Complex{T}}}}(r12,t12,r21,t21,fieldl,fieldr)
	end
end

function get_scatteringmatrixtype(lens::Lens, field::FieldSpaceScalarRadialSymmetric{T,D,Y,B}) where {T,D,Y,B}
	(fieldl, fieldr) = getfields_lr(lens, field)
	r12 = nothing
	r21 = nothing
	t12 = Diagonal(Vector{Complex{T}}(undef, length(field.e_SXY)))
	t21 = Diagonal(Vector{Complex{T}}(undef, length(field.e_SXY)))
	if dir(field) > 0
    	return ScatteringMatrix{T,FieldSpaceScalarRadialSymmetric{T,-1,Y,B}, FieldAngularSpectrumScalarRadialSymmetric{T,1,Y,B}, Nothing, Diagonal{Complex{T},Vector{Complex{T}}}}(r12,t12,r21,t21,fieldl,fieldr)
	else
    	return ScatteringMatrix{T,FieldAngularSpectrumScalarRadialSymmetric{T,-1,Y,B}, FieldSpaceScalarRadialSymmetric{T,1,Y,B}, Nothing, Diagonal{Complex{T},Vector{Complex{T}}}}(r12,t12,r21,t21,fieldl,fieldr)
	end
end

"""
	(field_back, field_forward)= coefficient_general(lens, angspe)

Calculates the field after propagating through the lens assuming the incident field direction of propagation. The model propagates from a focal plane to the other focal plane. The order is based on the light direction of propagation.
- **Type:** Transmission matrices is diagonal. No reflection matrix is computed (the model does not assume light reflection by the lens)
- **Time:** very short; scales with `length(angspe.nsx_X) length(angspe.nsy_Y)`
- **RAM:** very small; scales with `length(angspe.nsx_X)` `length(angspe.nsy_Y)`
- **Convergence:** sampling of `angspe.nsx_X` and `angspe.nsy_Y`
"""
function coefficient_general(lens::Lens, fieldi::Union{AbstractFieldAngularSpectrum, AbstractFieldSpace})
	checkapplicability(lens, fieldi) || errorToDo()

	scat = get_scatteringmatrixtype(lens, fieldi)
	cart = CartesianIndices(fieldi)
	@inbounds for i in iterator_index(fieldi)
		if dir(fieldi) > 0
			(scat.t₁₂.diag[i], scat.t₂₁.diag[i]) = t(lens, fieldi, cart[i][2], cart[i][3])
		else
			(scat.t₂₁.diag[i], scat.t₁₂.diag[i]) = t(lens, fieldi, cart[i][2], cart[i][3])
		end
	end
	correctscatteringmatrix_referenceframes!(scat, lens, fieldi)
	return scat
end

"""
	(field_back, field_forward) = lightinteraction(lens, angspe)

Calculates the field after propagating through the lens assuming the incident field direction of propagation. The model propagates from a focal plane to the other focal plane. The order is based on the light direction of propagation.
- **Type:** Transmission matrices is diagonal. No reflection matrix is computed (the model does not assume light reflection by the lens)
- **Time:** very short; scales with `length(angspe.nsx_X) length(angspe.nsy_Y)`
- **RAM:** very small; scales with `length(angspe.nsx_X)` `length(angspe.nsy_Y)`
- **Convergence:** sampling of `angspe.nsx_X` and `angspe.nsy_Y`
"""
# function lightinteractionvectorial(lens::Lens, angspe::FieldAngularSpectrum)
# 	@show("This must be remade. Do not trust result results")
# 	e_SXY = Array{Complex{Float64}, 3}(undef, size(angspe.e_SXY));
# 	(abs(imag(angspe.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
# 	sx_X = angspe.nsx_X;
# 	sy_Y = reshape(angspe.nsy_Y, 1,:)
# 	θ = acos.(.√(complex.(1 .- sx_X.^2 .- sy_Y.^2)))
# 	ϕ = atan.(sx_X, sy_Y)
#
# 	cosθ = cos.(θ); sinθ = sin.(θ); cosϕ = cos.(ϕ); sinϕ = sin.(ϕ);
#
# 	e_SXY[1,:,:] .= cosϕ .* (angspe.e_SXY[1,:,:] .* cosϕ .* cosθ .- angspe.e_SXY[3,:,:] .* sinθ) .+ angspe.e_SXY[2,:,:] .* (cosθ .- 1) .* cosϕ .* sinϕ .+ angspe.e_SXY[1,:,:] .* sinϕ.^2;
# 	e_SXY[2,:,:] .= angspe.e_SXY[2,:,:] .* cosϕ.^2 .+ angspe.e_SXY[1,:,:] .* (cosθ .- 1) .* cosϕ .* sinϕ .+ sinϕ .* (-angspe.e_SXY[3,:,:] .* sinθ .+ angspe.e_SXY[2,:,:] .* cosθ .* sinϕ);
# 	e_SXY[3,:,:] .= angspe.e_SXY[3,:,:] .* cosθ .+ sinθ .* (angspe.e_SXY[1,:,:] .* cosϕ .+ angspe.e_SXY[2,:,:] .* sinϕ);
#
# 	k = 2 .* π .* angspe.n ./ angspe.λ;
# 	e_SXY .= im ./ lens.f(angspe.λ) * k * 2 * π ./ .√(adddims(cosθ, (1,))) .* e_SXY;
#
# 	e_SXY .= e_SXY .* (adddims(real.(sinθ), (1,)) .< lens.na);
#
# 	x_X = sx_X * lens.f(angspe.λ);
# 	y_Y = angspe.nsy_Y * lens.f(angspe.λ);
#
# 	return (x_X, y_Y, e_SXY)
# end
#
# function lightinteractionscalar(lens::Lens, angspe::AbstractFieldAngularSpectrum{T}) where T
# 	(abs(imag(angspe.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
#
# 	sx_X = angspe.nsx_X / real(angspe.n)
# 	sy_Y = angspe.nsy_Y / real(angspe.n)
#
# 	x_X = sx_X * lens.f(angspe.λ)
# 	y_Y = sy_Y * lens.f(angspe.λ)
#
# 	k = 2 * π * angspe.n / angspe.λ
# 	e_SXY = Array{Complex{T}, 3}(undef, 1, length(sx_X), length(sy_Y))
# 	@inbounds Threads.@threads for iY in eachindex(sy_Y)
# 		for iX in eachindex(sx_X)
# 			cosθ2 = 1 - sx_X[iX]^2 - sy_Y[iY]^2 # cosθ2 is cosθ^2
# 			if cosθ2 < 0 # Evasnecent waves
# 				e_SXY[1,iX,iY] = zero(Complex{T})
# 			else # Plane waves
# 				(1 - cosθ2 > lens.na^2) && (e_SXY[1,iX,iY] = zero(Complex{T}); continue)
# 				e_SXY[1,iX,iY] = im / lens.f(angspe.λ) * k * 2π / cosθ2^(1/4) * angspe.e_SXY[1,iX,iY]
# 			end
# 		end
# 	end
# 	return (x_X, y_Y, e_SXY)
# end
#
# function lightinteractionscalar(lens::Lens, space::AbstractFieldSpace{T}) where T
# 	(abs(imag(space.n)) < @tol) || error("To apply a lens the medium cannot absorb light")
# 	f = lens.f(space.λ)
# 	sx_X = -space.x_X / f
# 	sy_Y = -space.y_Y / f
#
# 	k = 2 * π * space.n / space.λ
# 	e_SXY = Array{Complex{T},3}(undef, 1, length(space.x_X), length(space.y_Y))
# 	@inbounds Threads.@threads for iY in eachindex(space.y_Y)
# 		for iX in eachindex(space.x_X)
# 			cosθ2 = (f^2 - space.x_X[iX]^2 - space.y_Y[iY]^2)
# 			if cosθ2 < 0
# 				e_SXY[1,iX,iY] = zero(Complex{T})
# 			else
# 				cosθ2 = √(cosθ2) / f
# 				(1 - cosθ2 > lens.na^2) && (e_SXY[1,iX,iY] = zero(Complex{T}); continue)
# 				e_SXY[1,iX,iY] = -im * f / k / 2 / π * √(cosθ2) * space.e_SXY[1,iX,iY];
# 			end
# 		end
# 	end
# 	sx_X *= real(space.n) # is nsx
# 	sy_Y *= real(space.n) # is nsy
# 	return (sx_X, sy_Y, e_SXY)
# end
#
# function lightinteractionvectorial(lens::Lens, space::FieldSpace)
# 	@show "This must be remade. Do not trust results"
# 	e_SXY = Array{Complex{Float64}, 3}(undef, size(space.e_SXY));
# 	x_X = space.x_X;
# 	y_Y = reshape(space.y_Y, 1, :);
# 	θ = acos.(.√(complex.(lens.f(space.λ).^2 .- x_X.^2 .- y_Y.^2)) ./ lens.f(space.λ));
# 	ϕ = π .+ atan.(y_Y, x_X);
#
# 	cosθ = cos.(θ); sinθ = sin.(θ); cosϕ = cos.(ϕ); sinϕ = sin.(ϕ);
#
# 	e_SXY[1,:,:] .= (cosθ .* cosϕ .* cosϕ .+ sinϕ.^2) .* space.e_SXY[1,:,:] .+ cosϕ .* sinθ .* space.e_SXY[3,:,:] .+ (cosθ .- 1) .* cosϕ .* sinϕ .* space.e_SXY[2,:,:];
# 	e_SXY[2,:,:] .= (cosϕ.^2 .+ cosθ .* sinϕ.^2) .* space.e_SXY[2,:,:] + (cosθ .- 1) .* cosϕ .* sinϕ .* space.e_SXY[1,:,:] .+ sinθ .* sinϕ .* space.e_SXY[3,:,:];
# 	e_SXY[3,:,:] .= .- sinθ .* cosϕ .* space.e_SXY[1,:,:] .- sinθ .* sinϕ .* space.e_SXY[2,:,:] .+ cosθ .* space.e_SXY[3,:,:];
#
# 	k = 2 .* π .* space.n ./ space.λ;
# 	e_SXY .= -im .* lens.f(space.λ) ./ k ./ 2 ./ π .* .√(adddims(cosθ, (1,))) .* e_SXY;
# 	e_SXY .= e_SXY .* (adddims(real.(sinθ), (1,)) .< lens.na);
#
# 	sx_X = .-space.x_X ./ lens.f(space.λ);
# 	sy_Y = .-space.y_Y ./ lens.f(space.λ);
#
# 	return (sx_X, sy_Y, e_SXY)
# end

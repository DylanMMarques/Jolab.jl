abstract type AbstractPropagationCoefficient{T,X<:Union{JolabFunction1D{T,Complex{T}}, AbstractArray{Complex{T}}}} <: AbstractCoefficient{T} end

struct PropagationCoefficientScalar{T,X} <: AbstractPropagationCoefficient{T,X}
	r₁₂::X
	t₁₂::X
	r₂₁::X
	t₂₁::X
	λ::T
	n₁::Complex{T}
	ref₁::ReferenceFrame{T}
	n₂::Complex{T}
	ref₂::ReferenceFrame{T}
	PropagationCoefficientScalar{T,X}(r₁₂, t₁₂, r₂₁, t₂₁, λ, n₁, ref₁, n₂, ref₂) where {T,X} = new{T,X}(r₁₂, t₁₂, r₂₁, t₂₁, λ, n₁, ref₁, n₂, ref₂);
	PropagationCoefficientScalar{T}(r₁₂, t₁₂, r₂₁, t₂₁, λ, n₁, ref₁, n₂, ref₂) where T = new{T,JolabFunction1D{T,Complex{T}}}(r₁₂, t₁₂, r₂₁, t₂₁, λ, n₁, ref₁, n₂, ref₂);
end
# PropagationCoefficientScalar(r₁₂, t₁₂, r₂₁, t₂₁, λ, n₁, ref₁, n₂, ref₂) = PropagationCoefficientScalar{Float64}(r₁₂, t₁₂, r₂₁, t₂₁, λ, n₁, ref₁, n₂, ref₂)

struct PropagationCoefficientVectorial{T,X} <: AbstractPropagationCoefficient{T,X}
	rpp₁₂::X
	rss₁₂::X
	tpp₁₂::X
	tss₁₂::X
	rpp₂₁::X
	rss₂₁::X
	tpp₂₁::X
	tss₂₁::X
	λ::T
	n₁::Complex{T}
	ref₁::ReferenceFrame{T}
	n₂::Complex{T}
	ref₂::ReferenceFrame{T}
	function PropagationCoefficientVectorial{T,X}(rpp₁₂, rss₁₂, tpp₁₂, tss₁₂, rpp₂₁, rss₂₁, tpp₂₁, tss₂₁, λ, n₁, ref₁, n₂, ref₂) where {T,X}
		return new{T,X}(rpp₁₂, rss₁₂, tpp₁₂, tss₁₂, rpp₂₁, rss₂₁, tpp₂₁, tss₂₁, λ, n₁, ref₁, n₂, ref₂)
	end
end

function checkmergeble(coefs::AbstractArray{<:AbstractPropagationCoefficient})
	@inbounds for i in 1:length(coefs)-1
		(coefs[i].λ - coefs[i+1].λ < @tol) || return false;
		checkorientation(coefs[i].ref₁, coefs[i+1].ref₁) || return false
		(coefs[i+1].ref₁.z - coefs[i].ref₂.z > -@tol) || return false
		checkinline(coefs[i].ref₁, coefs[i+1].ref₁) || return false
		(abs(coefs[i+1].n₁ - coefs[i].n₂) < @tol) || return false
	end
	return true
end

function checkmergeble(coefs::AbstractArray{<:AbstractPropagationCoefficient{T, A}}) where {T <: Real, A <: AbstractArray{T}}
	@inbounds for i in 1:length(coefs)-1
		(coefs[i].λ - coefs[i+1].λ < @tol) || return false;
		checkorientation(coefs[i].ref₁, coefs[i+1].ref₁) || return false
		(coefs[i+1].ref₁.z - coefs[i].ref₂.z > -@tol) || return false
		checkinline(coefs[i].ref₁, coefs[i+1].ref₁) || return false
		(abs(coefs[i+1].n₁ - coefs[i].n₂) < @tol) || return false
		size(coefs[i]) == size(coefs[i+1]) || return false
	end
	return true
end

function lightinteraction(coefs::AbstractVector{<:AbstractPropagationCoefficient}, angspe::AbstractFieldAngularSpectrum{T}) where T
	checkapplicability(coefs, angspe);
	(sizeS, sizeX, sizeY) = size(angspe.e_SXY);

	nsx_XY = repeat(angspe.nsx_X, 1, sizeY);
	nsy_XY = repeat(angspe.nsy_Y', sizeX, 1);
	nsz_XY = angspe.dir * .√(angspe.n^2 .- nsx_XY.^2 .- nsy_XY.^2)

	refcoef = angspe.dir > 0 ? coefs[1].ref₁ : coefs[2].ref₂
	changereferential!(nsx_XY, nsy_XY, nsz_XY, ei_SXY, angspe.λ, angspe.ref, refcoef)
	checkorientation(angspe.ref, refcoef) || error("Cannot do that yet. Resampling of nsx must be done: TODO")

	coef = mergeorientated_propagationcoefficient(coefs);
	(nsrmin, nsrmax) = rextrema(angspe.nsx_X, angspe.nsy_Y);
	nsr = range(nsrmin, stop = nsrmax, length = round(Int, max(sizeX, sizeY) / √2))
	coef_itp = coefficient_itp(coef, nsr)

	lightinteraction!(ebackward_SCD, eforward_SCD, coef_itp, nsx_XY, nsy_XY, ei_SXY, angspe.λ, angspe.dir);

	angspeforward = FieldAngularSpectrum{T}(angspe.nsx_X, angspe.nsy_Y, eforward_SCD, angspe.λ, coefs[end].n₂, 1, coefs[end].ref₂)
	angspebackward = FieldAngularSpectrum{T}(angspe.nsx_X, angspe.nsy_Y, ebackward_SCD, angspe.λ, coefs[1].n₁, -1, coefs[1].ref₁)
	return (angspebackward, angspeforward)
end

function lightinteraction!(ebackward_SXY::AbstractArray{<:Number, 3}, eforward_SXY::AbstractArray{<:Number, 3}, coef::PropagationCoefficientScalar, nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number, 2}, ei_SXY::AbstractArray{<:Number, 3}, dir::Integer, tmpnₛ=[1.], tmpnₚ=[1.])
	length(nsx_XY) == length(nsy_XY) == length(ei_SXY) == length(ebackward_SXY) == length(eforward_SXY)  || error("Wrong sizes")
	@inbounds @simd for i in eachindex(nsx_XY)
		nsr = real(√(nsx_XY[i]^2 + nsy_XY[i]^2));
		if dir > 0
			ebackward_SXY[i] = ei_SXY[i] * coef.r₁₂(nsr);
			eforward_SXY[i] = ei_SXY[i] * coef.t₁₂(nsr);
		else
			ebackward_SXY[i] = ei_SXY[i] * coef.t₂₁(nsr);
			eforward_SXY[i] = ei_SXY[i] * coef.r₂₁(nsr);
		end
	end
end

function lightinteraction!(ebackward_SXY::AbstractArray{<:Number, 3}, eforward_SXY::AbstractArray{<:Number, 3}, coef::PropagationCoefficientScalar{T,X}, nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number, 2}, ei_SXY::AbstractArray{<:Number, 3}, dir::Integer, tmpnₛ=[1.], tmpnₚ=[1.]) where {T,X<:AbstractArray}
	if dir > 0
		vec(ebackward_SXY) .= vec(ei_SXY) .* vec(coef.r₁₂)
		vec(eforward_SXY) .= vec(ei_SXY) .* vec(coef.t₁₂)
	else
		vec(ebackward_SXY) .= vec(ei_SXY) .* vec(coef.t₂₁)
		vec(eforward_SXY) .= vec(ei_SXY) .* vec(coef.r₂₁)
	end
end

function lightinteraction!(ebackward_SXY::AbstractArray{Complex{T}, 3}, eforward_SXY::AbstractArray{Complex{T}, 3}, coef::PropagationCoefficientVectorial, nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number, 2}, ei_SXY::AbstractArray{<:Number, 3}, dir::Integer, tmpnₚ = MVector{3,Complex{T}}(undef)::AbstractVector{<:Number}, tmpnₛ = MVector{3,Complex{T}}(undef)::AbstractVector{<:Number}) where {T<:Real}
	n = dir > 0 ? coef.n₁ : coef.n₂
	tmpei_SXY = reshape(ei_SXY, size(ei_SXY, 1), :);
	length(nsx_XY) == length(nsy_XY) == size(tmpei_SXY, 2) || error("Wrong sizes")
	length(ebackward_SXY) == length(eforward_SXY) == length(ei_SXY) || error("Wrong sizes")
	length(tmpnₛ) == length(tmpnₚ) > 2 || error("Wrong sizes")

	tmpeforward_SXY = reshape(eforward_SXY, size(eforward_SXY, 1), :);
	tmpebackward_SXY = reshape(ebackward_SXY, size(ebackward_SXY, 1), :);

	@inbounds @simd for i in eachindex(nsx_XY)
		nsr = real(√(nsx_XY[i]^2 + nsy_XY[i]^2));
		if dir > 0
			(Ep, Es) = polarizedcomponents(tmpei_SXY[:,i], nsx_XY[i], nsy_XY[i], 1, coef.n₁);
			tmpEp = Ep * coef.rpp₁₂(nsr);
			tmpEs = Es * coef.rss₁₂(nsr);

			polarizeddirections!(tmpnₚ, tmpnₛ, nsx_XY[i], nsy_XY[i], -1, coef.n₁);
			tmpebackward_SXY[1,i] = tmpEp * tmpnₚ[1] + tmpEs * tmpnₛ[1];
			tmpebackward_SXY[2,i] = tmpEp * tmpnₚ[2] + tmpEs * tmpnₛ[2];
			tmpebackward_SXY[3,i] = tmpEp * tmpnₚ[3];

			tmpEp = Ep * coef.tpp₁₂(nsr);
			tmpEs = Es * coef.tss₁₂(nsr);

			polarizeddirections!(tmpnₚ, tmpnₛ, nsx_XY[i], nsy_XY[i], 1, coef.n₂);
			tmpeforward_SXY[1,i] = tmpEp * tmpnₚ[1] + tmpEs * tmpnₛ[1];
			tmpeforward_SXY[2,i] = tmpEp * tmpnₚ[2] + tmpEs * tmpnₛ[2];
			tmpeforward_SXY[3,i] = tmpEp * tmpnₚ[3];
		else
			(Ep, Es) = polarizedcomponents(tmpei_SXY[:,i], nsx_XY[i], nsy_XY[i], -1, coef.n₂);
			tmpEp = Ep * coef.rpp₂₁(nsr);
			tmpEs = Es * coef.rss₂₁(nsr);

			polarizeddirections!(tmpnₚ, tmpnₛ, nsx_XY[i], nsy_XY[i], 1, coef.n₂);
			tmpeforward_SXY[1,i] = tmpEp * tmpnₚ[1] + tmpEs * tmpnₛ[1];
			tmpeforward_SXY[2,i] = tmpEp * tmpnₚ[2] + tmpEs * tmpnₛ[2];
			tmpeforward_SXY[3,i] = tmpEp * tmpnₚ[3];

			tmpEp = Ep * coef.tpp₂₁(nsr);
			tmpEs = Es * coef.tss₂₁(nsr);

			polarizeddirections!(tmpnₚ, tmpnₛ, nsx_XY[i], nsy_XY[i], -1, coef.n₁);
			tmpebackward_SXY[1,i] = tmpEp * tmpnₚ[1] + tmpEs * tmpnₛ[1];
			tmpebackward_SXY[2,i] = tmpEp * tmpnₚ[2] + tmpEs * tmpnₛ[2];
			tmpebackward_SXY[3,i] = tmpEp * tmpnₚ[3];
		end
	end
end

function lightinteraction(coef::AbstractPropagationCoefficient, angspe::AbstractFieldAngularSpectrum{T}) where T
	checkapplicability(coef, angspe);
	(sizeS, sizeX, sizeY) = size(angspe.e_SXY);

	ebackward_SXY = Array{Complex{T}, 3}(undef, sizeS, sizeX, sizeY);
	eforward_SXY = Array{Complex{T}, 3}(undef, sizeS, sizeX, sizeY);
	ei_SXY = copy(angspe.e_SXY);

	nsx_XY = repeat(angspe.nsx_X, 1, sizeY);
	nsy_XY = repeat(angspe.nsy_Y', sizeX, 1);
	nsz_XY = angspe.dir .* .√(angspe.n^2 .- nsx_XY.^2 .- nsy_XY.^2)

	refcoef = angspe.dir > 0 ? coef.ref₁ : coef.ref₂
	changereferential!(nsx_XY, nsy_XY, nsz_XY, ei_SXY, angspe.λ, angspe.ref, refcoef)

	tmpnₛ = MVector{3,Complex{T}}(undef)
	tmpnₚ = MVector{3,Complex{T}}(undef)

	lightinteraction!(ebackward_SXY, eforward_SXY, coef, nsx_XY, nsy_XY, ei_SXY, angspe.dir, tmpnₚ, tmpnₛ)

	if checkorientation(angspe.ref, refcoef)
		nsx_X = nsx_XY[:,1]
		nsy_Y = nsy_XY[1,:]
	else
		(nsx_X, nsy_Y, eforward_SXY) = interpolatetogrid(nsx_XY, nsy_XY, eforward_SXY)
		ebackward_SXY = interpolatetogrid(nsx_XY, nsy_XY, ebackward_SXY, nsx_X, nsy_Y)
	end

	angspeforward = FieldAngularSpectrum{T}(nsx_X, nsy_Y, eforward_SXY, angspe.λ, coef.n₂, 1, coef.ref₂)
	angspebackward = FieldAngularSpectrum{T}(nsx_X, nsy_Y, ebackward_SXY, angspe.λ, coef.n₁, -1, coef.ref₁)
	return (angspebackward, angspeforward)
end

@inline function checkapplicability(coef::AbstractCoefficient, angspe::AbstractFieldAngularSpectrum)
	angspe.dir > 0 ? ((abs(angspe.n - coef.n₁) < @tol) || error("Field and coefficient refractive index are different") ) : ((abs(angspe.n - coef.n₂) < @tol) || error("Field and coefficient refractive index are different"))
	return true
end

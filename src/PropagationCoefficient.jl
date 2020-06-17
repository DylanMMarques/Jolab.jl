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

function coefficient_matrixform(coef::PropagationCoefficientScalar{T}, nsx_XY::AbstractArray{<:Real}, nsy_XY::AbstractArray{<:Real,2}) where {T}
	r₁₂_XY = Array{Complex{T}, 2}(undef, size(nsx_XY))
	t₁₂_XY = Array{Complex{T}, 2}(undef, size(nsx_XY))
	r₂₁_XY = Array{Complex{T}, 2}(undef, size(nsx_XY))
	t₂₁_XY = Array{Complex{T}, 2}(undef, size(nsx_XY))
	size(nsx_XY) == size(nsy_XY) || error("Wrong sizes")
	@inbounds Threads.@threads for i in eachindex(nsx_XY)
		nsr = √(nsx_XY[i]^2 + nsy_XY[i]^2)
		r₁₂_XY[i] = coef.r₁₂(nsr)
		t₁₂_XY[i] = coef.t₁₂(nsr)
		r₂₁_XY[i] = coef.r₂₁(nsr)
		t₂₁_XY[i] = coef.t₂₁(nsr)
	end
	return PropagationCoefficientScalar{T,Array{Complex{T},2}}(r₁₂_XY, t₁₂_XY, r₂₁_XY, t₂₁_XY, coef.λ, coef.n₁, coef.ref₁, coef.n₂, coef.ref₂)
end

function coefficient_itpform(coef::PropagationCoefficientScalar{T}, x::AbstractVector) where T
	r₁₂ = extrapolation(coef.r₁₂, x)
	t₁₂ = extrapolation(coef.t₁₂, x)
	r₂₁ = extrapolation(coef.r₂₁, x)
	t₂₁ = extrapolation(coef.t₂₁, x)
	PropagationCoefficientScalar{T}(r₁₂, t₁₂, r₂₁, t₂₁, coef.λ, coef.n₁, coef.ref₁, coef.n₂, coef.ref₂)
end

function coefficient_itpform(coef::PropagationCoefficientVectorial{T}, x::AbstractVector) where T
	rpp₁₂ = extrapolation(coef.rpp₁₂, x)
	tpp₁₂ = extrapolation(coef.tpp₁₂, x)
	rpp₂₁ = extrapolation(coef.rpp₂₁, x)
	tpp₂₁ = extrapolation(coef.tpp₂₁, x)
	rss₁₂ = extrapolation(coef.rss₁₂, x)
	tss₁₂ = extrapolation(coef.tss₁₂, x)
	rss₂₁ = extrapolation(coef.rss₂₁, x)
	tss₂₁ = extrapolation(coef.tss₂₁, x)
	return PropagationCoefficientVectorial{T}(rpp₁₂, rss₁₂, tpp₁₂, tss₁₂, rpp₂₁, rss₂₁, tpp₂₁, tss₂₁, coef.λ, coef.n₁, coef.ref₁, coef.n₂, coef.ref₂)
end

@inline function intensity(e_SXY::AbstractArray{Complex{T}}) where {T<:Number}
	val = zero(T)
	@inbounds @simd for i in eachindex(e_SXY)
		val += abs2(e_SXY[i])
	end
	return val
end

function checkmergeble(coefs::AbstractArray{<:AbstractPropagationCoefficient})
	@inbounds for i in 1:length(coefs)-1
		(coefs[i].λ - coefs[i+1].λ < @tol) || return false;
		checkorientation(coefs[i].ref₁, coefs[i+1].ref₁) || return false
		(coefs[i+1].ref₁.z - coefs[i].ref₂.z > -@tol) || return false
		checkinline(coefs[i].ref₁, coefs[i+1].ref₁) || return false
		(abs(coefs[i+1].n₁ - coefs[i].n₂) < @tol) || return false
	end
	@inbounds for i in eachindex(coefs)
		checkinline(coefs[i].ref₁, coefs[i].ref₂) || return false
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
	@inbounds for i in eachindex(coefs)
		checkinline(coefs[i].ref₁, coefs[i].ref₂) || return false
	end
	return true
end

function r₁₂(nsr::Real, coefs::AbstractVector{PropagationCoefficientScalar{T,X}}) where {T,X}
 	k = 2π / coefs[1].λ;

	ri = coefs[end].r₁₂(nsr)
	xref = coefs[end].ref₁.x;
	yref = coefs[end].ref₁.y;
	zref = coefs[end].ref₁.z;
	for i in length(coefs)-1:-1:1
		d = √((xref - coefs[i].ref₂.x)^2 + (yref - coefs[i].ref₂.y)^2 + (zref - coefs[i].ref₂.z)^2)
		tmpphase = exp(im * k * d * √(coefs[i].n₂^2 - nsr^2));
		ri = coefs[i].r₁₂(nsr) + coefs[i].t₁₂(nsr) * coefs[i].t₂₁(nsr) * ri * tmpphase^2 / (1 - coefs[i].r₂₁(nsr) * ri * tmpphase^2)
		xref = coefs[i].ref₁.x;
		yref = coefs[i].ref₁.y;
		zref = coefs[i].ref₁.z;
	end
	return ri
end

function t₁₂(nsr::Real, coefs::AbstractVector{PropagationCoefficientScalar{T,X}}) where {T,X}
 	k = 2π / coefs[1].λ;
	ri = coefs[end].r₁₂(nsr)
	ti = coefs[end].t₁₂(nsr)
	xref = coefs[end].ref₁.x;
	yref = coefs[end].ref₁.y;
	zref = coefs[end].ref₁.z;
	for i in length(coefs)-1:-1:1
		d = √((xref - coefs[i].ref₂.x)^2 + (yref - coefs[i].ref₂.y)^2 + (zref - coefs[i].ref₂.z)^2)
		tmpphase = exp(im * k * d * √(coefs[i].n₂^2 - nsr^2));
		ti = coefs[i].t₁₂(nsr) * ti * tmpphase / (1 - coefs[i].r₂₁.(nsr) * ri * tmpphase^2)
		ri = coefs[i].r₁₂(nsr) + coefs[i].t₁₂(nsr) * coefs[i].t₂₁(nsr) * ri * tmpphase^2 / (1 - coefs[i].r₂₁(nsr) * ri * tmpphase^2)
		xref = coefs[i].ref₁.x;
		yref = coefs[i].ref₁.y;
		zref = coefs[i].ref₁.z;
	end
	return ti
end

function t₂₁(nsr::Real, coefs::AbstractVector{PropagationCoefficientScalar{T,X}}) where {T,X}
	ri = coefs[1].r₂₁(nsr)
	ti = coefs[1].t₂₁(nsr)
 	k = 2π / coefs[1].λ;
	xref = coefs[1].ref₂.x;
	yref = coefs[1].ref₂.y;
	zref = coefs[1].ref₂.z;
	@inbounds for i in 2:1:length(coefs)
		d = √((xref - coefs[i].ref₁.x)^2 + (yref - coefs[i].ref₁.y)^2 + (zref - coefs[i].ref₁.z)^2)
		tmpphase = exp(im * k * d * √(coefs[i].n₁^2 - nsr^2));

		ti = coefs[i].t₂₁(nsr) * ti * tmpphase / (1 - coefs[i].r₁₂(nsr) * ri * tmpphase^2)
		ri = coefs[i].r₂₁(nsr) + coefs[i].t₂₁(nsr) * coefs[i].t₁₂(nsr) * ri * tmpphase^2 / (1 - coefs[i].r₁₂(nsr) * ri * tmpphase^2)
		xref = coefs[i].ref₂.x;
		yref = coefs[i].ref₂.y;
		zref = coefs[i].ref₂.z;
	end
	return ti
end

function r₂₁(nsr::Real, coefs::AbstractVector{PropagationCoefficientScalar{T,X}}) where {T,X}
	ri = coefs[1].r₂₁(nsr)
 	k = 2π / coefs[1].λ;
	xref = coefs[1].ref₂.x;
	yref = coefs[1].ref₂.y;
	zref = coefs[1].ref₂.z;
	@inbounds for i in 2:1:length(coefs)
		d = √((xref - coefs[i].ref₁.x)^2 + (yref - coefs[i].ref₁.y)^2 + (zref - coefs[i].ref₁.z)^2)
		tmpphase = exp(im * k * d * √(coefs[i].n₁^2 - nsr^2));

		ri = coefs[i].r₂₁(nsr) + coefs[i].t₂₁(nsr) * coefs[i].t₁₂(nsr) * ri * tmpphase^2 / (1 - coefs[i].r₁₂(nsr) * ri * tmpphase^2)
		xref = coefs[i].ref₂.x;
		yref = coefs[i].ref₂.y;
		zref = coefs[i].ref₂.z;
	end
	return ri
end

function mergeorientated_propagationcoefficient(coefs::AbstractVector{<:PropagationCoefficientScalar{T}}) where {T<:Real}
	for i in 1:length(coefs)-1
		(abs(coefs[i+1].n₁ - coefs[i].n₂) < @tol) || error("There is a discontinuty on the refractive index.")
	end
 	checkmergeble(coefs) || error("The structures are not parallel to each other. Use lightinteraction_recursive instead.")
	r12(nsr) = r₁₂(nsr, coefs)
	t12(nsr) = t₁₂(nsr, coefs)
	r21(nsr) = r₂₁(nsr, coefs)
	t21(nsr) = t₂₁(nsr, coefs)
	return PropagationCoefficientScalar{T}(r12, t12, r21, t21, coefs[1].λ, coefs[1].n₁, coefs[1].ref₁, coefs[end].n₂, coefs[end].ref₂)
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

function lightinteraction_recursive!(ebackward_SCD::AbstractArray{<:T, 3}, eforward_SCD::AbstractArray{<:T, 3}, coefs::AbstractVector{<:AbstractPropagationCoefficient}, nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number,2}, ei_SXY::AbstractArray{<:Number, 3}, λ::Real, dir::Integer, nsxOut::AbstractVector{<:Real}, nsyOut::AbstractVector{<:Real}, thresold::Real, sizeM::Integer) where T<:Number

	numberInterfaces = length(coefs);
	(numberInterfaces < 3) && (sizeM = 5)
	(sizeS, sizeX, sizeY) = size(ei_SXY);
	sizeC, sizeD = length(nsxOut), length(nsyOut);
	(sizeX, sizeY) == size(nsx_XY) == size(nsy_XY) || error("Wrong sizes")
	sizeS == size(ebackward_SCD, 1) == size(eforward_SCD, 1) || error("wrong sizes")
	sizeC == size(ebackward_SCD, 2) == size(eforward_SCD, 2) || error("wrong sizes")
	sizeD == size(ebackward_SCD, 3) == size(eforward_SCD, 3) || error("wrong sizes")

	ESave_SXYM = Array{T, 4}(undef, sizeS, sizeX, sizeY, sizeM);
	nsxSave_XYM = Array{T, 3}(undef, sizeX, sizeY, sizeM);
	nsySave_XYM = Array{T, 3}(undef, sizeX, sizeY, sizeM);
	nszSave_XYM = Array{T, 3}(undef, sizeX, sizeY, sizeM);
	mlsSave_M = Vector{Int64}(undef, sizeM);

	EOutMLSForward_SXY = Array{T, 3}(undef, sizeS, sizeX, sizeY);
	nsxOutMLSForward_XY = Array{T, 2}(undef, sizeX, sizeY);
	nsyOutMLSForward_XY = Array{T, 2}(undef, sizeX, sizeY);
	nszOutMLSForward_XY = Array{T, 2}(undef, sizeX, sizeY);
	EtoFileForward_SXY = Array{T, 3}(undef, sizeS, sizeX, sizeY);
	nsxtoFileForward_XY = Array{T, 2}(undef, sizeX, sizeY);
	nsytoFileForward_XY = Array{T, 2}(undef, sizeX, sizeY);
	nsztoFileForward_XY = Array{T, 2}(undef, sizeX, sizeY);

	EOutMLSBackward_SXY = Array{T, 3}(undef, sizeS, sizeX, sizeY);
	nsxOutMLSBackward_XY = Array{T, 2}(undef, sizeX, sizeY);
	nsyOutMLSBackward_XY = Array{T, 2}(undef, sizeX, sizeY);
	nszOutMLSBackward_XY = Array{T, 2}(undef, sizeX, sizeY);

	EtoFileBackward_SXY = Array{T, 3}(undef, sizeS, sizeX, sizeY);
	nsxtoFileBackward_XY = Array{T, 2}(undef, sizeX, sizeY);
	nsytoFileBackward_XY = Array{T, 2}(undef, sizeX, sizeY);
	nsztoFileBackward_XY = Array{T, 2}(undef, sizeX, sizeY);

	EfromFile_SXY = Array{T, 3}(undef, sizeS, sizeX, sizeY);
	nsxfromFile_XY = Array{T, 2}(undef, sizeX, sizeY);
	nsyfromFile_XY = Array{T, 2}(undef, sizeX, sizeY);
	nszfromFile_XY = Array{T, 2}(undef, sizeX, sizeY);
	tmpnₛ_3 = MVector{3,T}(undef)
	tmpnₚ_3 = MVector{3,T}(undef)

	EsaveGrided_SCD = Array{T, 3}(undef, sizeS, sizeC, sizeD);
	ebackward_SCD .= zero(T);
	eforward_SCD .= zero(T);
	boolOutForward = false; boolOutBackward = false;

	iSaveMem = 1; iReadMem = 0;
	iSaveFile = 0; iReadFile = 0;
	isForward = dir > 0;

	ESave_SXYM[:,:,:,1] .= ei_SXY
	mlsSave_M[1] = isForward ? 1 : (numberInterfaces + 1)
	incidentbeamintensity = 0.
	@inbounds @simd for i in eachindex(ei_SXY)
		incidentbeamintensity += abs2(ei_SXY[i]);
	end
	thresold *= thresold

	n = dir > 0 ? coefs[1].n₁ : coefs[end].n₂;
	@inbounds @simd for iXY in eachindex(nsx_XY)
			nsxSave_XYM[iXY] = nsx_XY[iXY]
			nsySave_XYM[iXY] = nsy_XY[iXY]
			nszSave_XYM[iXY] = dir * √(n^2 - nsx_XY[iXY]^2 - nsy_XY[iXY]^2);
	end

	filename = tempdir() * "/JolabTMP" * randstring(48) * "_";

	i = 1
	booltoFileForward = false; booltoFileBackward = false;
	@inbounds while (iReadMem < iSaveMem) || (iReadFile < iSaveFile)
		iWrapSaveMem = (iSaveMem % sizeM) + 1; iWrapReadMem = (iReadMem % sizeM) + 1;
		if iReadMem < iSaveMem
			mlsᵢ = mlsSave_M[iWrapReadMem];
			Eᵢ_SXY = @view ESave_SXYM[:,:,:,iWrapReadMem];
			nsxᵢ_XY = @view nsxSave_XYM[:,:,iWrapReadMem];
			nsyᵢ_XY = @view nsySave_XYM[:,:,iWrapReadMem];
			nszᵢ_XY = @view nszSave_XYM[:,:,iWrapReadMem];
			loadFromMemory = true
			iReadMem += 1;
		else
			mlsᵢ = readdata(filename * string(iReadFile), EfromFile_SXY, nsxfromFile_XY, nsyfromFile_XY, nszfromFile_XY)
			Eᵢ_SXY = EfromFile_SXY;
			nsxᵢ_XY = nsxfromFile_XY;
			nsyᵢ_XY = nsyfromFile_XY;
			nszᵢ_XY = nszfromFile_XY;
			loadFromMemory = false
			iReadFile += 1;
		end

		roundtripintensity = 0.;
		@simd for i in eachindex(Eᵢ_SXY)
			roundtripintensity += abs2(Eᵢ_SXY[i]);
		end
		roundtripintensity < incidentbeamintensity * thresold && continue
		roundtripintensity > incidentbeamintensity && (println("Error on the simulation. Power increased on recursive algorithm"); break)

		saveTableMemFull = (iWrapSaveMem == iWrapReadMem - 2 || iWrapSaveMem == iWrapReadMem - 1);
		isForward = real(nszᵢ_XY[Integer(round(sizeX * sizeY / 2))]) > 0;
		if (mlsᵢ == 1 && isForward) || (mlsᵢ == 2 && !isForward)
			ESaveBackward_SXY = EOutMLSBackward_SXY;
			nsxSaveBackward_XY = nsxOutMLSBackward_XY;
			nsySaveBackward_XY = nsyOutMLSBackward_XY;
			nszSaveBackward_XY = nszOutMLSBackward_XY;
			boolOutBackward = true;
		else
			if saveTableMemFull
				ESaveBackward_SXY = EtoFileBackward_SXY;
				nsxSaveBackward_XY = nsxtoFileBackward_XY;
				nsySaveBackward_XY = nsytoFileBackward_XY;
				nszSaveBackward_XY = nsztoFileBackward_XY;
				mlstoFileBackward = isForward ? mlsᵢ : mlsᵢ - 1
				booltoFileBackward = true;
			else
				ESaveBackward_SXY = @view ESave_SXYM[:,:,:,iWrapSaveMem];
				nsxSaveBackward_XY = @view nsxSave_XYM[:,:,iWrapSaveMem];
				nsySaveBackward_XY = @view nsySave_XYM[:,:,iWrapSaveMem];
				nszSaveBackward_XY = @view nszSave_XYM[:,:,iWrapSaveMem];
				mlsSave_M[iWrapSaveMem] = isForward ? mlsᵢ : mlsᵢ - 1
				iSaveMem += 1
			end
		end
		if (mlsᵢ == numberInterfaces && isForward) || (mlsᵢ == numberInterfaces + 1 && !isForward)
			ESaveForward_SXY = EOutMLSForward_SXY;
			nsxSaveForward_XY = nsxOutMLSForward_XY;
			nsySaveForward_XY = nsyOutMLSForward_XY;
			nszSaveForward_XY = nszOutMLSForward_XY;
			boolOutForward = true;
		else
			if saveTableMemFull
				ESaveForward_SXY = EtoFileForward_SXY;
				nsxSaveForward_XY = nsxtoFileForward_XY;
				nsySaveForward_XY = nsytoFileForward_XY;
				nszSaveForward_XY = nsztoFileForward_XY;
				mlstoFileForward = isForward ? mlsᵢ + 1 : mlsᵢ
				booltoFileForward = true;
			else
				iWrapSaveMem = (iSaveMem % sizeM) + 1;
				ESaveForward_SXY = @view ESave_SXYM[:,:,:,iWrapSaveMem];
				nsxSaveForward_XY = @view nsxSave_XYM[:,:,iWrapSaveMem];
				nsySaveForward_XY = @view nsySave_XYM[:,:,iWrapSaveMem];
				nszSaveForward_XY = @view nszSave_XYM[:,:,iWrapSaveMem];
				mlsSave_M[iWrapSaveMem] = isForward ? mlsᵢ + 1 : mlsᵢ
				iSaveMem += 1
			end
		end
		if isForward
			mlsᵢ > 1 && changereferential!(nsxᵢ_XY, nsyᵢ_XY, nszᵢ_XY, Eᵢ_SXY, λ, coefs[mlsᵢ-1].ref₂, coefs[mlsᵢ].ref₁)
			lightinteraction!(ESaveBackward_SXY, ESaveForward_SXY, coefs[mlsᵢ], nsxᵢ_XY, nsyᵢ_XY, Eᵢ_SXY, 1, tmpnₚ_3, tmpnₛ_3)
			nsxSaveBackward_XY .= nsxᵢ_XY;
			nsySaveBackward_XY .= nsyᵢ_XY;
			nszSaveBackward_XY .= .-nszᵢ_XY;

			nsxSaveForward_XY .= nsxᵢ_XY;
			nsySaveForward_XY .= nsyᵢ_XY;
			nszSaveForward_XY .= .√(coefs[mlsᵢ].n₂.^2 .- nsxᵢ_XY.^2 .- nsyᵢ_XY.^2);
		else
			mlsᵢ< (numberInterfaces + 1) && changereferential!(nsxᵢ_XY, nsyᵢ_XY, nszᵢ_XY, Eᵢ_SXY, λ, coefs[mlsᵢ].ref₁, coefs[mlsᵢ-1].ref₂);
			lightinteraction!(ESaveBackward_SXY, ESaveForward_SXY, coefs[mlsᵢ-1], nsxᵢ_XY, nsyᵢ_XY, Eᵢ_SXY, -1, tmpnₚ_3, tmpnₛ_3)
			nsxSaveForward_XY .= nsxᵢ_XY;
			nsySaveForward_XY .= nsyᵢ_XY;
			nszSaveForward_XY .= .-nszᵢ_XY;

			nsxSaveBackward_XY .= nsxᵢ_XY;
			nsySaveBackward_XY .= nsyᵢ_XY;
			nszSaveBackward_XY .= -1 .* .√(coefs[mlsᵢ].n₂.^2 .- nsxᵢ_XY.^2 .- nsyᵢ_XY.^2);
		end

		if boolOutForward
			interpolatetogrid!(EsaveGrided_SCD, real.(nsxSaveForward_XY), real.(nsySaveForward_XY), ESaveForward_SXY, nsxOut, nsyOut);
			eforward_SCD .+= EsaveGrided_SCD;
			boolOutForward = false
		end
		if boolOutBackward
			interpolatetogrid!(EsaveGrided_SCD, real.(nsxSaveBackward_XY), real.(nsySaveBackward_XY), ESaveBackward_SXY, nsxOut, nsyOut);
			ebackward_SCD .+= EsaveGrided_SCD;
			boolOutBackward = false
		end
		if booltoFileForward
			writedata(filename * string(iSaveFile), EtoFileForward_SXY, nsxtoFileForward_XY, nsytoFileForward_XY, nsztoFileForward_XY, mlstoFileForward)
			booltoFileForward = false
			iSaveFile += 1;
		end
		if booltoFileBackward
			writedata(filename * string(iSaveFile), EtoFileBackward_SXY, nsxtoFileBackward_XY, nsytoFileBackward_XY, nsztoFileBackward_XY, mlstoFileBackward)
			booltoFileBackward = false
			iSaveFile += 1;
		end
		i += 1
	end
end

function lightinteraction_recursive(coefs::AbstractVector{<:AbstractPropagationCoefficient}, angspe::AbstractFieldAngularSpectrum{T}, nsxOut::AbstractArray{<:Number}, nsyOut::AbstractArray{<:Number}; thresold = 1E-4::Real, sizeM = 100::Integer) where T
	checkapplicability(coefs, angspe);
	(sizeS, sizeX, sizeY) = size(angspe.e_SXY);
	sizeC, sizeD  = length(nsxOut), length(nsyOut)

	ebackward_SCD = Array{Complex{T},3}(undef, sizeS, sizeC, sizeD);
	eforward_SCD = Array{Complex{T},3}(undef, sizeS, sizeC, sizeD);
	ei_SXY = copy(angspe.e_SXY);

	nsx_XY = repeat(angspe.nsx_X, 1, sizeY);
	nsy_XY = repeat(angspe.nsy_Y', sizeX, 1);
	nsz_XY = angspe.dir * .√(angspe.n^2 .- nsx_XY.^2 .- nsy_XY.^2)

	refcoef = angspe.dir > 0 ? coefs[1].ref₁ : coefs[end].ref₂
	changereferential!(nsx_XY, nsy_XY, nsz_XY, ei_SXY, angspe.λ, angspe.ref, refcoef)

	(nsrmin, nsrmax) = rextrema(angspe.nsx_X, angspe.nsy_Y);
	(nsrminC, nsrmaxC) = rextrema(nsxOut, nsyOut);
	nsr = range(min(nsrmin, nsrminC), stop = max(nsrmax, nsrmaxC), length = round(Int, max(sizeX, sizeY, sizeC, sizeD) / √2))

	coefs_itp = [coefficient_itpform(coefi, nsr) for coefi in coefs]

	lightinteraction_recursive!(ebackward_SCD, eforward_SCD, coefs_itp, nsx_XY, nsy_XY, ei_SXY, angspe.λ, angspe.dir, nsxOut, nsyOut, thresold, sizeM);

	angspeforward = FieldAngularSpectrum{T}(nsxOut, nsyOut, eforward_SCD, angspe.λ, coefs[end].n₂, 1, coefs[end].ref₂)
	angspebackward = FieldAngularSpectrum{T}(nsxOut, nsyOut, ebackward_SCD, angspe.λ, coefs[1].n₁, -1, coefs[1].ref₁)
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

@inline function checkapplicability(coefs::AbstractVector{<:AbstractCoefficient})
	@inbounds @simd for i in 2:length(coefs)
		(abs(coefs[i-1].n₂ - coefs[i].n₁) < @tol) || error("Refractive index are different")
		coefs[i-1].ref₂.z < coefs[i].ref₁.z || error("The order of the referential along the optical axis is incorrect")
	end
	return true;
end

@inline function checkapplicability(coefs::AbstractVector{<:AbstractCoefficient}, angspe::AbstractFieldAngularSpectrum)
	angspe.dir > 0 ? ((abs(angspe.n - coefs[1].n₁) < @tol) || error("refractive index of coefs and field are differente")) : ((abs(angspe.n - coefs[end].n₂) < @tol) || error("Refractive index of coefs and field are different"))
	@inbounds @simd for i in eachindex(coefs)
		(abs(angspe.λ - coefs[i].λ) < @tol) || error("λ of coefs and field are different")
	end
	return checkapplicability(coefs);
end

@inline function writedata(fname::String, E::AbstractArray, nsx::AbstractArray, nsy::AbstractArray, nsz::AbstractArray, mls::Integer)
	f = open(fname * ".bin", "w");
	for i in E
		write(f, i)
	end
	for i in nsx
		write(f, i)
	end
	for i in nsy
		write(f, i)
	end
	for i in nsz
		write(f, i)
	end
	write(f, mls)
	close(f);
end

@inline function readdata(fname::String, E::AbstractArray, nsx::AbstractArray, nsy::AbstractArray, nsz::AbstractArray)
	f = open(fname * ".bin", "r")
	read!(f, E)
	read!(f, nsx)
	read!(f, nsy)
	read!(f, nsz)
	mls = reinterpret(Int64, read(f))[1]
	close(f)
	rm(fname * ".bin")
	return mls
end

function rextrema(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
	min = x[1]^2 + y[1]^2
	max = x[1]^2 + y[1]^2
	@inbounds @simd for xi in x
		for yi in y
			r = xi^2 + yi^2
			(min > r) && (min = r)
			(max < r) && (max = r)
		end
	end
	return (√min, √max)
end

function mergeorientated_propagationcoefficient(coefs::AbstractVector{<:PropagationCoefficientScalar{T, A}}, nsx_XY::AbstractArray{<:Real}, nsy_XY::AbstractArray{<:Real}) where {T<:Real, A <: AbstractArray{T}}
	checkmergeble(coefs) || error("Coefficients cannot be merged")
	sizeA = length(coefs)
	r12 = coefs[sizeA].r₁₂
	t12 = coefs[sizeA].t₁₂
	nsz_XY = Array{Complex{T}, 1}(undef, length(nsx_XY))
	propMatrix = Diagonal(Vector{Complex{T}}(undef, length(nsx_XY)))

	for i in length(coefs)-1:-1:1
		nsz_XY .= .√(coefs[i].n₂^2 .- nsx_XY.^2 .- nsy_XY.^2)
		propagationmatrix!(propMatrix, nsx_XY, nsy_XY, nsz_XY, coefs[1].λ, coefs[i].ref₂, coefs[i+1].ref₁)
		lmul!(r12, propMatrix)
		rmul!(propMatrix, r12)
		# r12 is r23 inside the equation
		r12 = coefs[i].r₁₂ + coefs[i].t₂₁ * inv(I - r12 * coefs[i].r₂₁) * r12 * coefs[i].t₁₂
		@show typeof(r12)
	end
	return r12
end

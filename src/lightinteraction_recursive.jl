function lightinteraction_recursive!(ebackward_SCD::AbstractArray{<:T, 3}, eforward_SCD::AbstractArray{<:T, 3}, coefs::AbstractVector{<:AbstractCoefficient}, nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number,2}, ei_SXY::AbstractArray{<:Number, 3}, λ::Real, dir::Integer, nsxOut::AbstractVector{<:Real}, nsyOut::AbstractVector{<:Real}, thresold::Real, sizeM::Integer) where T<:Number

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
			mlsᵢ > 1 && changereferenceframe!(nsxᵢ_XY, nsyᵢ_XY, nszᵢ_XY, Eᵢ_SXY, λ, coefs[mlsᵢ-1].ref₂, coefs[mlsᵢ].ref₁)
			lightinteraction!(ESaveBackward_SXY, ESaveForward_SXY, coefs[mlsᵢ], nsxᵢ_XY, nsyᵢ_XY, Eᵢ_SXY, 1, tmpnₚ_3, tmpnₛ_3)
			nsxSaveBackward_XY .= nsxᵢ_XY;
			nsySaveBackward_XY .= nsyᵢ_XY;
			nszSaveBackward_XY .= .-nszᵢ_XY;

			nsxSaveForward_XY .= nsxᵢ_XY;
			nsySaveForward_XY .= nsyᵢ_XY;
			nszSaveForward_XY .= .√(coefs[mlsᵢ].n₂.^2 .- nsxᵢ_XY.^2 .- nsyᵢ_XY.^2);
		else
			mlsᵢ< (numberInterfaces + 1) && changereferenceframe!(nsxᵢ_XY, nsyᵢ_XY, nszᵢ_XY, Eᵢ_SXY, λ, coefs[mlsᵢ].ref₁, coefs[mlsᵢ-1].ref₂);
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

function lightinteraction_recursive(coefs::AbstractVector{<:AbstractCoefficient}, angspe::AbstractFieldAngularSpectrum{T}, nsxOut::AbstractArray{<:Number}, nsyOut::AbstractArray{<:Number}; thresold = 1E-4::Real, sizeM = 100::Integer) where T
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
	changereferenceframe!(nsx_XY, nsy_XY, nsz_XY, ei_SXY, angspe.λ, angspe.ref, refcoef)

	(nsrmin, nsrmax) = rextrema(angspe.nsx_X, angspe.nsy_Y);
	(nsrminC, nsrmaxC) = rextrema(nsxOut, nsyOut);
	nsr = range(min(nsrmin, nsrminC), stop = max(nsrmax, nsrmaxC), length = round(Int, max(sizeX, sizeY, sizeC, sizeD) / √2))

	coefs_itp = [coefficient_itpform(coefi, nsr) for coefi in coefs]

	lightinteraction_recursive!(ebackward_SCD, eforward_SCD, coefs_itp, nsx_XY, nsy_XY, ei_SXY, angspe.λ, angspe.dir, nsxOut, nsyOut, thresold, sizeM);

	angspeforward = FieldAngularSpectrum{T}(nsxOut, nsyOut, eforward_SCD, angspe.λ, coefs[end].n₂, 1, coefs[end].ref₂)
	angspebackward = FieldAngularSpectrum{T}(nsxOut, nsyOut, ebackward_SCD, angspe.λ, coefs[1].n₁, -1, coefs[1].ref₁)
	return (angspebackward, angspeforward)
end

@inline function checkapplicability(coefs::AbstractVector{<:AbstractCoefficient})
	@inbounds @simd for i in 2:length(coefs)
		(abs(coefs[i-1].n₂ - coefs[i].n₁) < @tol) || error("Refractive index are different")
		(coefs[i-1].ref₂.z - coefs[i].ref₁.z < @tol) || error("The order of the referenceframe along the optical axis is incorrect")
		(abs(coefs[i-1].λ - coefs[i].λ) < @tol) || error("The wavelength of the coefficients are different. Cannot be merged")
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

function lightinteraction_recursive(components::AbstractVector{T}, angspe::AbstractFieldAngularSpectrum; nsxOut = angspe.nsx_X::AbstractVector{<:Real}, nsyOut = angspe.nsy_Y::AbstractVector{<:Real}, thresold = 1E-4, sizeM=100::Integer) where T
	vectorial = (size(angspe.e_SXY, 1) == 3)
	vectorial && tobedone()

	coefs = [coefficientscallar(compi, angspe.λ) for compi in components]
	return lightinteraction_recursive(coefs, angspe, nsxOut, nsyOut, thresold = thresold, sizeM = sizeM);
end

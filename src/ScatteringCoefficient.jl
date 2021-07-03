struct RoughInterfaceConvolutionCoefficient{T,L,R, Y<:AbstractVector{Complex{T}}, X<:AbstractArray{Complex{T},2},F <: AbstractFFTs.Plan}  <: AbstractCoefficient{T,L,R}
	r₁₂::Y
	t₁₂::Y
	r₂₁::Y
	t₂₁::Y
	ir₁₂::X
	sr₁₂::X
	it₁₂::X
	st₁₂::X
	ir₂₁::X
	sr₂₁::X
	it₂₁::X
	st₂₁::X
	Δ::X
	fieldl::L
	fieldr::R
	planfft::F
	tmp1::X
	tmp2::X
end

@inbounds function lightinteraction!(fieldl::L, fieldr::R, coef::RoughInterfaceConvolutionCoefficient{T,L,R}, fieldi::Union{L,R}) where {T,L,R}
	samedefinitions(fieldi, dir(fieldi) > 0 ? fieldl : fieldr) || tobedone()
	sizeX, sizeY = length(fieldi.nsx_X), length(fieldi.nsy_Y)
	sizeM = size(coef.ir₁₂, 2)
	if dir(fieldi) > 0
		coef.tmp1 .= fieldi.e_SXY .* coef.ir₁₂
		Threads.@threads for iM in 1:sizeM
			mul!(reshape(view(coef.tmp2,:,iM), sizeX, sizeY), coef.planfft, reshape(view(coef.tmp1,:,iM), sizeX, sizeY))
		end
		coef.tmp2 .*= coef.Δ
		Threads.@threads for iM in 1:sizeM
			ldiv!(reshape(view(coef.tmp1,:,iM), sizeX, sizeY), coef.planfft, reshape(view(coef.tmp2,:,iM), sizeX, sizeY))
		end
		coef.tmp1 .*= coef.sr₁₂
		@simd for i in iterator_index(fieldl)
			fieldl.e_SXY[i] = coef.tmp1[i,1]
		end
		for iM in 2:sizeM
			for i in iterator_index(fieldl)
				fieldl.e_SXY[i] += coef.tmp1[i,iM]
			end
		end
		coef.tmp1 .= fieldi.e_SXY .* coef.it₁₂;
		Threads.@threads for iM in 1:sizeM
			mul!(reshape(view(coef.tmp2,:,iM),sizeX, sizeY), coef.planfft, reshape(view(coef.tmp1,:,iM), sizeX, sizeY))
		end
		coef.tmp2 .*= coef.Δ
		Threads.@threads for iM in 1:sizeM
			ldiv!(reshape(view(coef.tmp1,:,iM), sizeX, sizeY), coef.planfft, reshape(view(coef.tmp2,:,iM), sizeX, sizeY))
		end
		coef.tmp1 .*= coef.st₁₂
		@simd for i in iterator_index(fieldr)
			fieldr.e_SXY[i] = coef.tmp1[i,1]
		end
		@simd for iM in 2:sizeM
			for i in iterator_index(fieldr)
				fieldr.e_SXY[i] += coef.tmp1[i,iM]
			end
		end
	else
		coef.tmp1 .= fieldi.e_SXY .* coef.ir₂₁
		Threads.@threads for iM in 1:sizeM
			mul!(reshape(view(coef.tmp2,:,iM), sizeX, sizeY), coef.planfft, reshape(view(coef.tmp1,:,iM), sizeX, sizeY))
		end
		@simd for iM in 1:sizeM
			for i in 1:sizeX * sizeY
				coef.tmp2[i,iM] *= coef.Δ[i,sizeM-iM+1]
			end
		end
		Threads.@threads for iM in 1:sizeM
			ldiv!(reshape(view(coef.tmp1,:,iM), sizeX, sizeY), coef.planfft, reshape(view(coef.tmp2,:,iM), sizeX, sizeY))
		end
		coef.tmp1 .*= coef.sr₂₁
		@simd for i in iterator_index(fieldr)
			fieldr.e_SXY[i] = coef.tmp1[i,1]
		end
		@simd for iM in 2:sizeM
			for i in iterator_index(fieldr)
				fieldr.e_SXY[i] += coef.tmp1[i,iM]
			end
		end
		coef.tmp1 .= fieldi.e_SXY .* coef.it₂₁

		Threads.@threads for iM in 1:sizeM
			mul!(reshape(view(coef.tmp2,:,iM), sizeX, sizeY), coef.planfft, reshape(view(coef.tmp1,:,iM), sizeX, sizeY))
		end
		@simd for iM in 1:sizeM
			for i in 1:sizeX * sizeY
				coef.tmp2[i,iM] *= coef.Δ[i,sizeM-iM+1]
			end
		end
		Threads.@threads for iM in 1:sizeM
			ldiv!(reshape(view(coef.tmp1,:,iM), sizeX, sizeY), coef.planfft, reshape(view(coef.tmp2,:,iM), sizeX, sizeY))
		end
		coef.tmp1 .*= coef.st₂₁
		@simd for i in iterator_index(fieldl)
			fieldl.e_SXY[i] = coef.tmp1[i,1]
		end
		@simd for iM in 2:sizeM
			for i in iterator_index(fieldl)
				fieldl.e_SXY[i] += coef.tmp1[i,iM]
			end
		end
	end
	if dir(fieldi) > 0
		@simd for i in iterator_index(fieldl)
			fieldl.e_SXY[i] += fieldi.e_SXY[i] * coef.r₁₂[i]
		end
		@simd for i in iterator_index(fieldr)
			fieldr.e_SXY[i] += fieldi.e_SXY[i] * coef.t₁₂[i]
		end
	else
		@simd for i in iterator_index(fieldr)
			fieldr.e_SXY[i] += fieldi.e_SXY[i] * coef.r₂₁[i]
		end
		@simd for i in iterator_index(fieldl)
			fieldl.e_SXY[i] += fieldi.e_SXY[i] * coef.t₂₁[i]
		end
	end
end

function correctscatteringmatrix_referenceframes!(scat::RoughInterfaceConvolutionCoefficient, comp::AbstractOpticalComponent, fieldi::AbstractFieldMonochromatic)
	ref = dir(fieldi) > 0 ? ref1(comp) : ref2(comp)
	if !checkposition(fieldi.ref, ref)
		propM = propagationmatrix(fieldi, ref)
		if dir(fieldi) > 0
			scat.r₁₂ .*= propM.diag .* propM.diag
			scat.t₁₂ .*= propM.diag
			scat.t₂₁ .*= propM.diag
			scat.ir₁₂ .*= propM.diag
			scat.sr₁₂ .*= propM.diag
			scat.it₁₂ .*= propM.diag
			scat.st₂₁ .*= propM.diag
			scat.fieldl.ref = copy(fieldi.ref)
		else
			# conj!(propM.diag)
			scat.r₂₁ .*= propM.diag .* propM.diag
			scat.t₁₂ .*= propM.diag
			scat.t₂₁ .*= propM.diag
			scat.ir₂₁ .*= propM.diag
			scat.sr₂₁ .*= propM.diag
			scat.st₁₂ .*= propM.diag
			scat.it₂₁ .*= propM.diag
			scat.fieldr.ref = copy(fieldi.ref)
		end
	end
end

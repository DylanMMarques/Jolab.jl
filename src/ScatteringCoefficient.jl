struct RoughInterfaceConvolutionCoefficient{T,L,R, Y<:AbstractArray{Complex{T}}, X<:AbstractArray{Complex{T}},F <: AbstractFFTs.Plan}  <: AbstractCoefficient{T,L,R}
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
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	sizeM = size(coef.ir₁₂, 3)
	if dir(fieldi) > 0
		Threads.@threads for iM in 1:sizeM
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp1[iX,iY,iM] = fieldi.e_SXY[1,iX,iY] * coef.ir₁₂[iX,iY,iM]
				end
			end
			mul!(view(coef.tmp2,:,:,iM), coef.planfft, view(coef.tmp1,:,:,iM))
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp2[iX,iY,iM] *= coef.Δ[iX,iY,iM]
				end
			end
			ldiv!(view(coef.tmp1,:,:,iM), coef.planfft, view(coef.tmp2,:,:,iM))
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					if iM == 1
						fieldl.e_SXY[1,iX,iY] = coef.tmp1[iX,iY,1] * coef.sr₁₂[iX,iY,1]
					else
						coef.tmp1[iX,iY,iM] *= coef.sr₁₂[iX,iY,iM]
					end
				end
			end
		end
		@simd for iM in 2:sizeM
			for iY in 1:sizeY
				for iX in 1:sizeX
					fieldl.e_SXY[1,iX,iY] += coef.tmp1[iX,iY,iM]
				end
			end
		end
		Threads.@threads for iM in 1:sizeM
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp1[iX,iY,iM] = fieldi.e_SXY[1,iX,iY] * coef.it₁₂[iX,iY,iM]
				end
			end
			mul!(view(coef.tmp2,:,:,iM), coef.planfft, view(coef.tmp1,:,:,iM))
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp2[iX,iY,iM] *= coef.Δ[iX,iY,iM]
				end
			end
			ldiv!(view(coef.tmp1,:,:,iM), coef.planfft, view(coef.tmp2,:,:,iM))
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					if iM == 1
						fieldr.e_SXY[1,iX,iY] = coef.tmp1[iX,iY,1] * coef.st₁₂[iX,iY,1]
					else
						coef.tmp1[iX,iY,iM] *= coef.st₁₂[iX,iY,iM]
					end
				end
			end
		end
		@simd for iM in 2:sizeM
			for iY in 1:sizeY
				for iX in 1:sizeX
					fieldr.e_SXY[1,iX,iY] += coef.tmp1[iX,iY,iM]
				end
			end
		end
	else
		Threads.@threads for iM in 1:sizeM
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp1[iX,iY,iM] = fieldi.e_SXY[1,iX,iY] * coef.ir₂₁[iX,iY,iM]
				end
			end
			mul!(view(coef.tmp2,:,:,iM), coef.planfft, view(coef.tmp1,:,:,iM))
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp2[iX,iY,iM] *= coef.Δ[iX,iY,iM]
				end
			end
			ldiv!(view(coef.tmp1,:,:,iM), coef.planfft, view(coef.tmp2,:,:,iM))
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					if iM == 1
						fieldr.e_SXY[1,iX,iY] = coef.tmp1[iX,iY,1] * coef.sr₂₁[iX,iY,1]
					else
						coef.tmp1[iX,iY,iM] *= coef.sr₂₁[iX,iY,iM]
					end
				end
			end
		end
		@simd for iM in 2:sizeM
			for iY in 1:sizeY
				for iX in 1:sizeX
					fieldr.e_SXY[1,iX,iY] += coef.tmp1[iX,iY,iM]
				end
			end
		end
		Threads.@threads for iM in 1:sizeM
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp1[iX,iY,iM] = fieldi.e_SXY[1,iX,iY] * coef.it₂₁[iX,iY,iM]
				end
			end
			mul!(view(coef.tmp2,:,:,iM), coef.planfft, view(coef.tmp1,:,:,iM))
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp2[iX,iY,iM] *= coef.Δ[iX,iY,iM]
				end
			end
			ldiv!(view(coef.tmp1,:,:,iM), coef.planfft, view(coef.tmp2,:,:,iM))
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					if iM == 1
						fieldl.e_SXY[1,iX,iY] = coef.tmp1[iX,iY,iM] * coef.st₂₁[iX,iY,iM]
					else
						coef.tmp1[iX,iY,iM] *= coef.st₂₁[iX,iY,iM]
					end
				end
			end
		end
		@simd for iM in 2:sizeM
			for iY in 1:sizeY
				for iX in 1:sizeX
					fieldl.e_SXY[1,iX,iY] += coef.tmp1[iX,iY,iM]
				end
			end
		end
	end
	if dir(fieldi) > 0
		@simd for i in 1:sizeY*sizeX
			fieldl.e_SXY[i] += fieldi.e_SXY[i] * coef.r₁₂[i]
			fieldr.e_SXY[i] += fieldi.e_SXY[i] * coef.t₁₂[i]
		end
	else
		@simd for i in 1:sizeY*sizeX
			fieldr.e_SXY[i] += fieldi.e_SXY[i] * coef.r₂₁[i]
			fieldl.e_SXY[i] += fieldi.e_SXY[i] * coef.t₂₁[i]
		end
	end
end

function correctscatteringmatrix_referenceframes!(scat::RoughInterfaceConvolutionCoefficient, comp::AbstractOpticalComponent, fieldi::AbstractFieldMonochromatic)
	ref = dir(fieldi) > 0 ? ref1(comp) : ref2(comp)
	(sizeX, sizeY) = (length(fieldi.nsx_X), length(fieldi.nsy_Y))
	if !checkposition(fieldi.ref, ref)
		propM = propagationmatrix(fieldi, ref)
		ind = LinearIndices((sizeX, sizeY))
		if dir(fieldi) > 0
			vec(scat.r₁₂) .*= propM.diag .* propM.diag
			vec(scat.t₁₂) .*= propM.diag
			vec(scat.t₂₁) .*= propM.diag
			@inbounds @simd for iM in 1:size(scat.ir₁₂, 3)
				for iY in 1:sizeY
					for iX in 1:sizeX
						i = ind[iX,iY]
						scat.ir₁₂[iX,iY,iM] *= propM.diag[i]
						scat.sr₁₂[iX,iY,iM] *= propM.diag[i]
						scat.it₁₂[iX,iY,iM] *= propM.diag[i]
						scat.st₂₁[iX,iY,iM] *= propM.diag[i]
					end
				end
			end
		else
			conj!(propM.diag)
			vec(scat.r₂₁) .*= propM.diag .* propM.diag
			vec(scat.t₁₂) .*= propM.diag
			vec(scat.t₂₁) .*= propM.diag
			@inbounds @simd for iM in 1:size(scat.ir₁₂, 3)
				for iY in 1:sizeY
					for iX in 1:sizeX
						i = ind[iX,iY]
						scat.ir₂₁[iX,iY,iM] *= propM.diag[i]
						scat.sr₂₁[iX,iY,iM] *= propM.diag[i]
						scat.st₁₂[iX,iY,iM] *= propM.diag[i]
						scat.it₂₁[iX,iY,iM] *= propM.diag[i]
					end
				end
			end
		end
	end
end

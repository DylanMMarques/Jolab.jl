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
	tmp1::Y
	tmp2::Y
end

function lightinteraction!(fieldl::L, fieldr::R, coef::RoughInterfaceConvolutionCoefficient{T,L,R}, fieldi::Union{L,R}) where {T,L,R}
	samedefinitions(fieldi, fieldi.dir > 0 ? fieldl : fieldr) || tobedone()
	(sizeX, sizeY) = size(fieldi.e_SXY)[2:3]
	@inbounds for iM in 1:size(coef.ir₁₂, 3)
		if fieldi.dir > 0
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp1[iX,iY] = fieldi.e_SXY[1,iX,iY] * coef.ir₁₂[iX,iY,iM]
				end
			end
			mul!(coef.tmp2, coef.planfft, coef.tmp1)
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp2[iX,iY] *= coef.Δ[iX,iY,iM]
				end
			end
			ldiv!(coef.tmp1, coef.planfft, coef.tmp2)
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					if iM > 1
						fieldl.e_SXY[1,iX,iY] += coef.tmp1[iX,iY] * coef.sr₁₂[iX,iY,iM]
					else
						fieldl.e_SXY[1,iX,iY] = coef.tmp1[iX,iY] * coef.sr₁₂[iX,iY,iM]
					end
				end
			end
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp1[iX,iY] = fieldi.e_SXY[1,iX,iY] * coef.it₁₂[iX,iY,iM]
				end
			end
			mul!(coef.tmp2, coef.planfft, coef.tmp1)
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp2[iX,iY] *= coef.Δ[iX,iY,iM]
				end
			end
			ldiv!(coef.tmp1, coef.planfft, coef.tmp2)
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					if iM > 1
						fieldr.e_SXY[1,iX,iY] += coef.tmp1[iX,iY] * coef.st₁₂[iX,iY,iM]
					else
						fieldr.e_SXY[1,iX,iY] = coef.tmp1[iX,iY] * coef.st₁₂[iX,iY,iM]
					end
				end
			end
		else
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp1[iX,iY] = fieldi.e_SXY[1,iX,iY] * coef.ir₂₁[iX,iY,iM]
				end
			end
			mul!(coef.tmp2, coef.planfft, coef.tmp1)
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp2[iX,iY] *= coef.Δ[iX,iY,iM]
				end
			end
			ldiv!(coef.tmp1, coef.planfft, coef.tmp2)
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					if iM > 1
						fieldr.e_SXY[1,iX,iY] += coef.tmp1[iX,iY] * coef.sr₂₁[iX,iY,iM]
					else
						fieldr.e_SXY[1,iX,iY] = coef.tmp1[iX,iY] * coef.sr₂₁[iX,iY,iM]
					end
				end
			end
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp1[iX,iY] = fieldi.e_SXY[1,iX,iY] * coef.it₂₁[iX,iY,iM]
				end
			end
			mul!(coef.tmp2, coef.planfft, coef.tmp1)
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					coef.tmp2[iX,iY] *= coef.Δ[iX,iY,iM]
				end
			end
			ldiv!(coef.tmp1, coef.planfft, coef.tmp2)
			@simd for iY in 1:sizeY
				for iX in 1:sizeX
					if iM > 1
						fieldl.e_SXY[1,iX,iY] += coef.tmp1[iX,iY] * coef.st₂₁[iX,iY,iM]
					else
						fieldl.e_SXY[1,iX,iY] = coef.tmp1[iX,iY] * coef.st₂₁[iX,iY,iM]
					end
				end
			end
		end
	end
	if fieldi.dir > 0
		@inbounds @simd for i in 1:sizeY*sizeX
			fieldl.e_SXY[i] += fieldi.e_SXY[i] * coef.r₁₂[i]
			fieldr.e_SXY[i] += fieldi.e_SXY[i] * coef.t₁₂[i]
		end
	else
		@inbounds @simd for i in 1:sizeY*sizeX
			fieldr.e_SXY[i] += fieldi.e_SXY[i] * coef.r₂₁[i]
			fieldl.e_SXY[i] += fieldi.e_SXY[i] * coef.t₂₁[i]
		end
	end
end

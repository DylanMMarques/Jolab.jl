struct RoughInterfaceConvolutionCoefficient{T,L,R,X<:AbstractArray{Complex{T}},F <: AbstractFFTs.Plan}  <: AbstractCoefficient{T,L,R}
	r₁₂::X
	t₁₂::X
	r₂₁::X
	t₂₁::X
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
	tmp::X
end

function lightinteraction!(fieldl::L, fieldr::R, coef::RoughInterfaceConvolutionCoefficient{T,L,R}, fieldi::Union{L,R}) where {T,L,R}
	samedefinitions(fieldi, fieldi.dir > 0 ? fieldl : fieldr) || tobedone()
	fieldle_SXY = reshape(fieldl.e_SXY, size(fieldl.e_SXY)[2:3]) # to apply 2D FFT
	fieldre_SXY = reshape(fieldr.e_SXY, size(fieldr.e_SXY)[2:3]) # to apply 2D FFT
	if fieldi.dir > 0
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			fieldle_SXY[i] = fieldi.e_SXY[i] * coef.ir₁₂[i]
			fieldre_SXY[i] = fieldi.e_SXY[i] * coef.it₁₂[i]
		end
		mul!(coef.tmp, coef.planfft, fieldle_SXY)
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			coef.tmp[i] *= coef.Δ[i]
		end
		ldiv!(fieldle_SXY, coef.planfft, coef.tmp)
		mul!(coef.tmp, coef.planfft, fieldre_SXY)
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			coef.tmp[i] *= coef.Δ[i]
		end
		ldiv!(fieldre_SXY, coef.planfft, coef.tmp)
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			fieldle_SXY[i] *= coef.sr₁₂[i]
			fieldle_SXY[i] += fieldi.e_SXY[i] * coef.r₁₂[i]

			fieldre_SXY[i] *= coef.st₁₂[i]
			fieldre_SXY[i] += fieldi.e_SXY[i] * coef.t₁₂[i]
		end
	else
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			fieldle_SXY[i] = fieldi.e_SXY[i] * coef.it₂₁[i]
			fieldre_SXY[i] = fieldi.e_SXY[i] * coef.ir₂₁[i]
		end
		mul!(coef.tmp, coef.planfft, fieldle_SXY)
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			coef.tmp[i] *= coef.Δ[i]
		end
		ldiv!(fieldle_SXY, coef.planfft, coef.tmp)
		mul!(coef.tmp, coef.planfft, fieldre_SXY)
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			coef.tmp[i] *= coef.Δ[i]
		end
		ldiv!(fieldre_SXY, coef.planfft, coef.tmp)
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			fieldle_SXY[i] *= coef.st₂₁[i]
			fieldle_SXY[i] += fieldi.e_SXY[i] * coef.t₂₁[i]

			fieldre_SXY[i] *= coef.sr₂₁[i]
			fieldre_SXY[i] += fieldi.e_SXY[i] * coef.r₂₁[i]
		end
	end
end

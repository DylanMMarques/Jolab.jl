struct RoughInterfaceConvolutionCoefficient{T,L,R,X<:AbstractArray{Complex{T}},F <: AbstractFFTs.Plan}  <: AbstractCoefficient{T,L,R}
	r₁₂::X
	t₁₂::X
	r₂₁::X
	t₂₁::X
	it₁₂::X
	st₁₂::X
	Δt₁₂::X
	ir₂₁::X
	sr₂₁::X
	Δr₂₁::X
	it₂₁::X
	st₂₁::X
	Δt₂₁::X
	fieldl::L
	fieldr::R
	planfft::F
	tmp::X
end

function lightinteraction!(fieldl::L, fieldr::R, coef::RoughInterfaceConvolutionCoefficient{T,L,R}, fieldi::Union{L,R}) where {T,L,R}
	fieldle_SXY = reshape(fieldl.e_SXY, size(fieldl.e_SXY)[2:3]) # to apply 2D FFT
	fieldre_SXY = reshape(fieldr.e_SXY, size(fieldr.e_SXY)[2:3]) # to apply 2D FFT
	if fieldi.dir > 0
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			fieldle_SXY[i] = fieldi.e_SXY[i] * coef.ir₁₂[i]
			fieldre_SXY[i] = fieldi.e_SXY[i] * coef.it₁₂[i]
		end
		mul!(coef.tmp, coef.planfft, fieldle_SXY)
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			coef.tmp[i] *= coef.Δr₁₂[i]
		end
		ldiv!(fieldle_SXY, fftPlan, coef.tmp_XY)
		mul!(coef.tmp, coef.planfft, fieldr.e_SXY)
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			coef.tmp[i] *= coef.Δt₁₂[i]
		end
		ldiv!(fieldre_SXY, fftPlan, tmp2_XY)
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
			coef.tmp[i] *= coef.Δt₂₁[i]
		end
		ldiv!(fieldle_SXY, fftPlan, coef.tmp_XY)
		mul!(coef.tmp, coef.planfft, fieldr.e_SXY)
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			coef.tmp[i] *= coef.Δr₂₁[i]
		end
		ldiv!(fieldre_SXY, fftPlan, tmp2_XY)
		@inbounds @simd for i in eachindex(fieldi.e_SXY)
			fieldle_SXY[i] *= coef.st₂₁[i]
			fieldle_SXY[i] += fieldi.e_SXY[i] * coef.t₂₁[i]

			fieldre_SXY[i] *= coef.sr₂₁[i]
			fieldre_SXY[i] += fieldi.e_SXY[i] * coef.r₂₁[i]
		end
	end
end

# function lightinteraction!(ebackward_SXY::AbstractArray{Complex{T}, 3}, eforward_SXY::AbstractArray{Complex{T}, 3}, coef::ScatteringConvolutionCoefficientScalar, nsx_X::AbstractRange{<:Number}, nsy_Y::AbstractRange{<:Number}, ei_SXY::AbstractArray{<:Number, 3}, dir::Integer, fftPlan=plan_fft(reshape(ebackward_SXY, size(ebackward_SXY)[2:3]))::AbstractFFTs.Plan{Complex{T}}, tmp_SXY=copy(ebackward_SXY)::AbstractArray{Complex{T},3}) where {T<:Real}
# 	(sizeS, sizeX, sizeY) = size(ebackward_SXY)
# 	size(eforward_SXY, 1) == sizeS == 1 || error("Field must be scallar to apply this theory")
#
# 	ebackward_XY = reshape(ebackward_SXY, sizeX, sizeY)
# 	eforward_XY = reshape(eforward_XY, sizeX, sizeY)
# 	tmp_XY = reshape(tmp_SXY, sizeX, sizeY)
# 	#replace coef by the respective matrices.
# 	#create a new for not recursive optimization
# 	x_X = FFTW.fftfreq(sizeX, nsx_X[2] - nsx_X[1])
# 	y_Y = FFTW.fftfreq(sizeY, nsy_Y[2] - nsy_Y[1])
# 	if dir > 0
# 	else
# 		@inbounds @simd for iY in eachindex(nsy_Y)
# 			for iX in eachindex(nsx_X)
# 				nsr = √(nsx_X[iX]^2 + nsy_Y[iY]^2)
# 				ebackward_XY[iX,iY] = ei_SXY[1,iX,iY] * coef.it₂₁(nsr)
# 				eforward_XY[iX,iY] = ei_SXY[1,iX,iY] * coef.ir₂₁(nsr)
# 			end
# 		end
# 		mul!(tmp_XY, fftPlan, ebackward_XY)
# 		@inbounds @simd for iY in eachindex(nsy_Y)
# 			for iX in eachindex(nsx_X)
# 				tmp_XY[iX,iY] *= coef.Δt₂₁(x_X[iX], y_Y[iY])
# 			end
# 		end
# 		ldiv!(ebackward_XY, fftPlan, tmp_XY)
# 		mul!(tmp_XY, fftPlan, eforward_XY)
# 		@inbounds @simd for iY in eachindex(nsy_Y)
# 			for iX in eachindex(nsx_X)
# 				tmp_XY[iX,iY] *= coef.Δr₂₁(x_X[iX], y_Y[iY])
# 			end
# 		end
# 		ldiv!(eforward_XY, fftPlan, tmp_XY)
# 		@inbounds @simd for iY in eachindex(nsy_Y)
# 			for iX in eachindex(nsx_X)
# 				nsr = √(nsx_X[iX]^2 + nsy_Y[iY]^2)
# 				eforward_XY[iX,iY] *= coef.sr₂₁(nsr)
# 				ebackward_XY[iX,iY] *= coef.st₂₁(nsr)
# 			end
# 		end
# 	end
# end

# function lightinteraction!(ebackward_SXY::AbstractArray{Complex{T}, 3}, eforward_SXY::AbstractArray{Complex{T}, 3}, coef::ScatteringConvolutionCoefficientScalar{A,X1,X2}, nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number, 2}, ei_SXY::AbstractArray{<:Number, 3}, dir::Integer, fftPlan=plan_fft(reshape(ebackward_SXY, size(ebackward_SXY)[2:3]))::AbstractFFTs.Plan{Complex{T}}, tmp_SXY=copy(ebackward_SXY)::AbstractArray{Complex{T},3}) where {T<:Real, A, X1<:AbstractArray, X2<:AbstractArray}
# 	(sizeS, sizeX, sizeY) = size(ebackward_SXY)
# 	size(eforward_SXY, 1) == sizeS == 1 || error("Field must be scallar to apply this theory")
#
# 	ebackward_XY = reshape(ebackward_SXY, sizeX, sizeY)
# 	eforward_XY = reshape(eforward_SXY, sizeX, sizeY)
# 	tmp_XY = reshape(tmp_SXY, sizeX, sizeY)
# 	iterator = eachindex(ebackward_XY)
# 	if dir > 0
# 		@inbounds @simd for i in iterator
# 			ebackward_XY[i] = ei_SXY[i] * coef.ir₁₂[i]
# 			eforward_XY[i] = ei_SXY[i] * coef.it₁₂[i]
# 		end
# 		ldiv!(tmp_XY, fftPlan, ebackward_XY)
# 		@inbounds @simd for i in iterator
# 			tmp_XY[i] *= coef.Δr₁₂[i]
# 		end
# 		mul!(ebackward_XY, fftPlan, tmp_XY)
#
# 		ldiv!(tmp_XY, fftPlan, eforward_XY)
# 		@inbounds @simd for i in iterator
# 		 tmp_XY[i] *= coef.Δt₁₂[i]
# 		end
# 		mul!(eforward_XY, fftPlan, tmp_XY)
# 		@inbounds @simd for i in iterator
# 			ebackward_XY[i] *= coef.sr₁₂[i]
# 			eforward_XY[i] *= coef.st₁₂[i]
# 		end
# 	else
# 		@inbounds @simd for i in iterator
# 			ebackward_XY[i] = ei_SXY[i] * coef.it₂₁[i]
# 			eforward_XY[i] = ei_SXY[i] * coef.ir₂₁[i]
# 		end
# 		ldiv!(tmp_XY, fftPlan, ebackward_XY)
# 		@inbounds @simd for i in iterator
# 			tmp_XY[i] *= coef.Δt₂₁[i]
# 		end
# 		mul!(ebackward_XY, fftPlan, tmp_XY)
#
# 		ldiv!(tmp_XY, fftPlan, eforward_XY)
# 		@inbounds @simd for i in iterator
# 			tmp_XY[i] *= coef.Δr₂₁[i]
# 		end
# 		mul!(eforward_XY, fftPlan, tmp_XY)
# 		@inbounds @simd for i in iterator
# 			ebackward_XY[i] *= coef.st₂₁[i]
# 			eforward_XY[i] *= coef.sr₂₁[i]
# 		end
# 	end
# end
#
# function lightinteraction!(ebackward_SXY::AbstractArray{Complex{T}, 3}, eforward_SXY::AbstractArray{Complex{T}, 3}, coef::PropagationScatteringConvolutionCoefficientScalar, nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number, 2}, ei_SXY::AbstractArray{<:Number, 3}, dir::Integer, tmpₚ_3 = MVector{3,Complex{T}}(undef)::AbstractVector{Complex{T}}, fftPlan=plan_fft(reshape(ebackward_SXY, size(ebackward_SXY)[2:3]))::AbstractFFTs.Plan{Complex{T}}, tmp_SXY=copy(ebackward_SXY)::AbstractArray{Complex{T},3}, tmp2_SXY=copy(eforward_SXY)::AbstractArray{Complex{T},3}, tmpₛ_3 = MVector{3, Complex{T}}(undef)::AbstractVector{Complex{T}}) where T<:Real
# 	lightinteraction!(ebackward_SXY, eforward_SXY, coef.scatConvCoef, nsx_XY, nsy_XY, ei_SXY, dir, fftPlan, tmp_SXY);
# 	ebackward_SXY .*= 4π^2
# 	eforward_SXY .*= 4π^2
# 	lightinteraction!(tmp_SXY, tmp2_SXY, coef.propCoef, nsx_XY, nsy_XY, ei_SXY, dir, tmpₛ_3, tmpₚ_3)
# 	ebackward_SXY .+= tmp_SXY
# 	eforward_SXY .+= tmp2_SXY
# end

# function lightinteraction_recursivegridded!(ebackward_SXY::AbstractArray{Complex{T}, 3}, eforward_SXY::AbstractArray{Complex{T}, 3}, coefs::AbstractVector{<:AbstractCoefficient}, nsx_XY::AbstractArray{<:Number, 2}, nsy_XY::AbstractArray{<:Number,2}, ei_SXY::AbstractArray{<:Number, 3}, λ::Real, dir::Integer, thresold::Real) where T<:Real
#
# 	sizeL = length(coefs) + 1;
# 	(sizeS, sizeX, sizeY) = size(ei_SXY);
# 	(sizeX, sizeY) == size(nsx_XY) == size(nsy_XY) == size(ebackward_SXY)[2:3] == size(eforward_SXY)[2:3] || error("Wrong sizes")
# 	sizeS == size(ebackward_SXY, 1) == size(eforward_SXY, 1) || error("wrong sizes")
#
# 	eForwardLayers_SXYL = zeros(Complex{T}, sizeS, sizeX, sizeY, sizeL)
# 	eBackwardLayers_SXYL = zeros(Complex{T}, sizeS, sizeX, sizeY, sizeL)
# 	eSaveBackward_SXY = Array{Complex{T},3}(undef, sizeS, sizeX, sizeY)
# 	eSaveForward_SXY = Array{Complex{T},3}(undef, sizeS, sizeX, sizeY)
# 	nsz_XYL = Array{Complex{T},3}(undef, sizeX, sizeY, sizeL)
# 	intForward_L = @MVector zeros(T, sizeL)
# 	intBackward_L = @MVector zeros(T, sizeL)
# 	tmpnₛ_3 = MVector{3, Complex{T}}(undef)
# 	tmpnₚ_3 = MVector{3, Complex{T}}(undef)
#
# 	needFFT = false
# 	for i in eachindex(coefs)
# 		view(nsz_XYL,:,:,i) .= .√(coefs[i].n₁^2 .- nsx_XY.^2 .- nsy_XY.^2);
# 		if (coefs[i] isa AbstractScatteringConvolutionCoefficient || coefs[i] isa AbstractPropagationScatteringConvolutionCoefficient)
# 			needFFT = true
# 		end
# 	end
# 	if needFFT
# 		fftPlan = plan_fft(view(ebackward_SXY,1,:,:)) # Pre calculates the fourrier transform
# 		inv(fftPlan) # Pre calculates the inverse fourrier transform
# 		tmp_SXY = Array{Complex{T}, 3}(undef, sizeS, sizeX, sizeY) # Initialize temporary arrays to avoid allocations inside loop
# 		tmp2_SXY = Array{Complex{T}, 3}(undef, sizeS, sizeX, sizeY) # Initialize temporary arrays to avoid allocations inside loop
# 	end
#
# 	if dir > 0
# 		view(eForwardLayers_SXYL,:,:,:,1) .= ei_SXY;
# 		intForward_L[1] = intensity(ei_SXY[1,:,:] .* .√(nsz_XYL[:,:,1]))
# 		rtol = intForward_L[1] * thresold^2
# 	else
# 		view(eBackwardLayers_SXYL,:,:,:,sizeL) .= ei_SXY;
# 		intBackward_L[sizeL] = intensity(ei_SXY[1,:,:] .* √(nsz_XYL[:,:,sizeL]))
# 		rtol = iBackward_L[sizeL] * thresold^2
# 	end
# 	view(nsz_XYL,:,:,sizeL) .= .√(coefs[sizeL-1].n₂^2 .- nsx_XY.^2 .- nsy_XY.^2)
#
# 	i = 1
# 	mls = 1
# 	isForward = true
#
# 	while true
# 		# Select the next field to consider
# 		(maxBack, argBack) = findmax(view(intBackward_L, 2:sizeL)) # Removes the first as it is the field going out
# 		(maxForw, argForw) = findmax(view(intForward_L, 1:sizeL-1)) # Removes the last as it is the field going out
# 		# @show isForward
# 		# @show intBackward_L
# 		# @show intForward_L
# 		# @show (sum(intBackward_L) + sum(intForward_L)) / rtol * thresold^2
# 		# readline()
# 		# (maxForw + maxBack) / rtol * thresold^2 > 10 && error("Energy is being created.")
# 		isForward = maxForw > maxBack
# 		if isForward
# 			maxForw < rtol && break
# 			mls = argForw
# 			e_SXY = @view eForwardLayers_SXYL[:,:,:,mls]
# 			nsz_XY = @view nsz_XYL[:,:,mls]
# 		else
# 			maxBack < rtol && break
# 			mls = argBack + 1 # need to add +1 because length start from 2
# 			e_SXY = @view eBackwardLayers_SXYL[:,:,:,mls]
# 			nsz_XY = @view nsz_XYL[:,:,mls]
# 		end
#
# 		if isForward
# 			(mls > 1) && changereferenceframe!(nsx_XY, nsy_XY, nsz_XY, e_SXY, λ, coefs[mls-1].ref₂, coefs[mls].ref₁)
# 			if coefs[mls] isa AbstractScatteringConvolutionCoefficient
# 				lightinteraction!(eSaveBackward_SXY, eSaveForward_SXY, coefs[mls], nsx_XY, nsy_XY, e_SXY, 1, fftPlan, tmp_SXY, tmp2_SXY)
# 			elseif coefs[mls] isa AbstractPropagationCoefficient
# 				lightinteraction!(eSaveBackward_SXY, eSaveForward_SXY, coefs[mls], nsx_XY, nsy_XY, e_SXY, 1, tmpnₚ_3, tmpnₛ_3)
# 			elseif coefs[mls] isa  AbstractPropagationScatteringConvolutionCoefficient
# 				lightinteraction!(eSaveBackward_SXY, eSaveForward_SXY, coefs[mls], nsx_XY, nsy_XY, e_SXY, 1, tmpnₚ_3, fftPlan, tmp_SXY, tmp2_SXY)
# 			else
# 				error("Coefficient type not known")
# 			end
# 			view(eBackwardLayers_SXYL,:,:,:,mls) .+= eSaveBackward_SXY
# 			intBackward_L[mls] = intensity(view(eBackwardLayers_SXYL,1,:,:,mls) .* .√(view(nsz_XYL,:,:,mls)))
# 			view(eForwardLayers_SXYL,:,:,:,mls+1) .+= eSaveForward_SXY
# 			intForward_L[mls+1] = intensity(view(eForwardLayers_SXYL,1,:,:,mls+1) .* .√(view(nsz_XYL,:,:,mls+1)))
# 			intForward_L[mls] = zero(T)
# 		else
# 			nsz_XY .*= -1
# 			(mls < sizeL) && changereferenceframe!(nsx_XY, nsy_XY, nsz_XY, e_SXY, λ, coefs[mls].ref₁, coefs[mls-1].ref₂)
# 			nsz_XY .*= -1
# 			if coefs[mls-1] isa AbstractScatteringConvolutionCoefficient
# 				lightinteraction!(eSaveBackward_SXY, eSaveForward_SXY, coefs[mls-1], nsx_XY, nsy_XY, e_SXY, -1, fftPlan, tmp_SXY, tmp2_SXY)
# 			elseif coefs[mls-1] isa AbstractPropagationCoefficient
# 				lightinteraction!(eSaveBackward_SXY, eSaveForward_SXY, coefs[mls-1], nsx_XY, nsy_XY, e_SXY, -1, tmpnₚ_3, tmpnₛ_3)
# 			elseif coefs[mls-1] isa  AbstractPropagationScatteringConvolutionCoefficient
# 				lightinteraction!(eSaveBackward_SXY, eSaveForward_SXY, coefs[mls-1], nsx_XY, nsy_XY, e_SXY, -1, tmpnₚ_3, fftPlan, tmp_SXY, tmp2_SXY)
# 			else
# 				error("Coefficient type not known")
# 			end
# 			vec(view(eBackwardLayers_SXYL,:,:,:,mls-1)) .+= vec(eSaveBackward_SXY)
# 			intBackward_L[mls-1] = intensity(view(eBackwardLayers_SXYL,1,:,:,mls-1) .* .√(view(nsz_XYL,:,:,mls-1)))
# 			vec(view(eForwardLayers_SXYL,:,:,:,mls)) .+= vec(eSaveForward_SXY)
# 			intForward_L[mls] = intensity(view(eForwardLayers_SXYL,1,:,:,mls) .* .√(view(nsz_XYL,:,:,mls)))
# 			intBackward_L[mls] = zero(T)
# 		end
# 		e_SXY .= zero(Complex{T})
# 		i = i + 1;
# 		i > 5000 && break
# 	end
# 	ebackward_SXY .= view(eBackwardLayers_SXYL,:,:,:,1)
# 	eforward_SXY .= view(eForwardLayers_SXYL,:,:,:,sizeL)
# 	return (ebackward_SXY, eforward_SXY)
# end
#
# function lightinteraction_recursivegridded(coefs::AbstractVector{<:AbstractCoefficient{A}}, angspe::FieldAngularSpectrum{T,X}; thresold = 1E-4::Real) where {T,X<:AbstractRange, A}
# 	checkapplicability(coefs, angspe);
# 	(sizeS, sizeX, sizeY) = size(angspe.e_SXY);
#
# 	ebackward_SCD = Array{Complex{T},3}(undef, sizeS, sizeX, sizeY);
# 	eforward_SCD = Array{Complex{T},3}(undef, sizeS, sizeX, sizeY);
# 	ei_SXY = copy(angspe.e_SXY);
#
# 	nsx_XY = repeat(angspe.nsx_X, 1, sizeY);
# 	nsy_XY = repeat(angspe.nsy_Y', sizeX, 1);
# 	nsz_XY = angspe.dir * .√(complex.(angspe.n^2 .- nsx_XY.^2 .- nsy_XY.^2))
#
# 	refcoef = angspe.dir > 0 ? coefs[1].ref₁ : coefs[end].ref₂
# 	checkorientation(refcoef, angspe.ref) || tobedone()
# 	changereferenceframe!(nsx_XY, nsy_XY, nsz_XY, ei_SXY, angspe.λ, angspe.ref, refcoef)
#
# 	coefs_matrix = [coefficient_matrixform(coefi, nsx_XY, nsy_XY) for coefi in coefs]
#
# 	lightinteraction_recursivegridded!(ebackward_SCD, eforward_SCD, coefs_matrix, nsx_XY, nsy_XY, ei_SXY, angspe.λ, angspe.dir, thresold);
#
# 	angspeforward = FieldAngularSpectrum{T}(angspe.nsx_X, angspe.nsy_Y, eforward_SCD, angspe.λ, coefs[end].n₂, 1, coefs[end].ref₂)
# 	angspebackward = FieldAngularSpectrum{T}(angspe.nsx_X, angspe.nsy_Y, ebackward_SCD, angspe.λ, coefs[1].n₁, -1, coefs[1].ref₁)
# 	return (angspebackward, angspeforward)
# end
#
# function lightinteraction_recursivegridded(components::AbstractVector{<:AbstractOpticalComponent{T}}, angspe::AbstractFieldAngularSpectrum; thresold = 1E-4) where T
# 	vectorial = (size(angspe.e_SXY, 1) == 3)
# 	vectorial && tobedone()
#
# 	coefs = [coefficientscallar(compi, angspe.λ) for compi in components]
# 	return lightinteraction_recursivegridded(coefs, angspe, thresold = thresold);
# end

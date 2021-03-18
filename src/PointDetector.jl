struct PointDetector{T} <: AbstractOpticalComponent{T}
	ref::ReferenceFrame{T}
end
PointDetector(ref) = PointDetector{Float64}(ref)

function signal(pointdet::PointDetector, angspe::AbstractFieldAngularSpectrum{T}) where T
	# Return the signal measured on the point detector
	# i = ¦∫∫ E k² dsx dsy¦²
	angspe_new = changereferenceframe(angspe, pointdet.ref)
	aux = zeros(T,1)
	fourier = FourierTransform{T, Vector{T}, ReferenceFrame{T}}(aux, aux, aux, aux, pointdet)
	(fieldl, fieldr) = lightinteraction(fourier, angspe)
	return dir(angspe) > 0 ? abs2(fieldr.e_SXY[1]) : abs2(fieldl.e_SXY[1])
end

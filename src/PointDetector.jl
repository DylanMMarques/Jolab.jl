struct PointDetector{T} <: AbstractOpticalComponent{T}
	ref::ReferenceFrame{T}
end
PointDetector(ref) = PointDetector{Float64}(ref)

function signal(pointdet::PointDetector, angspe::AbstractFieldAngularSpectrum)
	# Return the signal measured on the point detector
	# i = ¦∫∫ E k² dsx dsy¦²
	e_SXY = angspeto3Dspace(angspe, [pointdet.ref.x], [pointdet.ref.y], [pointdet.ref.z]);
	length(e_SXY) == 1 ? i = abs2.(e_SXY) : i = dotdim(e_SXY, e_SXY, 1);
	return i[1];
end

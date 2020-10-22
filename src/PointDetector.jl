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

function coefficient_general(pointdet::PointDetector{T}, fieldi::FieldAngularSpectrum{T}) where T
	tobedone()
	sizeXY = length(fieldi.nsx_X) * length(fieldi.nsy_Y)
	r12 = zeros(Complex{T}, sizeXY, sizeXY)
	r21 = zeros(Complex{T}, sizeXY, sizeXY)

	if dir(fieldi) > 0
		t12 = Matrix{Complex{T}}(undef, 1, sizeXY)
		t21 = zeros(Complex{T}, 1, sizeXY)
	else
		t21 = Matrix{Complex{T}}(undef, 1, sizeXY)
		t12 = zeros(Complex{T}, 1, sizeXY)
	end

	k = 2π / fieldi.λ
	Δkx = k * ΔIntegrationTrap(fieldi.nsx_X)
	Δky = k * ΔIntegrationTrap(fieldi.nsy_Y)

	xdet = pointdet.ref.x - fieldi.ref.x;
	ydet = pointdet.ref.y - fieldi.ref.y;
	zdet = pointdet.ref.z - fieldi.ref.z;

	(xdet, ydet, zdet) = rotatecoordinatesfrom(xdet, ydet, zdet, fieldi.ref.θ, fieldi.ref.ϕ);
	i = 1
	@inbounds for iY in eachindex(fieldi.nsy_Y)
		for iX in eachindex(fieldi.nsx_X)
			nsz = √(fieldi.n^2 - (fieldi.nsx_X[iX]^2 + fieldi.nsy_Y[iY]^2))
			if dir(fieldi) > 0
				t12[i] = exp(im * k * (xdet * fieldi.nsx_X[iX] + ydet * fieldi.nsy_Y[iY] + zdet * nsz)) * Δkx[iX] * Δky[iY]
			else
				t21[i] = exp(im * k * (xdet * fieldi.nsx_X[iX] + ydet * fieldi.nsy_Y[iY] - zdet * nsz)) * Δkx[iX] * Δky[iY]
			end
			i +=1
		end
	end

	if dir(fieldi) > 0
		fieldl = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, -1, fieldi.ref)
		fieldr = FieldAngularSpectrum{T}(zeros(T,1), zeros(T, 1), zeros(Complex{T}, 1,1,1), fieldi.λ, fieldi.n, 1, pointdet.ref)
	else
		fieldr = FieldAngularSpectrum{T}(zeros(T,1), zeros(T, 1), zeros(Complex{T}, 1,1,1), fieldi.λ, fieldi.n, -1, pointdet.ref)
		fieldr = FieldAngularSpectrum{T}(copy(fieldi.nsx_X), copy(fieldi.nsy_Y), fieldi.e_SXY, fieldi.λ, fieldi.n, 1, fieldi.ref)
	end
	return ScatteringMatrix{T, typeof(fieldl), typeof(fieldr), Matrix{Complex{T}}}(r12, t12, r21, t21, fieldl, fieldr)
end

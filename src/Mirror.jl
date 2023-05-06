struct Mirror{T} <: AbstractPropagationComponent{T}
    R::JolabFunction{T}
    n₁::Material{T}
    n₂::Material{T}
    ref::ReferenceFrame{T}
	Mirror{T}(R, n1, n2, ref) where T = new{T}(R, n1, n2, ref)
end

errorToDo() = error("Update in next versions")

Mirror(R, n₁, n₂, ref) = Mirror{Float64}(R,n₁,n₂,ref)
R(mirror::Mirror, λ) = mirror.R(λ)

n1(mirror::Mirror, λ) = n(mirror.n₁, λ)
n2(mirror::Mirror, λ) = n(mirror.n₂, λ)

function coefficient_general(mirror::Mirror, field::AbstractFieldSpace)
	checkapplicability(mirror, field)

	n1_val, n2_val = √n1(mirror, field.λ), √n2(mirror, field.λ)
    tλ = complex(√(1 - R(mirror, field.λ)));
    rλ = complex(√R(mirror, field.λ));
	scat = get_scatteringmatrixtype(mirror, field)
	scat.r₁₂.diag .= rλ
	scat.t₁₂.diag .= tλ * n1_val / n2_val
	scat.r₂₁.diag .= -rλ
	scat.t₂₁.diag .= tλ * n2_val / n1_val

	return scat
end

function coefficient_general(mirror::Mirror, field::FieldAngularSpectrumScalar)
	checkapplicability(mirror, field)

	n1_val, n2_val = n1(mirror, field.λ), n2(mirror, field.λ)

    tλ = complex(√(1 - R(mirror, field.λ)));
    rλ = complex(√R(mirror, field.λ));
	scat = get_scatteringmatrixtype(mirror, field)

	scat.r₁₂.diag .= rλ
	scat.r₂₁.diag .= -rλ
	cart = CartesianIndices(field)
	@inbounds @simd for i in eachindex(scat.t₁₂.diag)
		nsz1 = √nsz(n1_val, field.nsx_X[cart[i][2]], field.nsy_Y[cart[i][3]])
		nsz2 = √nsz(n2_val, field.nsx_X[cart[i][2]], field.nsy_Y[cart[i][3]])
		scat.t₁₂.diag[i] = tλ * nsz1 / nsz2
		scat.t₂₁.diag[i] = tλ * nsz2 / nsz1
	end
	correctscatteringmatrix_referenceframes!(scat, mirror, field)
	return scat
end

function coefficient_general(mirror::Mirror, field::FieldAngularSpectrumScalarRadialSymmetric)
	checkapplicability(mirror, field)

	n1_val, n2_val = n1(mirror, field.λ), n2(mirror, field.λ)

    tλ = complex(√(1 - R(mirror, field.λ)));
    rλ = complex(√R(mirror, field.λ));
	scat = get_scatteringmatrixtype(mirror, field)

	scat.r₁₂.diag .= rλ
	scat.r₂₁.diag .= -rλ
	@inbounds @simd for i in eachindex(scat.t₁₂.diag)
		nsz1 = √nsz(n1_val, field.nsr_R[i], 0)
		nsz2 = √nsz(n2_val, field.nsr_R[i], 0)
		scat.t₁₂.diag[i] = tλ * nsz1 / nsz2
		scat.t₂₁.diag[i] = tλ * nsz2 / nsz1
	end
	correctscatteringmatrix_referenceframes!(scat, mirror, field)
	return scat
end

function checkapplicability(mirror::Mirror, field::AbstractFieldAngularSpectrum)
	isapprox(field.n, dir(field) > 0 ? n1(mirror, field.λ) : n2(mirror, field.λ), atol = @tol) || error("Field medium and mirror are different")
	checkorientation(field.ref, mirror.ref) || errorToDo()
end

function checkapplicability(mirror::Mirror, field::AbstractFieldSpace)
	isapprox(field.n, dir(field) > 0 ? n1(mirror, field.λ) : n2(mirror, field.λ), atol = @tol) || error("Field medium and mirror are different")
	checkinplane(field.ref, mirror.ref) || error("cannot do this")
end

function get_scatteringmatrixtype(mirror::Mirror, field::FieldAngularSpectrumScalar{T,D,X,Y}) where {T,D,X,Y}
	sizeXY = length(field.e_SXY)
	r12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))

	(fieldl, fieldr) = getfields_lr(mirror, field)

	return ScatteringMatrix{T, FieldAngularSpectrumScalar{T,-1,X,Y}, FieldAngularSpectrumScalar{T,1,X,Y}, Diagonal{Complex{T},Vector{Complex{T}}}, Diagonal{Complex{T},Vector{Complex{T}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

function get_scatteringmatrixtype(mirror::Mirror, field::FieldAngularSpectrumScalarRadialSymmetric{T,D,X,Y}) where {T,D,X,Y}
	sizeXY = length(field.e_SXY)
	r12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))

	(fieldl, fieldr) = getfields_lr(mirror, field)

	return ScatteringMatrix{T, FieldAngularSpectrumScalarRadialSymmetric{T,-1,X,Y}, FieldAngularSpectrumScalarRadialSymmetric{T,1,X,Y}, Diagonal{Complex{T},Vector{Complex{T}}}, Diagonal{Complex{T},Vector{Complex{T}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

function get_scatteringmatrixtype(mirror::Mirror, field::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y}
	sizeXY = length(field.e_SXY)
	r12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))

	(fieldl, fieldr) = getfields_lr(mirror, field)

	return ScatteringMatrix{T, FieldSpaceScalar{T,-1,X,Y}, FieldSpaceScalar{T,1,X,Y}, Diagonal{Complex{T},Vector{Complex{T}}}, Diagonal{Complex{T},Vector{Complex{T}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

function get_scatteringmatrixtype(mirror::Mirror, field::FieldSpaceScalarRadialSymmetric{T,D,X,Y}) where {T,D,X,Y}
	sizeXY = length(field.e_SXY)
	r12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t12 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	r21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))
	t21 = Diagonal(Vector{Complex{T}}(undef, sizeXY))

	(fieldl, fieldr) = getfields_lr(mirror, field)

	return ScatteringMatrix{T, FieldSpaceScalarRadialSymmetric{T,-1,X,Y}, FieldSpaceScalarRadialSymmetric{T,1,X,Y}, Diagonal{Complex{T},Vector{Complex{T}}}, Diagonal{Complex{T},Vector{Complex{T}}}}(r12, t12, r21, t21, fieldl, fieldr)
end

function getfields_lr(mirror::Mirror, field::FieldSpaceScalar{T,D,X,Y}) where {T,D,X,Y}
	fieldl = FieldSpaceScalar{T,-1,X,Y}(deepcopy(field.x_X), deepcopy(field.y_Y), deepcopy(field.e_SXY), field.λ, n1(mirror, field.λ), ref1(mirror))
	fieldr = FieldSpaceScalar{T,1,X,Y}(deepcopy(field.x_X), deepcopy(field.y_Y), deepcopy(field.e_SXY), field.λ, n2(mirror, field.λ), ref2(mirror))
	return (fieldl, fieldr)
end

function getfields_lr(mirror::Mirror, field::FieldSpaceScalarRadialSymmetric{T,D,X,Y}) where {T,D,X,Y}
	fieldl = FieldSpaceScalarRadialSymmetric{T,-1,X,Y}(deepcopy(field.r_R), deepcopy(field.e_SXY), field.λ, n1(mirror, field.λ), ref1(mirror))
	fieldr = FieldSpaceScalarRadialSymmetric{T,1,X,Y}(deepcopy(field.r_R), deepcopy(field.e_SXY), field.λ, n2(mirror, field.λ), ref2(mirror))
	return (fieldl, fieldr)
end

function getfields_lr(mirror::Mirror, field::FieldAngularSpectrumScalar{T,D,X,Y}) where {T,D,X,Y}
	fieldl = FieldAngularSpectrumScalar{T,-1,X,Y}(deepcopy(field.nsx_X), deepcopy(field.nsy_Y), deepcopy(field.e_SXY), field.λ, n1(mirror, field.λ), ref1(mirror))
	fieldr = FieldAngularSpectrumScalar{T,1,X,Y}(deepcopy(field.nsx_X), deepcopy(field.nsy_Y), deepcopy(field.e_SXY), field.λ, n2(mirror, field.λ), ref2(mirror))
	return (fieldl, fieldr)
end

function getfields_lr(mirror::Mirror, field::FieldAngularSpectrumScalarRadialSymmetric{T,D,X,Y}) where {T,D,X,Y}
	fieldl = FieldAngularSpectrumScalarRadialSymmetric{T,-1,X,Y}(deepcopy(field.nsr_R), deepcopy(field.e_SXY), field.λ, n1(mirror, field.λ), ref1(mirror))
	fieldr = FieldAngularSpectrumScalarRadialSymmetric{T,1,X,Y}(deepcopy(field.nsr_R), deepcopy(field.e_SXY), field.λ, n2(mirror, field.λ), ref2(mirror))
	return (fieldl, fieldr)
end

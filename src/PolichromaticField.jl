mutable struct FieldPolychromatic{T,D,F<:AbstractFieldMonochromatic{T,D}} <: AbstractField{T,D}
	fields::Vector{F}
end

function FieldPolychromatic_monochromaticpulse(fieldMono::F, t, It, λ_v) where {F<:AbstractFieldMonochromatic{T,D}} where {T,D}
	sizeλ = length(λ_v)
	It .= .√It
	length(It) == length(t) || ToDo()
	λ0 = fieldMono.λ

	fields = Vector{F}(undef, sizeλ)
	c = 299792458

	ω0 = 2π * c / λ0
	for iλ in eachindex(λ_v)
		Aλ = zero(Complex{T})
		ω = 2π * c / λ_v[iλ]
		Δω = (ω - ω0)
		for iT in 2:length(It)
			# Aλ[iλ] = √It[iT] * (exp(im * t[iT-1] * (λ0 - λ_v[iλ])) -
				# exp(im * t[iT] * (λ0 - λ_v[iλ]))) / im / (λ_v[iλ] - λ0)
			if abs(Δω) > 1E-15
				Aλ += It[iT] * (sin(t[iT] * Δω) -
					sin(t[iT-1] * Δω)) / Δω
			else
				Aλ += It[iT] * (t[iT] - t[iT-1])
			end
		end
		Aλ /= 2π
		fields[iλ] = deepcopy(fieldMono)
		vec(fields[iλ].e_SXY) .*= Aλ
		fields[iλ].λ = λ_v[iλ]
	end
	return FieldPolychromatic{T,D,F}(fields)
end

function lightinteraction(comp, poly::FieldPolychromatic)
	field = [lightinteraction(comp, poly.fields[i]) for i in eachindex(poly.fields)]
	lfield = first.(field)
	rfield = last.(field)
	lpoly = FieldPolychromatic(lfield)
	rpoly = FieldPolychromatic(rfield)
	return (lpoly, rpoly)
end

struct Mirror{T} <: AbstractPropagationComponent{T}
    R::JolabFunction1D{T,T}
    n₁::JolabFunction1D{T,Complex{T}}
    n₂::JolabFunction1D{T,Complex{T}}
    ref::ReferenceFrame{T}
    Mirror{T}(R, n₁, n₂, ref) where T = new{T}(R, n₁, n₂, ref)
end
Mirror(R, n₁, n₂, ref) = Mirror{Float64}(R,n₁,n₂,ref)

function PropagationCoefficientScalar(mirror::Mirror{T}, λ::Real) where T
    tλ = complex(√(1 - mirror.R(λ)));
    rλ = complex(√mirror.R(λ));
    r₁₂ = JolabFunction1D{T,Complex{T}}(rλ)
    r₂₁ = JolabFunction1D{T,Complex{T}}(-rλ)
    t₁₂ = JolabFunction1D{T,Complex{T}}(tλ)
    t₂₁ = JolabFunction1D{T,Complex{T}}(tλ)
	return PropagationCoefficientScalar{T}(rλ, tλ, -rλ, tλ, λ, mirror.n₁(λ), mirror.ref, mirror.n₂(λ), mirror.ref)
end

function PropagationCoefficientVectorial(mirror::Mirror{T}, λ::Real) where T
    tλ = complex(√(1 - mirror.R(λ)));
    rλ = complex(√mirror.R(λ));
    r₁₂ = JolabFunction1D{T,Complex{T}}(rλ)
    r₂₁ = JolabFunction1D{T,Complex{T}}(-rλ)
    t₁₂ = JolabFunction1D{T,Complex{T}}(tλ)
    t₂₁ = JolabFunction1D{T,Complex{T}}(tλ)
    return PropagationCoefficientVectorial{T}(-r₁₂, r₁₂, t₁₂, t₁₂, -r₂₁, r₂₁, t₂₁, t₂₁, λ, mirror.n₁(λ), mirror.ref, mirror.n₂(λ), mirror.ref)
end

coefficientscallar(mirror::Mirror, λ) = PropagationCoefficientScalar(mirror, λ)

function lightinteraction(mirror::Mirror, angspe::AbstractFieldAngularSpectrum)
	vectorial = size(angspe.e_SXY, 1) == 3;
	if vectorial
		coef = propagationcoefficientsvectorial(mirror, angspe.λ)
	else
		coef = propagationcoefficientscalar(mirror, angspe.λ)
	end
	(angspeforward, angspebackward) = lightinteraction(coef, angspe)
end

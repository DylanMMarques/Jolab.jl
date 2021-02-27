mutable struct Material{T <: Real}
    n::JolabFunction{T}
end

n(mat::Material, λ) = mat.n(λ)
(mat::Material)(λ) = n(mat, λ)

Base.convert(::Type{Material{T}}, a) where {T} = Material(convert(JolabFunction{T}, a))

export Medium

struct Medium{T, F} <: AbstractMedium{T}
    n::F
end

function Medium(::Type{T}, n::N) where {T, N}
    Medium{T,N}(n)
end
Medium(n) = Medium(Float64, n) 

const DefinedMedium{T} = Medium{T, <:Union{T, Complex{T}}}

(med::Medium{T, <:RealOrComplex})(wavelength::Real) where T = med.n

Base.isequal(m1::Medium, m2::Medium) = isequal(m1.n, m2.n)
Base.isapprox(m1::Medium, m2::Medium; kwarg...) = isapprox(m1.n, m2.n; kwarg...)

Base.convert(::Type{Medium{T, T}}, medium::Medium) where T = Medium{T, T}(medium.n) 
Base.convert(::Type{Medium{T, Complex{T}}}, medium::Medium) where T = Medium{T, Complex{T}}(medium.n)

is_complex_medium(med::Medium{<:Any, <:Complex}) = imag(med.n) > 0
is_complex_medium(med::Medium{<:Any, <:Real}) = false
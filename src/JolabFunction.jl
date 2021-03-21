mutable struct JolabFunction{T <: Real, A <: Union{T, Complex{T}, Function, AbstractExtrapolation{T}, AbstractInterpolation{T}}}
    y::A
end

JolabFunction(f::A) where {A<:Function} = JolabFunction{Float64,A}(f)
JolabFunction{T}(f::A) where {T,A<:Union{AbstractExtrapolation{T},AbstractInterpolation{T}}} = JolabFunction{T,A}(f)
function JolabFunction{T}(f::A) where {T,A<:Union{AbstractExtrapolation,AbstractInterpolation}}
    f2(x) =  T(f(x))
    JolabFunction{T,Function}(f2)
end

(fun::JolabFunction{T,A})(x) where {T,A<:Union{AbstractExtrapolation,AbstractInterpolation}} = fun.y(x)
(fun::JolabFunction{T,A})(x) where {T,A<:Function} = Complex{T}(fun.y(x))
(fun::JolabFunction{T,A})(x) where {T,A<:Number} = fun.y

Base.convert(::Type{JolabFunction{T}}, f::A) where {T<:Real, A<:Real} = JolabFunction{T,T}(T(f))
Base.convert(::Type{JolabFunction{T}}, f::A) where {T<:Real, A<:Number} = JolabFunction{T,Complex{T}}(Complex{T}(f))
Base.convert(::Type{JolabFunction{T}}, f::A) where {T<:Real, A<:Function} = JolabFunction{T,A}(f)
Base.convert(::Type{JolabFunction{T}}, f::A) where {T<:Real, A<:Union{AbstractExtrapolation, AbstractInterpolation}} = JolabFunction{T,A}(f)

Base.convert(::Type{JolabFunction{A}}, fun::JolabFunction{B,C}) where {A<:Real, B<:Real, C<:Real} = JolabFunction{A,A}(A(fun.y))
Base.convert(::Type{JolabFunction{A}}, fun::JolabFunction{B,C}) where {A<:Real, B<:Real, C<:Number} = JolabFunction{A,Complex{A}}(Complex{A}(fun.y))
Base.convert(::Type{JolabFunction{A}}, fun::JolabFunction{B,C}) where {A<:Real, B<:Real, C<:Function} = JolabFunction{A,C}(fun.y)
function Base.convert(::Type{JolabFunction{A}}, fun::JolabFunction{B,C}) where {A<:Number, B<:Number, C<:Union{AbstractExtrapolation, AbstractInterpolation}}
    f(x) = A(fun.y(x))
    return JolabFunction{A,Function}(f)
end

mutable struct JolabFunction2D{T <: Real, A <: Union{T, Complex{T}, Function, AbstractExtrapolation{T}, AbstractInterpolation{T}}}
    y::A
end

(fun::JolabFunction2D{T,A})(x, y) where {T,A<:Function} = T(fun.y(x,y))
(fun::JolabFunction2D{T,A})(x, y) where {T,A<:Union{AbstractExtrapolation,AbstractInterpolation}} = fun.y(x,y)
(fun::JolabFunction2D{T,A})(x, y) where {T,A<:Number} = fun.y

Base.convert(::Type{JolabFunction2D{T}}, f::A) where {T<:Real, A<:Real} = JolabFunction2D{T,T}(T(f))
Base.convert(::Type{JolabFunction2D{T}}, f::A) where {T<:Real, A<:Number} = JolabFunction2D{T,Complex{T}}(Complex{T}(f))
Base.convert(::Type{JolabFunction2D{T}}, f::A) where {T<:Real, A<:Function} = JolabFunction2D{T,A}(f)
Base.convert(::Type{JolabFunction2D{T}}, f::A) where {T<:Real, A<:Union{AbstractExtrapolation, AbstractInterpolation}} = JolabFunction2D{T,A}(f)

Base.convert(::Type{JolabFunction2D{A}}, fun::JolabFunction2D{B,C}) where {A<:Real, B<:Real, C<:Real} = JolabFunction2D{A,A}(A(fun.y))
Base.convert(::Type{JolabFunction2D{A}}, fun::JolabFunction2D{B,C}) where {A<:Real, B<:Real, C<:Number} = JolabFunction2D{A,Complex{A}}(Complex{A}(fun.y))
Base.convert(::Type{JolabFunction2D{A}}, fun::JolabFunction2D{B,C}) where {A<:Real, B<:Real, C<:Function} = JolabFunction2D{A,C}(fun.y)
function Base.convert(::Type{JolabFunction{A}}, fun::JolabFunction2D{B,C}) where {A<:Number, B<:Number, C<:Union{AbstractExtrapolation, AbstractInterpolation}}
    f(x,y) = A(fun.y(x,y))
    return JolabFunction2D{A,Function}(f)
end

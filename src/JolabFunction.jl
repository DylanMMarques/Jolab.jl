abstract type JolabFunction{X<:Number,Y<:Number} end

struct JolabFunction1D{X,Y} <: JolabFunction{X,Y}
    f::FunctionWrapper{Y, NTuple{1, X}}
end

struct JolabFunction2D{X,Y} <: JolabFunction{X,Y}
    f::FunctionWrapper{Y, NTuple{2, X}}
end

@inline (fun::JolabFunction1D{X,Y})(x) where {X,Y} = fun.f(x)
@inline (fun::JolabFunction2D{X,Y})(x, y) where {X,Y} = fun.f(x, y)

typeofx(f::JolabFunction1D{X,Y}) where {X,Y} = X
typeofx(f::JolabFunction2D{X,Y}) where {X,Y} = X
typeofy(f::JolabFunction1D{X,Y}) where {X,Y} = Y
typeofy(f::JolabFunction2D{X,Y}) where {X,Y} = Y

Base.convert(::Type{JolabFunction1D{X,Y}}, fnc::Function) where {X,Y} = JolabFunction1D{X,Y}(fnc)
Base.convert(::Type{JolabFunction2D{X,Y}}, fnc::Function) where {X,Y} = JolabFunction2D{X,Y}(fnc)

function JolabFunction1D{X,Y}(val::Number) where {X,Y}
    f(x) = val
    JolabFunction1D{X,Y}(FunctionWrapper{Y,NTuple{1,X}}(f))
end
Base.convert(::Type{JolabFunction1D{X,Y}}, val::Number) where {X,Y} = JolabFunction1D{X,Y}(val)

function JolabFunction2D{X,Y}(val::Number) where {X,Y}
    f(x,y) = val
    JolabFunction1D{X,Y}(FunctionWrapper{Y,NTuple{2,X}}(f))
end
Base.convert(::Type{JolabFunction2D{X,Y}}, val::Number) where {X,Y} = JolabFunction2D{X,Y}(val)

function JolabFunction2D(x::AbstractRange{X}, y::AbstractRange{Y}, z::AbstractArray{A,2}) where {X<:Real, Y<:Real, A<:Number}
    itp = interpolate(z, BSpline(Linear()))
    itp = scale(itp, x, y)
    itp = extrapolate(itp, zero(A))
    return JolabFunction2D{promote_type(X,Y),A}(itp)
end

function JolabFunction1D(x::AbstractRange{X}, z::AbstractVector{A}) where {X<:Real, A<:Number}
    itp = interpolate(z, BSpline(Linear()))
    itp = scale(itp, x)
    itp = extrapolate(itp, zero(A))
    return JolabFunction1D{X,A}(itp)
end

function JolabFunction2D(x::AbstractVector{X}, y::AbstractVector{Y}, z::AbstractArray{A,2}) where {X<:Real, Y<:Real, A<:Number}
    itp = interpolate((x,y), z, BSpline(Linear()))
    itp = extrapolate(itp, zero(A))
    return JolabFunction2D{promote_type(X,Y),A}(itp)
end

function JolabFunction1D(x::AbstractVector{T}, z::AbstractVector{A}) where {T<:Real, A<:Number}
    itp = interpolate((x,), z, BSpline(Linear()))
    itp = extrapolate(itp, zero(A))
    return JolabFunction1D{T,A}(itp)
end

function JolabFunction1D(itp::Interpolations.GriddedInterpolation{T,1,TCoefs,IT,K}) where {T,TCoefs,IT,K}
    f(x) = itp(x)
    return JolabFunction1D{eltype(itp.knots[1]),T}(f)
end

function JolabFunction2D(itp::Interpolations.GriddedInterpolation{T,2,TCoefs,IT,K}) where {T,TCoefs,IT,K}
    f(x,y) = itp(x,y)
    return JolabFunction2D{eltype(itp.knots[1]),T}(f)
end

function JolabFunction1D(etp::AbstractExtrapolation{T,1,ITPT,IT}) where {T,ITPT<:Interpolations.GriddedInterpolation,IT}
    f(x) = etp(x)
    return JolabFunction1D{eltype(etp.itp.knots[1]),T}(f)
end

function JolabFunction2D(etp::AbstractExtrapolation{T,2,ITPT,IT}) where {T,ITPT<:Interpolations.GriddedInterpolation,IT}
    f(x,y) = etp(x,y)
    return JolabFunction2D{eltype(etp.itp.knots[1]),T}(f)
end

function JolabFunction1D(itp::Interpolations.ScaledInterpolation{T,1,ITPT,IT,RT}) where {T,ITPT,IT,RT}
    f(x) = itp(x)
    return JolabFunction1D{eltype(itp.ranges[1]),T}(f)
end

function JolabFunction2D(itp::Interpolations.ScaledInterpolation{T,2,ITPT,IT,RT}) where {T,ITPT,IT,RT}
    f(x,y) = itp(x,y)
    return JolabFunction2D{eltype(itp.ranges[1]),T}(f)
end

function JolabFunction1D(etp::AbstractExtrapolation{T,1,ITPT,IT}) where {T,ITPT<:Interpolations.ScaledInterpolation,IT}
    f(x) = etp(x)
    return JolabFunction1D{eltype(etp.itp.ranges[1]),T}(f)
end

function JolabFunction2D(etp::AbstractExtrapolation{T,2,ITPT,IT}) where {T,ITPT<:Interpolations.ScaledInterpolation,IT}
    f(x,y) = etp(x,y)
    return JolabFunction2D{eltype(etp.itp.ranges[1]),T}(f)
end

function Base.convert(::Type{JolabFunction1D{X,Y}}, fnc::JolabFunction1D{A,B}) where {X,Y,A,B}
    f(x) = convert(Y, fnc.f(convert(X, x)))
    return JolabFunction1D{X,Y}(FunctionWrapper{Y,NTuple{1,X}}(f))
end

function Base.convert(::Type{JolabFunction2D{X,Y}}, fnc::JolabFunction2D{A,B}) where {X,Y,A,B}
    f(x,y) = convert(Y, fnc.f(convert(X, x), convert(Y,y)))
    return JolabFunction2D{X,Y}(FunctionWrapper{Y,NTuple{2,X}}(f))
end

function extrapolation(fun::JolabFunction1D{X,Y}, x::AbstractRange) where {X,Y}
    val = [fun.f(xi) for xi in x]
    itp = interpolate(val, BSpline(Linear()))
    itp = scale(itp, x)
    itp = extrapolate(itp, zero(Y))
    f(x1) = itp(x1)
    return JolabFunction1D{X,Y}(f)
end

function extrapolation(fun::JolabFunction2D{X,Y}, x::AbstractRange, y::AbstractRange) where {X,Y}
    val = [fun.f(xi,yi) for xi in x, yi in y]
    itp = interpolate(val, BSpline(Linear()))
    itp = scale(itp, x, y)
    itp = extrapolate(itp, zero(Y))
    f(x1,y1) = itp(x1,y1)
    return JolabFunction2D{X,Y}(f)
end

function Base.:*(fun::JolabFunction1D{X,Y}, a::Real) where {X,Y}
    b(x) = fun.f(x) * a
    return JolabFunction1D{X,Y}(b)
end

using Enzyme
using Test, Enzyme
import .EnzymeRules: forward

function EnzymeRules.forward(func::Const{typeof(√)}, ::Type{<:Duplicated}, x::Duplicated{<:Complex{<:Real}})
    # println("Using custom rule_were2!")
    ret = func.val(x.val)
    dval = 1 / 2ret * x.dval
    return Duplicated(ret, dval)
end

function EnzymeRules.forward(func::Const{typeof(√)}, ::Type{<:Duplicated}, x::Const{<:Complex{<:Real}})
    # println("Using custom rule_were2!")
    ret = func.val(x.val)
    return Duplicated(ret, zero(ret))
end

f(x) = 1/2√(complex(x))

test_sqrt(x) = autodiff(Forward, √, Duplicated, Duplicated(x, one(x))) == (√(complex(x)), f(x))
@test test_sqrt(0.5)
@test test_sqrt(-0.5 - im)
@test test_sqrt(-0.5 + 0im)
@test test_sqrt(-0.5 + im)
@test test_sqrt(-0.5 - im)


f(x) = exp(im + (x + im))
autodiff(Forward, f, Duplicated, Duplicated(1.0, 1.0))
autodiff(Forward, f, Duplicated, Duplicated(1.0, 1.0))


f(nsr, n) = √(complex(1 - (nsr / n))^2)
autodiff(Forward, f, Duplicated, Duplicated(0.1, 1.0), Duplicated(1.0, 1.0))
autodiff(Forward, f, Duplicated, Duplicated(0.1, 1.0), Duplicated(1.0 +0im, 1.0+0im))

autodiff(Forward, i-> sqrt((i)), Duplicated, Duplicated(0.1, 1.0))

using BenchmarkTools
@btime √(complex(1))
@btime complex(1)^(1/2)
@btime √(1.0)
f1(x) = x .^ (1/2)
f2(x) = .√(x)
a = rand(ComplexF64, 1000);
@btime f1(a);
@btime f2(a);
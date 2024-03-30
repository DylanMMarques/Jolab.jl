module JolabEnzymeOptimExt

using Jolab, Enzyme, Optim, StaticArrays
import Jolab: modecondition

# Find the β that satisfies the mode condition, i.e. the mode condition == 0
# The algorithm performs better if the initial guess is close to the solution
function Jolab.wavefunction_solutions(profile, λ, β, m_i)
    f(β) = modecondition(profile, λ, m_i, β[1])^2
    df_dx(β) = autodiff_deferred(Enzyme.Forward, f, Duplicated, Duplicated(β, one(β)))
    df2_dx2(β) = autodiff(Enzyme.Forward, df_dx, Duplicated, Duplicated(β, one(β)))
    function fgh!(F, G, H, β)
        ((val, der), (_der, der2)) = df2_dx2(β[1])
        if G !== nothing
            G[1] = der
        end
        if H !== nothing 
            H[1] = der2
        end
        if F !== nothing
            return val
        end
    end
    inner_optimizer = NewtonTrustRegion()
    res = optimize(Optim.only_fgh!(fgh!), (@MArray [β]), inner_optimizer)
    Optim.minimizer(res)[1]
end
end
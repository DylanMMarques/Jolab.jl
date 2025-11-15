module JolabEnzymeOptimExt

using Jolab, Enzyme, Optim, StaticArrays
import Jolab: modecondition

# Find the β that satisfies the mode condition, i.e. the mode condition == 0
# The algorithm performs better if the initial guess is close to the solution
function Jolab.wavefunction_solutions(profile, λ, β, m_i)
    f(β) = modecondition(profile, λ, m_i, β[1])^2
    df_dx(β) = autodiff(Enzyme.ForwardWithPrimal, f, Duplicated, Duplicated(β, one(β)))
    # df2_dx2(β) = autodiff(Enzyme.Forward, df_dx, Duplicated, Duplicated(β, one(β)))
    function fg!(F, G, β)
        (der, val) = df_dx(β[1])
        if G !== nothing
            G[1] = der
        end
        if F !== nothing
            return val
        end
    end
    inner_optimizer = ConjugateGradient()
    res = optimize(Optim.only_fg!(fg!), (@MArray [β]), inner_optimizer)
    Optim.minimizer(res)[1]
end
end


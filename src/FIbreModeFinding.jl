using .Optim, .Enzyme
# Find the β that satisfies the mode condition, i.e. the mode condition == 0
# The algorithm performs better if the initial guess is close to the solution
function wavefunction_solutions(profile, β, m_i)
    f(β) = modecondition(profile, 1500E-9, m_i, β[1])^2
    g(β) = autodiff(Enzyme.Forward, f, Duplicated, Duplicated(β[1], one(β[1])))
    function fgh!(F, G, H, β)
        (val, der) = g(β[1])
        if G !== nothing
            G[1] = der
        end
        if H !== nothing # Calculate numerical derivative
            tol = β[1] * 1E-12
            H[1] = (g(β[1] + tol)[2] - der) / tol
        end
        if F !== nothing
            return val
        end
    end
    inner_optimizer = NewtonTrustRegion()
    res = optimize(Optim.only_fgh!(fgh!), (@MArray [β]), inner_optimizer)
    Optim.minimizer(res)[1]
end
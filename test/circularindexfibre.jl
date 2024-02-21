using Jolab, Test

profile = CircularStepIndexProfile(100E-6, 0.2, Medium(1.55))

fibre = Fibre(profile, 1, Medium.((1,1)), (ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0,0,1), (0,0,0))))

modes = Jolab.findmodes!(profile, 1500E-9)

(βmin, βmax) = (6.438349049607e6, 6.492624817418906e6)

using NonlinearSolve
f(u, p) = u * u - 2.0

uspan = (βmin, βmax) # brackets
function g(uspan, m)

    sols = zeros(Float64, 10000)
    i = 1
    function find_all_zeros(int)
        f(β, p) = Jolab.modecondition(profile, 1500E-9, m, β)
        prob_int = IntervalNonlinearProblem(f, int)
        sol = solve(prob_int, rtol = 1E-12)
        @show f(sol.u, 0), uspan
        if abs(f(sol.u, 0)) < 2
            sols[i] = sol.u
            i += 1 
            find_all_zeros((int[1], sol.u - 10))
            find_all_zeros((sol.u + 10, int[2]))
        end
    end
    find_all_zeros(uspan)
    (sols, i)

end

function g(uo, u2, m)
    f(β, m) = Jolab.modecondition(profile, 1500E-9, m, β)
    prob = IntervalNonlinearProblem(f, (uo, u2), m)
    solve(prob).u
end

function find_modex(m_i)
    m_i = 1
    (βmin, βmax) = (6.438349049607e6, 6.492624817418906e6)
    size_M = 10000
    function rough_zeros(size_M)
        x = range(βmin+5, βmax-5, length = size_M)
        initial_guess = Jolab.modecondition.(Ref(profile), 1500E-9, m_i, x)
        int_line = initial_guess[1:end-1] .* initial_guess[2:end] .< 0 .&& initial_guess[1:end-1].^2 .< 100
        ind = findall(int_line)
        (x[ind] .+ x[ind .+ 1]) ./ 2
    end
    xi = rough_zeros(size_M)
    old_length = 0
    while length(xi) != old_length
        old_length = length(xi)
        size_M = size_M * 2
        xi = rough_zeros(size_M)
    end
    g.(xi, xi .+ step(x), m_i)
end

using FiniteDiff
function g2(uo, m_i)
    f(β) = Jolab.modecondition(profile, 1500E-9, m_i, β[1])^2
    g(β) = autodiff(Enzyme.Forward, f, Duplicated, Duplicated(β, one(β)))
    function fgh!(F, G, H, β)
        (val, der) = g(β[1])
        if G !== nothing
            G[1] = der
        end
        if H !== nothing # Calculate numerical derivative
            tol = β * 1E-12
            H[1] = (g(β + tol) - der) / tol
        end
        if F !== nothing
            return val
        end
    end
    inner_optimizer = NewtonTrustRegion()
    res = optimize(Optim.only_fgh!(fgh!), (@MArray [uo]), inner_optimizer)
    Optim.minimizer(res)[1]
end

function find_local_minima(f, x, initial_size)
    local_minima = Vector{Float64}(undef, initial_size)
    ind_zeros = 1
    value = f(x[1])
    next_value = f(x[2])
    @inbounds for ix in view(x, 3:length(x))
        prev_value = value
        value = next_value
        next_value = f(ix)
        if value < prev_value && value < next_value
            local_minima[ind_zeros] = ix - step(x)
            ind_zeros += 1
        end
    end
    local_minima[1:(ind_zeros-1)]
end

function find_modex3(m_i)
    m_i = 1
    f(x) = Jolab.modecondition(profile, 1500E-9, m_i, x)^2
    (βmin, βmax) = (6.438349049607e6, 6.492624817418906e6)
    initial_guess = 100
    initial_diff = 500
    i = 1
    x = range(βmin + 1, βmax - 1, length = initial_diff * i)
    xi = find_local_minima(f, x, initial_guess)
    cur_len = length(xi)
    old_length = 0
    while cur_len != old_length
        old_length = cur_len
        x = range(βmin + 1, βmax - 1, length = initial_diff * i)
        xi = find_local_minima(f, x, initial_guess)
        cur_len = length(xi)
        i += 1
    end
    return xi
    g2.(xi, m_i)
end

function find_modex2(m)
    T = Float64
    condition(β) = Jolab.modecondition(profile, 1500E-9, m, β)

    (βmin, βmax) = (6.438349049607e6, 6.492624817418906e6)

    roots = Jolab.find_zeros(condition, βmin, βmax, k = 100)
    filter(i -> condition(i) < 1, roots)
end
@time a = find_modex(1)
@btime b = find_modex2(1);
@time c = find_modex3(1);


## Tests

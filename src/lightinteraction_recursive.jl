function _lightinteraction_recursivegridded!(fields_l, fields_r, coefs, fieldi::AbstractField{T,D}; rtol = 1E-3::Real, printBool = true, maximum_iterations = 10000) where {T,D}
    sizeL = length(coefs) + 1;
    length(fields_l) == length(fields_r) == sizeL || error()
    
    # int_l, int_r = zeros(T, sizeL), zeros(T, sizeL)
    
    fill_zeros!.(fields_r)
    fill_zeros!.(fields_l)
    fields_aux_r = deepcopy.(fields_r)
    fields_aux_l = deepcopy.(fields_l)
    fields2_l = (fields_l[1], deepcopy.(fields_l[2:sizeL])...)
    fields2_r = (deepcopy.(fields_r[1:sizeL-1])..., fields_r[sizeL])

    
    rtol = intensity(fieldi) * rtol^2
    
    if D == Forward
        copy!(fields_r[1].e, fieldi.e)
    else
        copy!(fields_l[sizeL].e, fieldi.e)
    end
    initial_int = intensity(fieldi)
    int_r = MVector(intensity.(fields_r))
    int_l = MVector(intensity.(fields_l))
    
    i = 1
    min_int = initial_int
    converge = false
    
    while true
        if isodd(i)
    		(toSave_l, toSave_r) = (fields2_l, fields2_r)
    		(iE_l, iE_r) = (fields_l, fields_r)
    	else
    		(toSave_l, toSave_r) = (fields_l, fields_r)
    		(iE_l, iE_r) = (fields2_l, fields2_r)
    	end
    
        for mls in 1:sizeL-1
    	    if int_r[mls] > 1E-15
    	    	# iE_r[mls].frame == coefs[mls].fieldl.frame || tobedone()
                _light_interaction!(fields_aux_l[mls], fields_aux_r[mls+1], coefs[mls], iE_r[mls])
    	    	_unchecked_add!(toSave_l[mls], fields_aux_l[mls])
    	    	_unchecked_add!(toSave_r[mls+1], fields_aux_r[mls+1])
    	    end
    	    fill_zeros!(iE_r[mls])
        end
        for mls in 2:sizeL
    	    if int_l[mls] > 1E-15
    	    	# iE_l[mls].frame == coefs[mls-1].fieldr.frame || tobedone()
    	    	_light_interaction!(fields_aux_l[mls-1], fields_aux_r[mls], coefs[mls-1], iE_l[mls])
    	    	_unchecked_add!(toSave_l[mls-1], fields_aux_l[mls-1])
    	    	_unchecked_add!(toSave_r[mls], fields_aux_r[mls])
    	    end
    	    fill_zeros!(iE_l[mls])
        end
        int_l .= intensity.(toSave_l)
        int_r .= intensity.(toSave_r)
    	now_int = sum(view(int_l,2:sizeL)) + sum(view(int_r, 1:sizeL-1))
    
    	if now_int < rtol 
            if printBool
                println(""); 
                println("Interactions until convergence: ", i)
            end
            converge = true
            break
        end
    
    	(now_int < min_int) && (min_int = now_int)
    
    	if sum(int_l) + sum(int_r) > 10 * initial_int || now_int > 10 * min_int
    		println("")
    		println("Light intensity propagating forward:", int_r)
    		println("Light intensity propagating backward:", int_l)
    		println("lightinteraction_recursivegridded is not converging. Current number of iterations:", i)
    		converge = false
    		# break
    	end
    
    	i > maximum_iterations && (println("Max number of iterations achieved. Current light intensity:", (sum(int_l) + sum(int_r)) / initial_int); converge = false; break)
    	i += 1
    	if printBool
                println("Iteration number: ", i)
    		println("Light intensity propagating forward:", int_r)
    		println("Light intensity propagating backward:", int_l)
    		println("convergence condition: ", sum(view(int_l,2:sizeL)) + sum(view(int_r, 1:sizeL-1)), " < ", rtol)
    	end
    end

    return fields_l[1], fields_r[end]
end

function lightinteraction_recursivegridded(comp, fieldi::AbstractField{T,D}; kwargs...) where {T,D}
    iter = D == Forward ? identity : Iterators.reverse
    field_transmitted(coef, field_i) = last(reverse_if_backward(D, forward_backward_field(coef, field_i)))

    out = accumulate((field_i, comp) -> begin 
            sol = solver(comp, first(field_i))
            field_t = field_transmitted(sol, first(field_i))
            (field_t, sol)
        end
        , iter(comp), init = (fieldi, nothing))
    solvers = last.(out)
    _fields_transmitted = first.(out)
    fields_transmitted = (deepcopy(fieldi), _fields_transmitted...)
    fields_reflected = reverse_direction.(deepcopy.(fields_transmitted))

    # _fields_transmitted = accumulate((field_i, coef) -> field_transmitted(forward_backward_field(coef, field_i)), solvers, init = fieldi)
    (fields_backward, fields_forward) = reverse_if_backward(D, (fields_reflected, fields_transmitted))
    _lightinteraction_recursivegridded!(fields_backward, fields_forward, solvers, fieldi; kwargs...)
end

function _lightinteraction_recursivegridded!(fields_l, fields_r, coefs, fieldi; rtol = 1E-3one(T)::Real, printBool = true)
    T = Float64
    sizeL = length(coefs) + 1;
    length(fields_l) == length(fields_r) == sizeL || error()
    
    int_l, int_r = zeros(T, sizeL), zeros(T, sizeL)
    
    fields_aux_r = similar.(fields_r)
    fields_aux_l = similar.(fields_l)
    fields2_l = similar.(fields_l)
    fields2_r = similar.(fields_r)
    fields_r[sizeL] = fields2_r[sizeL]
    fields_l[1] = fields2_l[1]
    
    rtol = intensity(fieldi) * rtol^2
    
    dir_fieldi = Forward
    if dir_fieldi == Forward
        fields_r[1].e .= fieldi.e
    else
        fields_l[sizeL].e .= fieldi.e
    end
    initial_int = intensity(fieldi)
    dir_fieldi == Forward ? int_r[1] = initial_int : int_l[sizeL] = initial_int
    
    i = 1
    toSave_l, toSave_r = fields_l, fields_r
    min_int = initial_int
    converge = false
    
    while true
    	if i % 2 == 1
    		(toSave_l, toSave_r) = (fields_l, fields_r)
    		(iE_l, iE_r) = (fields2_l, fields2_r)
    	else
    		(toSave_l, toSave_r) = (fields2_l, fields2_r)
    		(iE_l, iE_r) = (fields_l, fields_r)
    	end
    
        for mls in 1:sizeL-1
    	    if int_r[mls] > 1E-15
    	    	iE_r[mls].ref == coefs[mls].fieldl.ref || tobedone()
    	    	lightinteraction!(fields_aux_l[mls], fields_aux_r[mls+1], coefs[mls], iE_r[mls])
    	    	_unchecked_add!(toSave_l[mls], fields_aux_l[mls])
    	    	_unchecked_add!(toSave_r[mls+1], fields_aux_r[mls+1])
    	    end
    	    fill_zeros!(iE_r[mls])
        end
        for mls in 2:sizeL
    	    if int_l[mls] > 1E-15
    	    	iE_l[mls].ref == coefs[mls-1].fieldr.ref || tobedone()
    	    	lightinteraction!(fields_aux_l[mls-1], fields_aux_r[mls], coefs[mls-1], iE_l[mls])
    	    	_unchecked_add!(toSave_l[mls-1], fields_aux_l[mls-1])
    	    	_unchecked_add!(toSave_r[mls], fields_aux_r[mls])
    	    end
    	    fill_zeros!(iE_l[mls])
        end
    	int_l .= intensity.(toSave_l)
    	int_r .= intensity.(toSave_r)
    	now_int = sum(view(int_l,2:sizeL)) + sum(view(int_r, 1:sizeL-1))
    
    	now_int < rtol && (println(""); println("Interactions until convergence: ", i); converge = true; break)
    
    	(now_int < min_int) && (min_int = now_int)
    
    	if sum(int_l) + sum(int_r) > 10 * initial_int || now_int > 10 * min_int
    		println("")
    		println("Light intensity propagating forward:", int_r)
    		println("Light intensity propagating backward:", int_l)
    		println("lightinteraction_recursivegridded is not converging. Current number of iterations:", i)
    		converge = false
    		# break
    	end
    
    	i > 10000 && (println("Max number of iterations achieved. Current light intensity:", (sum(int_l) + sum(int_r)) / initial_int); converge = false; break)
    	i += 1
    	if (i % 100 == 99) && printBool
    		println("")
    		println("Light intensity propagating forward:", int_r)
    		println("Light intensity propagating backward:", int_l)
    		println("convergence condition: ", sum(view(int_l,2:sizeL)) + sum(view(int_r, 1:sizeL-1)), " < ", rtol)
    	end
    end
    return fields_l, fields_r, converge
end

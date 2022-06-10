struct Interpolation{T, X <: AbstractVector{T}} <: AbstractOpticalComponent{T}
    x::X
    y::X
end

function coefficient_specific(itp::Interpolation, field::AbstractFieldMonochromatic)
    (fieldl, fieldr) = getfields_lr(itp, field)
    coefficient_specific!(fieldl, fieldr, itp, field)
    return (fieldl, fieldr)
end

#function coefficient_specific!(fieldl::FieldSpaceScalar, fieldr::FieldSpaceScalar, itp::Interpolation, space::FieldSpaceScalar{T}) where T
#    checkapplicability(fieldl, fieldr, itp, space)
#    field = dir(space) > 0 ? fieldr : fieldl
#
#    e_itp = LinearInterpolation((space.x_X, space.y_Y), reshape(space.e_SXY, length(space.x_X), length(space.y_Y)), extrapolation_bc = zero(Complex{T}))
#
#	cart = CartesianIndices(field)
#	@inbounds @simd for i in iterator_index(field)
#		field.e_SXY[i] = e_itp(field.x_X[cart[i][2]], field.y_Y[cart[i][3]])
#	end
#    return (fieldl, fieldr)
#end

for (field_type, x, y) in ((:FieldAngularSpectrumScalar, :nsx_X, :nsy_Y), (:FieldSpaceScalar, :x_X, :y_Y))
    eval(quote
        function coefficient_specific!(fieldl::$field_type, fieldr::$field_type, itp::Interpolation, fieldi::$field_type{T}) where T
            checkapplicability(fieldl, fieldr, itp, fieldi)
            field = dir(fieldi) > 0 ? fieldr : fieldl

            e_itp = LinearInterpolation((fieldi.$x, fieldi.$y), reshape(fieldi.e_SXY, length(fieldi.$x), length(fieldi.$y)), extrapolation_bc = zero(Complex{T}))

        	cart = CartesianIndices(field)
        	@inbounds @simd for i in iterator_index(field)
        		field.e_SXY[i] = e_itp(field.$x[cart[i][2]], field.$y[cart[i][3]])
        	end
            return (fieldl, fieldr)
        end
        function getfields_lr(itp::Interpolation{T1, X1}, field::$field_type{T,D,X,Y}) where {T,D,X,Y, T1,X1}
            e_SXY = zeros(Complex{T}, length(itp.x) * length(itp.y))
            if dir(field) > 0 
        	    fieldl = $field_type{T,-1,X,Y}(deepcopy(field.$x), deepcopy(field.$y), 0deepcopy(field.e_SXY), field.λ, field.n, field.ref)
        	    fieldr = $field_type{T,1,X1,Vector{Complex{T}}}(deepcopy(itp.x), deepcopy(itp.y), e_SXY, field.λ, field.n, field.ref)
            else
        	    fieldl = $field_type{T,-1,X1,Vector{Complex{T}}}(deepcopy(itp.x), deepcopy(itp.y), e_SXY, field.λ, field.n, field.ref)
        	    fieldr = $field_type{T,1,X,Y}(deepcopy(field.$x), deepcopy(field.$y), 0deepcopy(field.e_SXY), field.λ, field.n, field.ref)
            end
            return (fieldl, fieldr)
        end
    end)
end


for (field_type, r) in ((:FieldAngularSpectrumScalarRadialSymmetric, :nsr_R), (:FieldSpaceScalarRadialSymmetric, :r_R))
    eval(quote
        function getfields_lr(itp::Interpolation{T1, X1}, field::$field_type{T,D,X,Y}) where {T,D,X,Y, T1,X1}
            e_SXY = zeros(Complex{T}, length(itp.x))
            if dir(field) > 0 
        	    fieldl = $field_type{T,-1,X,Y}(deepcopy(field.$r), 0deepcopy(field.e_SXY), field.λ, field.n, field.ref)
        	    fieldr = $field_type{T,1,X1,Vector{Complex{T}}}(deepcopy(itp.x), e_SXY, field.λ, field.n, field.ref)
            else                                
        	    fieldl = $field_type{T,-1,X1,Vector{Complex{T}}}(deepcopy(itp.x), e_SXY, field.λ, field.n, field.ref)
        	    fieldr = $field_type{T,1,X,Y}(deepcopy(field.$r), 0deepcopy(field.e_SXY), field.λ, field.n, field.ref)
            end
            return (fieldl, fieldr)
        end
        function coefficient_general(itp::Interpolation{TI, XI}, field::$field_type{T,D,X,B}) where {TI,XI,T,D,X,B}
            sizeIn, sizeOut = length(field.$r), length(itp.x)
            ind_in = ones(Int64, 2sizeOut)
            ind_out = ones(Int64, 2sizeOut)
            val = zeros(Complex{T}, 2sizeOut)
        
            e_itp = LinearInterpolation((field.$r,), field.e_SXY)
            for i in eachindex(itp.x)
                weight_itp = Interpolations.weightedindexes((Interpolations.value_weights,), Interpolations.itpinfo(e_itp)..., (itp.x[i],))
                @show Interpolations.indexes(weight_itp[1])
                continue
                if weight_itp.istart == 0
                    ind_in[2i] = weight_itp.istart + 1
                    ind_out[2i] = i
                    val[2i] = weight_itp.weights[2]
                elseif weight_itp.istart == sizeIn
                    ind_in[2i-1] = weight_itp.istart
                    ind_out[2i-1] = i
                    val[2i-1] = weight_itp.weights[1]
                elseif 1 <= weight_itp.istart <= sizeIn - 1
                    (ind_in[2i-1], ind_in[2i]) = (weight_itp.istart, weight_itp.istart + 1)
                    (ind_out[2i-1], ind_out[2i]) = (i, i)
                    (val[2i-1], val[2i]) = weight_itp.weights
                end
            end
            @show val
            t12 = sparse(ind_out, ind_in, val, sizeOut, sizeIn)
            @show t12
        
            ind_in = ones(Int64, 2sizeIn)
            ind_out = ones(Int64, 2sizeIn)
            val = zeros(Complex{T}, 2sizeIn)
            e_itp = LinearInterpolation((itp.x,), Vector{Float16}(undef, sizeOut))
            for i in eachindex(field.$r)
                weight_itp = Interpolations.weightedindexes((Interpolations.value_weights,), Interpolations.itpinfo(e_itp)..., (field.$r[i],))[1]
                if weight_itp.istart == 0
                    ind_out[2i] = weight_itp.istart + 1
                    ind_in[2i] = i
                    val[2i] = weight_itp.weights[2]
                elseif weight_itp.istart == sizeOut
                    ind_out[2i-1] = weight_itp.istart
                    ind_int[2i-1] = i
                    val[2i-1] = weight_itp.weights[1]
                elseif 1 <= weight_itp.istart <= sizeOut - 1
                    (ind_out[2i-1], ind_out[2i]) = (weight_itp.istart, weight_itp.istart + 1)
                    (ind_in[2i-1], ind_in[2i]) = (i, i)
                    (val[2i-1], val[2i]) = weight_itp.weights
                end
            end
            t21 = sparse(ind_in, ind_out, val, sizeIn, sizeOut)
            (fieldl, fieldr) = getfields_lr(itp, field)
            if dir(field) > 0
               return ScatteringMatrix{T, $field_type{T,-1,X,B}, $field_type{T,1,XI,Vector{Complex{T}}}, Nothing, SparseMatrixCSC{Complex{T},Int64}}(nothing, t12, nothing, t21, fieldl, fieldr)
            else
               return ScatteringMatrix{T, $field_type{T,-1,XI,Vector{Complex{T}}}, $field_type{T,1,X,B}, Nothing, SparseMatrixCSC{Complex{T},Int64}}(nothing, t12, nothing, t21, fieldl, fieldr)
            end
        end
        function coefficient_specific!(fieldl::$field_type, fieldr::$field_type, itp::Interpolation, fieldi::$field_type{T}) where T
            checkapplicability(fieldl, fieldr, itp, fieldi)
            field = dir(fieldi) > 0 ? fieldr : fieldl

            e_itp = LinearInterpolation((fieldi.$r,), fieldi.e_SXY, extrapolation_bc = zero(Complex{T}))

        	cart = CartesianIndices(field)
        	@inbounds @simd for i in iterator_index(field)
        		field.e_SXY[i] = e_itp(field.$r[cart[i][2]])
        	end
            return (fieldl, fieldr)
        end
    end)
end

function checkapplicability(fieldl, fieldr, itp::Interpolation, fieldi::AbstractFieldMonochromatic)
    field = dir(fieldi) > 0 ? fieldr : fieldl
    fieldi.ref == field.ref || error("The interpolated field must have the same reference frame as the original field")
    isapprox(fieldi.n, field.n, atol = @tol) || error("The interpolated field must be on the same refractive index that the original field")
    isapprox(fieldi.λ, field.λ, atol = @tol) || error("The interpolated field must have the same wavelength that the original field")
end


function Interpolations.interpolate(field::Union{FieldAngularSpectrumScalar{T}, FieldSpaceScalar{T}, FieldSpaceScalarRadialSymmetric{T}, FieldAngularSpectrumScalarRadialSymmetric{T}}, x_X::X1, y_Y::X2) where {T,X1,X2}
    itp = Interpolation{T, promote_type(X1,X2)}(x_X, y_Y)
    (fieldl, fieldr) = coefficient_specific(itp, field)
    return dir(field) > 0 ? fieldr : fieldl
end

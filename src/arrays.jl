@inline eachindex_nonzeros(x::AbstractArray) = eachindex(x)
@inline values_nonzeros(x::AbstractArray) = vec(x)
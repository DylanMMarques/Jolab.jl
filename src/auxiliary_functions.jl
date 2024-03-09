reverse_if_backward(::Type{Forward}, x) = x
reverse_if_backward(::Type{Backward}, x) = reverse(x)

import Base.!
Base.length(beam::MeshedBeam) = length(beam.e)
Base.size(beam::MeshedBeam) = size(beam.e)
(!)(::Type{Forward}) = Backward
(!)(::Type{Backward}) = Forward

function check_same_definition(mode1::M1, mode2::M2) where {M2<:AbstractFieldMode{<:Any,D2}, M1<:AbstractFieldMode{<:Any,D1}} where {D1, D2}
    same_mode_type(M1, M2) || return false
    D1 == D2 || return false
    fields = fieldnames(M1)
    for field in Iterators.filter(!=(:e), fields)
        getfield(mode1, field) ≈ getfield(mode2, field) || return false
    end
    return true
end

function check_same_definition(beam1::MeshedBeam{T1,D,C}, beam2::MeshedBeam{T2,D,C}) where {T1,T2,D, C}
    beam1.mesh ≈ beam2.mesh || return false
    beam1.frame ≈ beam2.frame || return false
    beam1.medium ≈ beam2.medium || return false
    return true
end
check_same_definition(beam1::MeshedBeam, beam2::MeshedBeam) = false

same_mode_type(::Type{<:PlaneWaveScalar}, ::Type{<:PlaneWaveScalar}) = true
same_mode_type(::Type{<:PlaneWaveVectorial}, ::Type{<:PlaneWaveVectorial}) = true
same_mode_type(::Type{<:AbstractFieldMode}, ::Type{<:AbstractFieldMode}) = false

function find_local_minima(f, x::AbstractVector{T}, initial_size) where T
    local_minima = Vector{T}(undef, initial_size)
    ind_zeros = 1
    value = f(x[1])
    next_value = f(x[2])
    @inbounds for i in eachindex(x)[2:end-1]
        prev_value = value
        value = next_value
        next_value = f(x[i+1])
        if value < prev_value && value < next_value
            if ind_zeros > initial_size
                initial_size *= 2
                resize!(local_minima, initial_size)
            end
            local_minima[ind_zeros] = x[i]
            ind_zeros += 1
        end
    end
    local_minima[1:(ind_zeros-1)]
end

function overlap_integral(array1::AbstractArray, array2::AbstractArray, mesh::Domain)
    @argcheck size(array1) == size(mesh) == size(array2) DimensionMismatch
    mapreduce((i) -> array1[i] * conj(array2[i]), +, eachindex(mesh))
end

function overlap_integral(array::AbstractArray, f2::Function, mesh::Domain)
    @argcheck size(array) == size(mesh) DimensionMismatch
    mapreduce((i) -> array[i] * conj(f2(centroid(mesh, i).coords)), +, eachindex(mesh))
end

function overlap_integral(f2::Function, array::AbstractArray, mesh::Domain)
    @argcheck size(array) == size(mesh) DimensionMismatch
    mapreduce((i) -> conj(array[i]) * f2(centroid(mesh, index).coords), +, eachindex(mesh))
end

function overlap_integral(f1::Function, f2::Function, mesh::Domain)
    function f(index)
        coord = centroid(mesh, index).coords
        f1(coord) * conj(f2(coord))
    end
    mapreduce(f, +, eachindex(mesh))
end
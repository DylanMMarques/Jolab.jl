using Meshes

struct MeshData{M,D}
    mesh::M
    data::D
end

n = 100
a = CartesianGrid((n, n, n), (10.0, 20.0, 10.0), (2.0, 2.5, 3.2))
field = MeshData(a, rand(size(a)...))

_volume(x, ind) = prod(x.spacing)
function a_f(data)
    mapreduce((a, ind) -> a * _volume(data.mesh, ind), +, data.data, eachindex(data.data))
end

function _a_f(data)
    res = 0.0
    @inbounds for i in eachindex(data.data)
        res += data.data[i] * _volume(data.mesh, i)
    end
    res
end
@code_warntype a_f(field)
@time a_f(field)
@time _a_f(field)

function a_f_(a, da)
    mapreduce((a, da) -> a * da, +, a, da)
end
a_

function integral_mesh(x, data)
    mapreduce((x, data) -> area(x) * data, +, x, data)
end

function integral(data)
    mapreduce(( data) -> data, +, data)
end

function integral_meshdata(meshdata)
    mapreduce((x) -> area(x), +, eachindex(meshdata.mesh))
end

function integral_meshdata_for(meshdata)
    res = 0.0
    for f in eachindex(meshdata.mesh)
        res += centroid(meshdata.mesh, f).coords[1]
    end
    res
end

data = rand(size(a)...)
@time integral_mesh(a, data)
da = map(area, a)
_area2(a::Meshes.Quadrangle) = 2.5
mesh_data = MeshData{typeof(a), typeof(data)}(a, data)
@time integral(data)
@time integral_meshdata(mesh_data)
@time integral_meshdata_for(mesh_data)

b = Vector{Quadrangle{Float64}}(undef, length(a))

@time vertices(a)
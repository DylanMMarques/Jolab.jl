using Jolab, Test
R = 0.9
mirror = Mirror((Medium(1.0), Medium(2.0)), ReferenceFrame((0,0,0), (0,0,0)), reflectivity = R)

grid_pw2 = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, 0.1, 0.0, 500e-9)
pw2 = Jolab.AngularSpectrum(Jolab.Forward, grid_pw2, reshape([1.0 + 0.0im], 1, 1), Medium(2), ReferenceFrame((0,0,0), (0,0,0)))

grid_pw = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, 0.1, 0.0, 500e-9)
pw = Jolab.AngularSpectrum(Jolab.Forward, grid_pw, reshape([1.0 + 0.0im], 1, 1), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))

(rpw, tpw) = light_interaction(mirror, pw)
@test intensity(rpw) ≈ R * intensity(pw)
@test intensity(tpw) ≈ (1-R) * intensity(pw)
@test tpw.medium == last(mirror.mat)
@test rpw.medium == first(mirror.mat)
@test rpw.frame == tpw.frame == pw.frame

nsx = range(0, 0.5, length = 50)

grid_beam = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, 1550E-9)
beam = Jolab.AngularSpectrum(Float64, Jolab.Forward, grid_beam, nsx .* nsx', Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))
i_i = intensity(beam)
(r_beam, t_beam) = light_interaction(mirror, beam)
@test_opt light_interaction(mirror, beam)
@test intensity(r_beam) ≈ i_i * R
@test intensity(t_beam) ≈ i_i * (1-R)


mirror = Mirror((Medium(2.0), Medium(2.0)), ReferenceFrame((0,0,0), (0,0,0)), reflectivity = R)
@test_throws ArgumentError light_interaction(mirror, beam)
@test_throws ArgumentError light_interaction(mirror, beam)

mirror = Mirror((Medium(1.0), Medium(1.0)), ReferenceFrame((0,0,0), (.2,0,0)), reflectivity = R)
@test_throws ArgumentError light_interaction(mirror, beam)
@test_throws ArgumentError light_interaction(mirror, beam)


## Plane wave case

mirror = Mirror((Medium(1.0), Medium(1.0)), ReferenceFrame((0,0,0), (0,0,0)), reflectivity = R)

grid_pw_f = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, 0.1, 0.1, 500e-9)
pw_f = Jolab.AngularSpectrum(Jolab.Forward, grid_pw_f, reshape([1.0 + 0.0im], 1, 1), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))

grid_pw_b = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, 0.1, 0.1, 500e-9)
pw_b = Jolab.AngularSpectrum(Jolab.Backward, grid_pw_b, reshape([1.0 + 0.0im], 1, 1), Medium(1), ReferenceFrame((0,0,0), (0,0,0)))

(lpw_f, rpw_f) = light_interaction(mirror, pw_f)
(lpw_b, rpw_b) = light_interaction(mirror, pw_b)
@test lpw_f.e == -rpw_b.e
@test rpw_f.e == lpw_b.e

mirror = Mirror((Medium(2.0), Medium(2.0)), ReferenceFrame((0,0,0), (0,0,0)), reflectivity = R)
@test_throws ArgumentError light_interaction(mirror, pw_b)
@test_throws ArgumentError light_interaction(mirror, pw_f)

mirror = Mirror((Medium(1.0), Medium(1.0)), ReferenceFrame((0,0,0), (.2,0,0)), reflectivity = R)
@test_throws ArgumentError light_interaction(mirror, pw_b)
@test_throws ArgumentError light_interaction(mirror, pw_f)

using StaticArrays
function test_scatmat_f(mls)
    nsx = range(0, 0.95, length = 50)
    grid_beam = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, 1500E-9)
    beam = Jolab.AngularSpectrum(Float64, Jolab.Forward, grid_beam, reshape(nsx .* nsx', (length(nsx), length(nsx))), Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))
    mat = ScatteringMatrix(mls, beam)
    (mat_r, mat_t) = light_interaction(mat, beam)
    (aux_r, aux_t) = light_interaction(mls, beam)
    mat_r ≈ aux_r && mat_t ≈ aux_t
end
mls = DielectricStack(Medium.((@SVector [1, 1.5, 1])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,0)))
mirror = Mirror((Medium(1.0), Medium(1.0)), ReferenceFrame((0,0,0), (0,0,0)), reflectivity = 0.9)
@test_opt test_scatmat_f(mls)
@test test_scatmat_f(mls)
@test test_scatmat_f(mirror)

function test_scatmat_b(mls)
    nsx = range(0, 0.95, length = 50)
    grid_beam = Jolab.CartesianGrid(Jolab.NSX_NSY_λ, nsx, nsx, 1550E-9)
    beam = Jolab.AngularSpectrum(Float64, Jolab.Backward, grid_beam, nsx .* nsx', Medium(1.0), last(mls.frames))
    mat = ScatteringMatrix(mls, beam)
    (mat_r, mat_t) = mat * beam
    (aux_r, aux_t) = light_interaction(mls, beam)
    mat_r ≈ aux_r && mat_t ≈ aux_t
end
@test test_scatmat_b(mls)
@test test_scatmat_b(mirror)

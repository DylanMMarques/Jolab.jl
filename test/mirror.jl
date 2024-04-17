using Jolab, Test
R = 0.9
mirror = Mirror((Medium(1.0), Medium(2.0)), ReferenceFrame((0,0,0), (0,0,0)), reflectivity = R)
pw = MeshedPlaneWaveScalar(Forward, 0.1, 0.1, 1, 1, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
(rpw, tpw) = light_interaction(mirror, pw)
@test intensity(rpw) ≈ R * intensity(pw)
@test intensity(tpw) ≈ (1-R) * intensity(pw)
@test tpw.medium == last(mirror.mat)
@test rpw.medium == first(mirror.mat)
@test rpw.frame == tpw.frame == pw.frame

nsx = range(0, 0.95, length = 100)

beam = MonochromaticAngularSpectrum(Float64, Forward, nsx, nsx, (nsx .* nsx'), 1550E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))
i_i = intensity(beam)
(r_beam, t_beam) = light_interaction(mirror, beam)
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
pw_f = MeshedPlaneWaveScalar(Forward, 0.1, 0.1, 1, 1, Medium(1), ReferenceFrame((0,0,0), (0,0,0)))
pw_b = MeshedPlaneWaveScalar(Backward, 0.1, 0.1, 1, 1, Medium(1), ReferenceFrame((0,0,0), (0,0,0))) 
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

function test_scatmat_f(mls)
    nsx = range(0, 0.95, length = 100)
    beam = MonochromaticAngularSpectrum(Float64, Forward, nsx, nsx, (nsx .* nsx'), 1500E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)));
    mat = ScatteringMatrix(mls, beam)
    (mat_r, mat_t) = light_interaction(mat, beam)
    (aux_r, aux_t) = light_interaction(mls, beam)
    mat_r ≈ aux_r && mat_t ≈ aux_t
end
mls = DielectricStack(Medium.((@SVector [1, 1.5, 1])), (@SVector [100E-9]), ReferenceFrame((0,0,0), (0,0,0)))
mirror = Mirror((Medium(1.0), Medium(1.0)), ReferenceFrame((0,0,0), (0,0,0)), reflectivity = 0.9)
@test test_scatmat_f(mls)
@test test_scatmat_f(mirror)

function test_scatmat_b(mls)
    nsx = range(0, 0.95, length = 100)
    beam = MonochromaticAngularSpectrum(Float64, Backward, nsx, nsx, (nsx .* nsx'), 1550E-9, Medium(1.0), last(mls.frames));
    mat = ScatteringMatrix(mls, beam)
    (mat_r, mat_t) = mat * beam
    (aux_r, aux_t) = light_interaction(mls, beam)
    mat_r ≈ aux_r && mat_t ≈ aux_t
end
@test test_scatmat_b(mls)
@test test_scatmat_b(mirror)
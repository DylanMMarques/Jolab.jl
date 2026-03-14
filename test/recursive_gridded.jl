using Jolab, Test

nsx = range(-0.1, 0.1, length=64)
nsy = range(-0.1, 0.1, length=64)
λ = 1550E-9
beam = MonochromaticAngularSpectrum_gaussian(Float64, Forward, nsx, nsx, 50E-6, λ, Medium(1.0), ReferenceFrame((0,0,0.0), (0,0,0.0)));

mls = DielectricStack(Medium.([1, 1.5, 2]), [100E-9], ReferenceFrame((0,0,0.0), (0,0,0.0)))
(r, t) = light_interaction(mls, beam)

prop = Propagation(mls.frames, Medium(1.5))
mls_layer = (DielectricStack(Medium.([1, 1.5]), zeros(0), ReferenceFrame((0,0,0.0), (0,0,0.0))), prop, DielectricStack(Medium.([1.5, 2]), zeros(0), ReferenceFrame((0,0,100E-9), (0,0,0.0))))

(r2, t2) = Jolab.lightinteraction_recursivegridded(mls_layer, beam, rtol = 1E-15, printBool = false)
@test isapprox(r2, r, rtol = 1E-7)
@test isapprox(t2, t, rtol = 1E-7)

mls = DielectricStack(Medium.([1, 1.5, 2, 2.5, 1]), [100E-9, 500E-9, 200E-9], ReferenceFrame((0,0,0.0), (0,0,0.0)))
(r, t) = light_interaction(mls, beam)

z(x,y) = 0.0
z1 = zeros(0)

a = true
mls_layer = (Jolab.DielectricStack(Medium.([1, 1.5]), z1, ReferenceFrame((0,0,0.0), (0,0,0.0))), 
    Propagation((ReferenceFrame((0,0,0.0), (0,0,0.0)), ReferenceFrame((0,0,100E-9), (0,0,0.0))), Medium(1.5)),
    Jolab.DielectricStack(Medium.([1.5, 2]), z1, ReferenceFrame((0,0,100E-9), (0,0,0.0))),
    Propagation((ReferenceFrame((0,0,100E-9), (0,0,0.0)), ReferenceFrame((0,0,600E-9), (0,0,0.0))), Medium(2.0)),
    Jolab.DielectricStack(Medium.([2, 2.5]), z1, ReferenceFrame((0,0,600E-9), (0,0,0.0))),
    Propagation((ReferenceFrame((0,0,600E-9), (0,0,0.0)), ReferenceFrame((0,0,800E-9), (0,0,0.0))), Medium(2.5)),
    Jolab.DielectricStack(Medium.([2.5, 1]), z1, ReferenceFrame((0,0,800E-9), (0,0,0.0)))
)

(r2, t2) = Jolab.lightinteraction_recursivegridded(mls_layer, beam, rtol = 1E-20, printBool = false)
@test isapprox(r2, r, rtol = 1E-7)
@test isapprox(t2, t, rtol = 1E-7)

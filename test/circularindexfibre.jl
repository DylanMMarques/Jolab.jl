using Jolab, Test, Enzyme, Optim

profile = CircularStepIndexProfile(100E-6, 0.2, Medium(1.55))

fibre = Fibre(profile, 1, Medium.((1,1)), (ReferenceFrame((0,0,0), (0,0,0)), ReferenceFrame((0,0,1), (0,0,0))))

Jolab.findmodes!(fibre, 1500E-9)
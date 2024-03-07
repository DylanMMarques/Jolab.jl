using Jolab, Test
nsx = range(-.5, 0.5, length = 100)
beam = MonochromaticAngularSpectrum(Float64, Forward, nsx, nsx, (nsx .* nsx'), 1550E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))

# @test intensity(beam) ≈ 1

## Check same definition
@test Jolab.check_same_definition(beam, beam)
@test (beam ≈ deepcopy(beam))

beam2 = MonochromaticAngularSpectrum(Float64, Forward, 2nsx, 2nsx, (nsx .* nsx'), 1550E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))
@test !Jolab.check_same_definition(beam, beam2)
@test !(beam2 ≈ beam)

beam2 = MonochromaticAngularSpectrum(Float64, Forward, 2nsx, 2nsx, (nsx .* nsx'), 1500E-9, Medium(1.0), ReferenceFrame((0,0,0), (0,0,0)))
@test !Jolab.check_same_definition(beam, beam2)
@test !(beam2 ≈ beam)

beam2 = MonochromaticAngularSpectrum(Float64, Forward, 2nsx, 2nsx, (nsx .* nsx'), 1500E-9, Medium(2.0), ReferenceFrame((0,0,0), (0,0,0)))
@test !Jolab.check_same_definition(beam, beam2)
@test !(beam2 ≈ beam)

beam2 = MonochromaticAngularSpectrum(Float64, Forward, 2nsx, 2nsx, (nsx .* nsx'), 1500E-9, Medium(2.0), ReferenceFrame((0,0,0), (0,0,1)))
@test !Jolab.check_same_definition(beam, beam2)
@test !(beam2 ≈ beam)

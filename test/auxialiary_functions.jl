function test_cartesian_with_radialsymmetric(ω, comp, medium, sampling_rad, sampling_cart; rtol = 1E-9)
    nsr = range(0, 0.5, length = sampling_rad)
    radial_i = MonochromaticAngularSpectrumRadialSymmetric_gaussian(Forward, nsr, ω, 1550E-9, medium, ReferenceFrame((0,0,0), (0,0,0)))
    nsx = range(-0.5, 0.5, length = sampling_cart)
    cart_i = MonochromaticAngularSpectrum_gaussian(Forward, nsx, nsx, ω, 1550E-9, medium, ReferenceFrame((0,0,0), (0,0,0)))
    

    (b_radial, f_radial) = light_interaction(comp, radial_i)
    (b_cart, f_cart) = light_interaction(comp, cart_i)
    isapprox(intensity(b_radial), intensity(b_cart); rtol = rtol) && isapprox(intensity(f_radial), intensity(f_cart), rtol = rtol)
end

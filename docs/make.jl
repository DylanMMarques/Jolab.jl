using Documenter, Jolab

makedocs(sitename="Jolab.jl",
        doctest=false, clean=true,
        authors="Dylan Marques",
        # pages = Any[
         # "Home" => "index.md",
         # "Sequential Modelling" => Any[
         #   "exampleComponent.md"
         # ],
         # "Angular Spectrum" => Any["exampleAngularSpectrum.md"],
         # "Fabry-Pérot etalon" => Any["exampleAiry.md"],
         # "Multimode fibre" => Any["exampleMultimodeFibre.md"]]
         )

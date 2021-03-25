using Documenter, Jolab

makedocs(sitename="Jolab.jl",
        doctest=false, clean=true,
        authors="Dylan Marques",
        pages = Any[
         "Home" => "index.md",
         "Optical components" => Any["MultilayerStructure.md", "MultilayerStructure.md"],
         "Examples" => Any["exampleComponent.md",
         "exampleAngularSpectrum.md",
         "exampleAiry.md",
         "exampleMultimodeFibre.md"]
         ]
        )

deploydocs(
    repo = "github.com/DylanMMarques/Jolab.jl.git",
)

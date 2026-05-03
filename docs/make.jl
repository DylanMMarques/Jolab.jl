using Documenter, DocumenterVitepress

using Jolab

makedocs(;
    modules = [Jolab],
    authors = "Dylan M. Marques",
    repo = "https://github.com/DylanMMarques/Jolab.jl",
    sitename = "Jolab.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/DylanMMarques/Jolab.jl",
        devurl = "dev",
        deploy_url = "DylanMMarques.github.io/Jolab.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started/first_simulation.md",
        "Examples" => "examples/microscopy.md",
    ],
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/DylanMMarques/Jolab.jl",
    push_preview = true,
)

using jlmie
using Documenter

makedocs(;
    modules=[jlmie],
    authors="Tatsuki Hinamoto",
    repo="https://github.com/Hinamoooon/jlmie.jl/blob/{commit}{path}#L{line}",
    sitename="jlmie.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Hinamoooon.github.io/jlmie.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Hinamoooon/jlmie.jl",
)

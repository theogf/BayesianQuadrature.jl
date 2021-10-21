using BayesianQuadrature
using Documenter

makedocs(;
    modules=[BayesianQuadrature],
    authors="Theo Galy-Fajou <theo.galyfajou@gmail.com> and contributors",
    repo="https://github.com/theogf/BayesianQuadrature.jl/blob/{commit}{path}#L{line}",
    sitename="BayesianQuadrature.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://theogf.github.io/BayesianQuadrature.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "User Guide" => "userguide.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/theogf/BayesianQuadrature.jl", devbranch="main", push_preview=true
)

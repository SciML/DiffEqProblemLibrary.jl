using Documenter
using DiffEqProblemLibrary
using BVProblemLibrary, DAEProblemLibrary, DDEProblemLibrary, JumpProblemLibrary
using NonlinearProblemLibrary, ODEProblemLibrary, SDEProblemLibrary

makedocs(;
    modules = [
        DiffEqProblemLibrary,
        BVProblemLibrary, DAEProblemLibrary, DDEProblemLibrary, JumpProblemLibrary,
        NonlinearProblemLibrary, ODEProblemLibrary, SDEProblemLibrary,
    ],
    authors = "SciML Contributors",
    sitename = "DiffEqProblemLibrary.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://docs.sciml.ai/DiffEqProblemLibrary/stable/",
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/SciML/DiffEqProblemLibrary.jl.git",
)

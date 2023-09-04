push!(LOAD_PATH, joinpath(@__DIR__,"../src/"))
using Documenter
using CatmapInterface
using Catalyst

function mkdocs()
    # generate html for notebooks

    DocMeta.setdocmeta!(CatmapInterface, :DocTestSetup, :(using CatmapInterface, Catalyst); recursive=true)

    makedocs(
        sitename    = "CatmapInterface.jl",
        modules     = [CatmapInterface],
        format      = Documenter.HTML(mathengine=MathJax3()),
        clean       = false,
        doctest     = true,
        draft       = false,
        authors     = "S. Maaß",
        repo        = "https://github.com/smaasz/CatmapInterface.jl",
        pages       = [
            "Home"  => "index.md",
            "API"   => "API.md"
        ]
    )
end

mkdocs()
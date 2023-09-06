push!(LOAD_PATH, joinpath(@__DIR__,"../src/"))
using Documenter
using CatmapInterface
using Catalyst
using PlutoStaticHTML
using Pkg

const NOTEBOOK_DIR  = joinpath(@__DIR__, "..", "notebooks") 
const NOTEBOOKS     = ["CO2R"]
const NOTEBOOKS_JL  = NOTEBOOKS .* ".jl"
const NOTEBOOKS_MD  = NOTEBOOKS .* ".md"

function build_all_notebooks()
    thisdir=pwd()
    Pkg.activate(NOTEBOOK_DIR)
    Pkg.develop(PackageSpec(path=pwd()))
    Pkg.instantiate()
    Pkg.activate(thisdir)
    println("Building notebooks in $NOTEBOOK_DIR")
    ENV["PLUTO_PROJECT"]=NOTEBOOK_DIR
    oopts = OutputOptions(; append_build_context=true)
    output_format = documenter_output
    bopts = BuildOptions(NOTEBOOK_DIR; output_format)
    build_notebooks(bopts,NOTEBOOKS_JL, oopts)
    return nothing
end

function mkdocs()
    # generate html for notebooks
    notebook_md_dir  = joinpath(@__DIR__,"src","notebooks")
    rm(notebook_md_dir,force=true,recursive=true)
    mkdir(notebook_md_dir)
    build_all_notebooks()
    for nb in NOTEBOOKS_MD
        mv(joinpath(NOTEBOOK_DIR,nb),joinpath(notebook_md_dir,nb))
    end
    notebooks=joinpath.("notebooks",NOTEBOOKS_MD)

    notebooks=[ nb*".jl"=> joinpath("notebooks",nb*".md") for nb in NOTEBOOKS ]

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
            "Home"      => "index.md",
            "Guide"     => "guide.md",
            "Public"    => "public.md",
            "Internal"  => "internal.md",
            "Notebooks" => notebooks,
        ]
    )
end

mkdocs()

if !isinteractive()
    deploydocs(repo = "github.com/smaasz/CatmapInterface.jl")
end
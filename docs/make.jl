using Documenter, SimplexGridFactory, ExtendableGrids
import PlutoSliderServer
using GridVisualize, ExampleJuggler
import CairoMakie
ExampleJuggler.verbose!(true)


###zzaccessing

function mkdocs()
    cleanexamples()
    exampledir = joinpath(@__DIR__, "..", "examples")
    notebookdir = joinpath(@__DIR__, "..", "notebooks")
    cairo_examples = @docscripts(exampledir, ["examples2d.jl"], Plotter = CairoMakie)
    cairo_3d_examples = @docscripts(exampledir, ["examples3d.jl"], Plotter = CairoMakie)

    generated_examples = [cairo_examples..., cairo_3d_examples...]
    notebook_examples = @docplutonotebooks(notebookdir, ["gridgenvis.jl", "cylinder.jl"], iframe = true, iframe_height = "2000px")


    makedocs(;
        sitename = "SimplexGridFactory.jl",
        modules = [SimplexGridFactory],
        doctest = false,
        clean = false,
        authors = "J. Fuhrmann, Ch. Merdon",
        repo = "https://github.com/WIAS-PDELib/SimplexGridFactory.jl",
        pages = [
            "Home" => "index.md",
            "Changes" => "changes.md",
            "API" => "api.md",
            "Examples" => generated_examples,
            "Notebooks" => notebook_examples,
            "Internals" => "internals.md",
            "allindex.md",
        ]
    )

    return cleanexamples()

end

mkdocs()

deploydocs(; repo = "github.com/WIAS-PDELib/SimplexGridFactory.jl.git")

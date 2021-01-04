push!(LOAD_PATH,"../src/")
using Documenter, ComputeIBZ

makedocs(sitename="ComputeIBZ",
         # format = Documenter.HTML(prettyurls = false),
         modules = [ComputeIBZ],
         authors = "Jeremy Jorgensen",
         doctest = true,
         pages=["index.md", "Documentation.md", "Usage.md"])

deploydocs(
    repo = "github.com/jerjorg/ComputeIBZ.jl.git",
    versions = ["stable" => "v#.#"],
    devurl = "docs")

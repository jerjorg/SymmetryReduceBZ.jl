push!(LOAD_PATH,"../src/")
using Documenter, SymmetryReduceBZ

makedocs(sitename="SymmetryReduceBZ",
         # format = Documenter.HTML(prettyurls = false),
         modules = [SymmetryReduceBZ],
         authors = "Jeremy Jorgensen",
         doctest = true,
         pages=["index.md", "Documentation.md", "Usage.md"])

deploydocs(
    repo = "github.com/jerjorg/SymmetryReduceBZ.jl.git",
    versions = ["stable" => "v#.#"],
    devurl = "docs")

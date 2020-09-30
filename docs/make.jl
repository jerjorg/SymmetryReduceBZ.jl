push!(LOAD_PATH,"../src/")
using Documenter, IBZ

makedocs(sitename="IBZ",
         # format = Documenter.HTML(prettyurls = false),
         modules = [IBZ],
         authors = "Jeremy Jorgensen",
         doctest = true,
         pages=["index.md", "Documentation.md", "Usage.md"])

deploydocs(
    repo = "github.com/jerjorg/IBZ.jl.git",
    versions = ["stable" => "v#.#"],
    devurl = "docs")

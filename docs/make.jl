# Generate documentation with this command:
# (cd docs && julia make.jl)

push!(LOAD_PATH, "..")

using Documenter
using GraviPet

makedocs(;
    authors="Erik Schnetter",
    format=Documenter.HTML(),
    modules=[GraviPet],
    pages=["index.md", "Domains" => "domains.md", "Categories" => "categories.md", "Helpers" => "helpers.md"],
    sitename="GraviPet",
)

deploydocs(; devbranch="main", push_preview=true, repo="github.com/eschnett/GraviPet.jl.git")

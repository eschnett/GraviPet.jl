# Generate documentation with this command:
# (cd docs && julia make.jl)

push!(LOAD_PATH, "..")

using Documenter
using GraviPet

makedocs(; sitename="GraviPet", format=Documenter.HTML(), modules=[GraviPet])

deploydocs(; repo="github.com/eschnett/GraviPet.jl.git", devbranch="main", push_preview=true)

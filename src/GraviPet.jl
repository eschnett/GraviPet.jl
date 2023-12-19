module GraviPet

using Reexport

# internal definitions
include("Defs.jl")

# public helpers
include("Funs.jl")

# abstract types
include("Domains.jl")
include("Categories.jl")

# concrete types for domains
include("Intervals.jl")
include("Boxes.jl")

# concrete types for functions
include("BlockFunctions.jl")
include("GridFunctions.jl")
include("GridFunctions1D.jl")
include("JuliaFunctions.jl")
include("ThreadedFunctions.jl")

@reexport using ..Funs

@reexport using ..Domains
@reexport using ..Categories

@reexport using ..Boxes
@reexport using ..Intervals

@reexport using ..BlockFunctions
@reexport using ..GridFunctions
@reexport using ..GridFunctions1D
@reexport using ..JuliaFunctions
@reexport using ..ThreadedFunctions

end

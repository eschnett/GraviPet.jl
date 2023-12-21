using Base.Threads
using Distributed
using Hwloc

addprocs(num_physical_cores())
println("Running tests:")
println("    - $(nworkers()) worker processes")
println("    - $(nthreads()) threads")
println("    - $(num_physical_cores()) cores")

@everywhere using DoubleFloats
@everywhere using GraviPet
@everywhere using KernelAbstractions
@everywhere using Random
@everywhere using StaticArrays
@everywhere using Test

################################################################################

@everywhere const BigRat = Rational{BigInt}

@everywhere function random_SVector(xmin::SVector{D,T}, xstep::SVector{D,T}, xmax::SVector{D,T}) where {D,T}
    return SVector{D,T}(rand(xmin[d]:xstep[d]:xmax[d]) for d in 1:D)
end

################################################################################

# Define tests
include("test_Domains.jl")
include("test_Categories.jl")

# Run domain tests
include("runtests_Intervals.jl")
include("runtests_Boxes.jl")

# Run category tests
include("runtests_GridFunctions1D.jl")
include("runtests_GridFunctions.jl")
include("runtests_BlockFunctions.jl")
include("runtests_ThreadedFunctions.jl")
include("runtests_DistributedFunctions.jl")
include("runtests_KernelFunctions.jl")

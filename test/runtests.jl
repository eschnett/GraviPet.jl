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
@everywhere using Random
@everywhere using StaticArrays
@everywhere using Test

################################################################################

@everywhere const BigRat = Rational{BigInt}

@everywhere function random_SVector(xmin::SVector{D,T}, xstep::SVector{D,T}, xmax::SVector{D,T}) where {D,T}
    return SVector{D,T}(rand(xmin[d]:xstep[d]:xmax[d]) for d in 1:D)
end

################################################################################

include("test_Domain.jl")
include("test_Category.jl")

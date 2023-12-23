using Base.Threads
using Distributed
using Hwloc

addprocs(num_physical_cores())
println("Running tests:")
println("    - $(nworkers()) worker processes")
println("    - $(nthreads()) threads")
println("    - $(num_physical_cores()) cores")

@everywhere using CUDA
@everywhere using DoubleFloats
@everywhere using GraviPet
@everywhere using KernelAbstractions
@everywhere using Random
@everywhere using StaticArrays
@everywhere using Test

@static if Sys.isapple() && VERSION >= v"1.8"
    # Metal requires macOS and at least Julia 1.8
    using Pkg
    Pkg.add("Metal")
end
@everywhere @static if Sys.isapple() && VERSION >= v"1.8"
    # Metal requires macOS and at least Julia 1.8
    using Metal
end

################################################################################

@everywhere const BigRat = Rational{BigInt}

@everywhere function random_SVector(xmin::SVector{D,T}, xstep::SVector{D,T}, xmax::SVector{D,T}) where {D,T}
    return SVector{D,T}(rand(xmin[d]:xstep[d]:xmax[d]) for d in 1:D)
end

################################################################################

alltests = [
    # Domain tests
    "Intervals",
    "Boxes",
    # Category tests
    "GridFunctions1D",
    "GridFunctions",
    "BlockFunctions",
    "ThreadedFunctions",
    "DistributedFunctions",
    "KernelFunctions_CPU",
    "KernelFunctions_CUDA",
    "KernelFunctions_Metal",
]
chosentests = split(get(ENV, "TESTS", "all"))
tests = []
for test in chosentests
    if test == "all"
        append!(tests, alltests)
    elseif test in alltests
        push!(tests, test)
    else
        error("Unknown test \"$test\"")
    end
end
unique!(tests)
println("Selected tests:")
for test in tests
    println("    - $test")
end
if isempty(tests)
    println("No tests selected, doing nothing")
end

################################################################################

# Define tests
include("test_Domains.jl")
include("test_Categories.jl")

println("Running tests:")
for test in tests
    include("runtests_$test.jl")
end
println("Done.")

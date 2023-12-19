using Base.Threads
using DoubleFloats
using GraviPet
using Random
using StaticArrays
using Test

################################################################################

const BigRat = Rational{BigInt}

function random_SVector(xmin::SVector{D,T}, xstep::SVector{D,T}, xmax::SVector{D,T}) where {D,T}
    return SVector{D,T}(rand(xmin[d]:xstep[d]:xmax[d]) for d in 1:D)
end

################################################################################

include("test_Domain.jl")
include("test_Category.jl")

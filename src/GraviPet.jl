module GraviPet

using FillArrays
using HCubature
using SparseArrays
using SparseArraysCOO
using StaticArrays

function lincom(x0::S, y0::T, x1::S, y1::T, x::S) where {S,T}
    return T(y0 .* ((x - x1) ./ (x0 - x1)) + y1 .* ((x - x0) ./ (x1 - x0)))
end

include("Domain.jl")
include("Category.jl")

end

module GridFunctions1D

using ..Categories
using ..Defs
using ..Intervals

struct GridFunction1D{S,T} <: Category{Interval{S},Interval{T}}
    name::AbstractString
    domain::Interval{S}
    codomain::Interval{T}
    grid::AbstractVector{T}
    function GridFunction1D{S,T}(
        name::AbstractString, domain::Interval{S}, codomain::Interval{T}, grid::AbstractVector{T}
    ) where {S,T}
        length(grid) >= 2 ||
            throw(ArgumentError("GridFunction1D length must be at least 2 so that linear functions can be represented"))
        return new{S,T}(name, domain, codomain, grid)
    end
end
function GridFunction1D(name::AbstractString, domain::Interval{S}, codomain::Interval{T}, grid::AbstractVector{T}) where {S,T}
    return GridFunction1D{S,T}(name, domain, codomain, grid)
end
export GridFunction1D

"""
    GridFunction1D{S,T}(name::AbstractString, domain::Interval{S}, codomain::Interval{T}, npoints::Int)
    GridFunction1D(name::AbstractString, domain::Interval{S}, codomain::Interval{T}, npoints::Int)

Create a zero-valued grid function, i.e. a grid function that is zero
everywhere. The `codomain` can be a zero-valued domain created via
[`Interval{T}()`](@ref).

`npoints` specifies the resolution of the grid function. This makes
it convenient to use zero-valued grid functions as "templates" or
"skeletons" when creating other grid functions.

This zero-valued grid function is represented efficiently and does
*not* store one element per grid point.
"""
function GridFunction1D{S,T}(name::AbstractString, domain::Interval{S}, codomain::Interval{T}, npoints::Int) where {S,T}
    return GridFunction1D(name, domain, Interval{T}(), Zeros{T}(npoints))
end
function GridFunction1D(name::AbstractString, domain::Interval{S}, codomain::Interval{T}, npoints::Int) where {S,T}
    return GridFunction1D{S,T}(name, domain, codomain, npoints)
end

# Metadata
Categories.name(gf::GridFunction1D) = gf.name
Categories.domain(gf::GridFunction1D) = gf.domain
Categories.codomain(gf::GridFunction1D) = gf.codomain

# Equality
function Base.:(==)(gf1::GridFunction1D, gf2::GridFunction1D)
    domain(gf1) == domain(gf2) || return false
    codomain(gf1) == codomain(gf2) || return false
    return gf1.grid == gf2.grid
end

# Collection
Base.isempty(x::GridFunction1D) = isempty(x.grid)
Base.length(x::GridFunction1D) = length(x.grid)
Base.getindex(x::GridFunction1D, i) = x.grid[i]
function Base.map(f, x::GridFunction1D)
    grid′ = map(f, gf.grid)
    return GridFunction1D(x.name, x.domain, Interval(extrema(grid′)...), grid′)
end
function Base.map(f, x::GridFunction1D, y::GridFunction1D)
    @assert domain(x) == domain(y)
    @assert axes(x.grid) == axes(y.grid)
    grid′ = map(f, x.grid, y.grid)
    return GridFunction1D(x.name, x.domain, Interval(extrema(grid′)...), grid′)
end

# Category
function Categories.make_identity(gf::GridFunction1D)
    dom = gf.domain
    T = eltype(dom)
    @assert length(gf.grid) >= 2
    grid = if isempty(gf.grid)
        T[]
    elseif length(gf.grid) == 1
        T[first(dom)]
    else
        range(T(first(dom)), T(last(dom)); length=length(gf.grid))
    end
    return GridFunction1D(gf.name, dom, dom, grid)
end
function Base.:∘(gf2::GridFunction1D, gf1::GridFunction1D)
    @assert domain(gf2) == codomain(gf1)
    return map(gf2, gf1)
end

# Vector space
Base.zero(x::GridFunction1D) = map(zero, x)
Base.:+(x::GridFunction1D) = map(+, x)
Base.:-(x::GridFunction1D) = map(-, x)
Base.:+(x::GridFunction1D, y::GridFunction1D) = map(+, x, y)
Base.:-(x::GridFunction1D, y::GridFunction1D) = map(-, x, y)
Base.:*(a, x::GridFunction1D) = map(w -> a * w, x)
Base.:*(x::GridFunction1D, a) = map(w -> w * a, x)
Base.:\(a, x::GridFunction1D) = map(w -> a \ w, x)
Base.:/(x::GridFunction1D, a) = map(w -> w / a, x)

# We sample instead of projecting...
function Categories.project(gf::GridFunction1D{S}, cat::Category) where {S}
    @assert domain(gf) == domain(cat)
    id = make_identity(gf)
    return map(cat, id)::GridFunction1D{S}
end

function Categories.evaluate(gf::GridFunction1D{<:Any,T}, x::T) where {T}
    isempty(gf.grid) && return zero(T)::T
    length(gf.grid) == 1 && return gf.grid[begin]::T

    dom = domain(gf)
    S = eltype(dom)
    ix = lincom(first(dom), S(firstindex(gf.grid)), last(dom), S(lastindex(gf.grid)), x)
    i = clamp(floor(Int, ix), firstindex(gf.grid):(lastindex(gf.grid) - 1))
    di = ix - i

    fx = lincom(S(0), gf.grid[i], S(1), gf.grid[i + 1], di)
    return fx::T
end
Categories.evaluate(gf::GridFunction1D{<:Any,T}, x) where {T} = evaluate(gf, T(x))

function Categories.integrate(gf::GridFunction1D)
    dom = domain(gf)
    cod = codomain(gf)
    S = eltype(dom)
    T = eltype(cod)
    imin = first(axes(gf.grid)[1])
    imax = last(axes(gf.grid)[1])
    h = (last(dom) - first(dom)) / (length(gf.grid) - 1)
    s = zero(T)
    for i in imin:imax
        w = i == imin || i == imax ? one(S) / 2 : one(S)
        s += T(w * gf.grid[i])
    end
    return T(h * s)
end

end

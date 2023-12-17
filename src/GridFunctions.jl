module GridFunctions

using FillArrays
using HCubature
using SparseArrays
using SparseArraysCOO
using StaticArrays

using ..Boxes
using ..Categories
using ..Defs

struct GridFunction{DS,S,DT,T} <: Category{Box{DS,S},Box{DT,T}}
    name::AbstractString
    domain::Box{DS,S}
    codomain::Box{DT,T}
    grid::AbstractArray{SVector{DT,T},DS}
    function GridFunction{DS,S,DT,T}(
        name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid::AbstractArray{SVector{DT,T},DS}
    ) where {DS,S,DT,T}
        all(size(grid) .>= 2) || throw(
            ArgumentError("GridFunction size must be at least 2 in each dimension so that linear functions can be represented")
        )
        return new{DS,S,DT,T}(name, domain, codomain, grid)
    end
end
function GridFunction(
    name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid::AbstractArray{SVector{DT,T},DS}
) where {DS,S,DT,T}
    return GridFunction{DS,S,DT,T}(name, domain, codomain, grid)
end
export GridFunction

"""
    GridFunction{DS,S,DT,T}(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_size::SVector{DS,Int})
    GridFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_size::SVector{DS,Int})

Create a zero-valued grid function, i.e. a grid function that is zero
everywhere. The `codomain` can be a zero-valued domain created via
[`Box{DT,T}()`](@ref).

`grid_size` specifies the resolution of the grid function. This makes
it convenient to use zero-valued grid functions as "templates" or
"skeletons" when creating other grid functions.

This zero-valued grid function is represented efficiently and does
*not* store one element per grid point.
"""
function GridFunction{DS,S,DT,T}(
    name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_size::SVector{DS,Int}
) where {DS,S,DT,T}
    return GridFunction(name, domain, Box{DT,T}(), Zeros{SVector{DT,T}}(Tuple(grid_size))::AbstractArray{SVector{DT,T},DS})
end
function GridFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_size::SVector{DS,Int}) where {DS,S,DT,T}
    return GridFunction{DS,S,DT,T}(name, domain, codomain, grid_size)
end

# Metadata
Categories.name(gf::GridFunction) = gf.name
Categories.domain(gf::GridFunction) = gf.domain
Categories.codomain(gf::GridFunction) = gf.codomain

# Equality
function Base.:(==)(gf1::GridFunction, gf2::GridFunction)
    domain(gf1) == domain(gf2) || return false
    codomain(gf1) == codomain(gf2) || return false
    return gf1.grid == gf2.grid
end

# Collection
Base.isempty(x::GridFunction) = isempty(x.grid)
Base.length(x::GridFunction) = length(x.grid)
Base.getindex(x::GridFunction, i) = x.grid[i]
Base.getindex(x::GridFunction{DS}, i::SVector{DS}) where {DS} = x.grid[i...]
function Base.map(f, x::GridFunction)
    grid′ = map(f, x.grid)
    VT = eltype(grid′)
    @assert VT <: SVector
    cod = Box(extrema(grid′)...)
    return GridFunction(x.name, x.domain, cod, grid′)
end
function Base.map(f, x::GridFunction{DS}, y::GridFunction{DS}) where {DS}
    @assert domain(x) == domain(y)
    @assert axes(x.grid) == axes(y.grid)
    grid′ = map(f, x.grid, y.grid)
    VT = eltype(grid′)
    @assert VT <: SVector
    cod = Box(extrema(grid′)...)
    return GridFunction1D(x.name, x.domain, cod, grid′)
end

# Category
function Categories.make_identity(gf::GridFunction)
    dom = gf.domain
    VS = eltype(dom)::Type{<:SVector}
    DS = ndims(dom)
    xmin = first(dom)
    xmax = last(dom)
    imin = SVector{DS,Int}(first.(axes(gf.grid)))
    imax = SVector{DS,Int}(last.(axes(gf.grid)))
    grid = VS[lincom(VS(imin), xmin, VS(imax), xmax, VS(Tuple(i))) for i in CartesianIndex(Tuple(imin)):CartesianIndex(Tuple(imax))]
    return GridFunction(gf.name, dom, dom, grid)
end
function Base.:∘(gf2::GridFunction{<:Any,<:Any,DT}, gf1::GridFunction{DT}) where {DT}
    @assert domain(gf2) == codomain(gf1)
    return map(gf2, gf1)
end

# Vector space
Base.zero(x::GridFunction) = map(zero, x)
Base.:+(x::GridFunction) = map(+, x)
Base.:-(x::GridFunction) = map(-, x)
Base.:+(x::GridFunction, y::GridFunction) = map(+, x, y)
Base.:-(x::GridFunction, y::GridFunction) = map(-, x, y)
Base.:*(a, x::GridFunction) = map(w -> a * w, x)
Base.:*(x::GridFunction, a) = map(w -> w * a, x)
Base.:\(a, x::GridFunction) = map(w -> a \ w, x)
Base.:/(x::GridFunction, a) = map(w -> w / a, x)

function Categories.project(gf::GridFunction, cat::Category)
    @assert domain(gf) == domain(cat)
    dom = domain(cat)
    cod = codomain(cat)
    VS = eltype(dom)
    DS = ndims(dom)
    S = eltype(VS)
    VT = eltype(cod)
    DT = ndims(cod)
    T = eltype(VT)

    imin = SVector{DS,Int}(first.(axes(gf.grid)))
    imax = SVector{DS,Int}(last.(axes(gf.grid)))

    Mshape = size(gf.grid)
    Msize = prod(Mshape)
    @assert all(first.(axes(gf.grid)) .== 1)
    function Mindex(i)
        off = 0
        str = 1
        for d in 1:DS
            off += str * (i[d] - 1)
            str *= Mshape[d]
        end
        return off + 1
    end
    M = SparseMatrixCOO{S}(Msize, Msize)

    for i0 in CartesianIndices(axes(gf.grid))
        i = SVector{DS,Int}(Tuple(i0))
        for di0 in CartesianIndices(ntuple(d -> -1:+1, DS))
            di = SVector{DS,Int}(Tuple(di0))
            j = i + di
            if all(j .∈ axes(gf.grid))
                Mvalue = prod(ifelse(di[d] == 0, one(S) / 3, one(S) / 6) for d in 1:DS; init=one(S))
                M[Mindex(i), Mindex(j)] = Mvalue
            end
        end
    end

    M = sparse(M)
    @show M

    grid = VT[
        let
            i = SVector{DS,Int}(Tuple(i0))
            s = zero(VT)
            for di0 in CartesianIndices(ntuple(d -> -1:+1, DS))
                di = SVector{DS,Int}(Tuple(di0))
                j = i + di
                if all(j .∈ axes(gf.grid))
                    Mvalue = prod(ifelse(di[d] == 0, one(S) / 3, one(S) / 6) for d in 1:DS; init=one(S))
                    M[Mindex(i), Mindex(j)] = Mvalue
                end
            end
            for di0 in CartesianIndices(ntuple(d -> -1:0, DS))
                di = SVector{DS,Int}(Tuple(di0))
                j = i + di
                if all((j .>= imin) .& (j .+ 1 .<= imax))
                    x0 = lincom(imin, first(dom), imax, last(dom), j .+ 0)
                    x1 = lincom(imin, first(dom), imax, last(dom), j .+ 1)

                    y0 = SVector{DS,S}(di[d] + 0 for d in 1:DS)
                    y1 = SVector{DS,S}(1 - di[d] for d in 1:DS)
                    basis(x) = prod(lincom(x0[d], y0[d], x1[d], y1[d], x[d]) for d in 1:DS; init=one(S))::S

                    s += hcubature(x -> VT(basis(x) * cat(x)), x0, x1)[1]::VT
                end
            end
            s
        end for i0 in CartesianIndices(axes(gf.grid))
    ]
    grid::Array{VT,DS}

    # grid = reshape(M \ reshape(grid, :), size(grid))
    @assert DT == 1
    grid = reinterpret(S, grid)
    grid = reshape(M \ reshape(grid, :), size(grid))
    grid = reinterpret(VS, grid)

    return GridFunction{DS,S,DT,T}(cat.name, dom, cod, grid)
end

function Categories.evaluate(gf::GridFunction{DS,S,DT,T}, x::SVector{DS,S}) where {DS,S,DT,T}
    isempty(gf.grid) && return zero(SVector{DT,T})::SVector{DT,T}
    length(gf.grid) == 1 && return gf.grid[begin]::SVector{DT,T}

    # x = first <=> i = firstindex
    # x = last  <=> i = lastindex
    dom = domain(gf)
    VS = eltype(dom)
    imin = SVector{DS,Int}(first.(axes(gf.grid)))
    imax = SVector{DS,Int}(last.(axes(gf.grid)))
    ix = lincom(first(dom), VS(imin), last(dom), VS(imax), x)
    i = SVector{DS,Int}(clamp(floor(Int, ix[d]), imin[d]:(imax[d] - 1)) for d in 1:DS)
    q = (ix - i)::VS

    fx = zero(SVector{DT,T})
    for di0 in CartesianIndex(ntuple(d -> 0, DS)):CartesianIndex(ntuple(d -> 1, DS))
        di = SVector{DS,Int}(Tuple(di0))
        w = one(S)
        for d in 1:DS
            w *= di[d] == 0 ? 1 - q[d] : q[d]
        end
        fx += SVector{DT,T}(w * gf.grid[CartesianIndex(Tuple(i + di))])
    end
    return fx::SVector{DT,T}
end
Categories.evaluate(gf::GridFunction{DS,S}, x::SVector{DS}) where {DS,S} = evaluate(gf, SVector{DS,S}(x))

function Categories.integrate(gf::GridFunction)
    dom = domain(gf)
    cod = codomain(gf)
    VS = eltype(dom)
    DS = ndims(VS)
    S = eltype(VS)
    T = eltype(cod)
    imin = SVector{DS,Int}(first(axes(gf.grid)[d]) for d in 1:DS)
    imax = SVector{DS,Int}(last(axes(gf.grid)[d]) for d in 1:DS)
    h = prod((last(dom) - first(dom)) ./ (size(gf.grid) .- 1))::S
    s = zero(T)
    for i0 in CartesianIndex(Tuple(imin)):CartesianIndex(Tuple(imax))
        i = SVector{DS,Int}(Tuple(i0))
        w = prod(i[d] == imin[d] || i[d] == imax[d] ? one(S) / 2 : one(S) for d in 1:DS)
        s += T(w * gf.grid[i0])
    end
    return T(h * s)
end

end

module KernelFunctions

# using HCubature
using KernelAbstractions
# using SparseArrays
# using SparseArraysCOO
using StaticArrays

using ..Boxes
using ..Categories
using ..Defs

struct KernelFunction{DS,S,DT,T} <: Category{Box{DS,S},Box{DT,T}}
    name::AbstractString
    domain::Box{DS,S}
    codomain::Box{DT,T}
    grid::AbstractArray{SVector{DT,T},DS}
    function KernelFunction{DS,S,DT,T}(
        name::AbstractString,
        domain::Box{DS,S},
        codomain::Box{DT,T},
        grid::AbstractArray{SVector{DT,T},DS},
    ) where {DS,S,DT,T}
        all(size(grid) .>= 2) || throw(
            ArgumentError("KernelFunction size must be at least 2 in each dimension so that linear functions can be represented"),
        )
        return new{DS,S,DT,T}(name, domain, codomain, grid)
    end
end
function KernelFunction(
    name::AbstractString,
    domain::Box{DS,S},
    codomain::Box{DT,T},
    grid::AbstractArray{SVector{DT,T},DS},
) where {DS,S,DT,T}
    return KernelFunction{DS,S,DT,T}(name, domain, codomain, grid)
end
export KernelFunction

"""
    KernelFunction{DS,S,DT,T}(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_size::SVector{DS,Int})
    KernelFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_size::SVector{DS,Int})

Create a zero-valued Kernel function, i.e. a Kernel function that is zero
everywhere. The `codomain` can be a zero-valued domain created via
[`Box{DT,T}()`](@ref).

`grid_size` specifies the resolution of the Kernel function. This makes
it convenient to use zero-valued Kernel functions as "templates" or
"skeletons" when creating other Kernel functions.
"""
function KernelFunction{DS,S,DT,T}(
    make_zeros,
    name::AbstractString,
    domain::Box{DS,S},
    codomain::Box{DT,T},
    grid_size::SVector{DS,Int},
) where {DS,S,DT,T}
    # grid = Metal.zeros(SVector{DT,T}, Tuple(grid_size))
    grid = make_zeros(SVector{DT,T}, Tuple(grid_size))
    return KernelFunction(name, domain, Box{DT,T}(), grid::AbstractArray{SVector{DT,T},DS})
end
function KernelFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_size::SVector{DS,Int}) where {DS,S,DT,T}
    return KernelFunction{DS,S,DT,T}(name, domain, codomain, grid_size)
end

# Metadata
Categories.name(kf::KernelFunction) = kf.name
Categories.domain(kf::KernelFunction) = kf.domain
Categories.codomain(kf::KernelFunction) = kf.codomain

# Equality
function Base.:(==)(kf1::KernelFunction, kf2::KernelFunction)
    domain(kf1) == domain(kf2) || return false
    codomain(kf1) == codomain(kf2) || return false
    return kf1.grid == kf2.grid
end

# Collection
# @kernel function getindex_kernel!(result::AbstractArray{0,T}, @Const(grid::AbstractArray{D,T}), @Const(i)) where {D,T}
#     result[] = grid[i]
#     if false end
# end

Base.isempty(x::KernelFunction) = isempty(x.grid)
Base.length(x::KernelFunction) = length(x.grid)
function Base.getindex(x::KernelFunction, i)
    # backend = KernelAbstractions.get_backend(x.grid)
    # kernel! = getindex_kernel!(backend)
    # result = x.undefs(SVector{DT,T}, ())
    # kernel!(result, x.grid, i)
    # return Array(result)[]
    return Array(view(x.grid, i))[]
end
Base.getindex(x::KernelFunction{DS}, i::SVector{DS}) where {DS} = getindex(x, LinearIndices(x.grid)[i...])
function Base.map(f, x::KernelFunction)
    grid′ = map(f, x.grid)
    VT = eltype(grid′)
    @assert VT <: SVector
    cod = Box(extrema(grid′)...)
    return KernelFunction(x.name, x.domain, cod, grid′)
end
function Base.map(f, x::KernelFunction{DS}, y::KernelFunction{DS}) where {DS}
    @assert domain(x) == domain(y)
    @assert axes(x.grid) == axes(y.grid)
    grid′ = map(f, x.grid, y.grid)
    VT = eltype(grid′)
    @assert VT <: SVector
    cod = Box(extrema(grid′)...)
    return KernelFunction1D(x.name, x.domain, cod, grid′)
end

# Category
@kernel function make_identity_kernel!(grid::AbstractArray{S,DS},
        imin::SVector{DS,Int}, imax::SVector{DS,Int}, 
        xmin::SVector{DS,S}, xmax::SVector{DS,S}) where {DS,S}
    i0 = @index(Global, NTuple)
    i = SVector{DS,Int}(i0)
    x = lincom(imin, xmin, imax, xmax, i)
    grid[i0...] = x
    if false
    end
end

function Categories.make_identity(kf::KernelFunction)
    dom = kf.domain
    VS = eltype(dom)::Type{<:SVector}
    DS = ndims(dom)
    xmin = first(dom)
    xmax = last(dom)
    imin = SVector{DS,Int}(first.(axes(kf.grid)))
    imax = SVector{DS,Int}(last.(axes(kf.grid)))
    backend = KernelAbstractions.get_backend(kf.grid)
    kernel! = make_identity_kernel!(backend)
    grid = similar(kf.grid, VS)
    ndrange = ndims(grid) == 0 ? (1,) : size(grid)
    kernel!(grid, imin, imax, xmin, xmax; ndrange=ndrange)
    return KernelFunction(kf.name, dom, dom, grid)
end
function Base.:∘(kf2::KernelFunction{<:Any,<:Any,DT}, kf1::KernelFunction{DT}) where {DT}
    @assert domain(kf2) == codomain(kf1)
    return map(kf2, kf1)
end

# Vector space
Base.zero(x::KernelFunction) = map(zero, x)
Base.:+(x::KernelFunction) = map(+, x)
Base.:-(x::KernelFunction) = map(-, x)
Base.:+(x::KernelFunction, y::KernelFunction) = map(+, x, y)
Base.:-(x::KernelFunction, y::KernelFunction) = map(-, x, y)
Base.:*(a, x::KernelFunction) = map(w -> a * w, x)
Base.:*(x::KernelFunction, a) = map(w -> w * a, x)
Base.:\(a, x::KernelFunction) = map(w -> a \ w, x)
Base.:/(x::KernelFunction, a) = map(w -> w / a, x)

function Categories.project(kf::KernelFunction, cat::Category)
    @assert domain(kf) == domain(cat)
    dom = domain(cat)
    cod = codomain(cat)
    VS = eltype(dom)
    DS = ndims(dom)
    S = eltype(VS)
    VT = eltype(cod)
    DT = ndims(cod)
    T = eltype(VT)

    imin = SVector{DS,Int}(first.(axes(kf.grid)))
    imax = SVector{DS,Int}(last.(axes(kf.grid)))
    dx = (last(dom) - first(dom)) / (imax - imin)

    # Create mass matrix (TODO: Cache this / implement this as structural matrix)
    Mshape = size(kf.grid)
    Msize = prod(Mshape)
    @assert all(first.(axes(kf.grid)) .== 1)
    function Mindex(i)
        off = 0
        str = 1
        for d = 1:DS
            off += str * (i[d] - 1)
            str *= Mshape[d]
        end
        return off + 1
    end
    M = SparseMatrixCOO{S}(Msize, Msize)
    for i0 in CartesianIndices(axes(kf.grid))
        i = SVector{DS,Int}(Tuple(i0))
        for di0 in CartesianIndices(ntuple(d -> -1:+1, DS))
            di = SVector{DS,Int}(Tuple(di0))
            j = i + di
            if all(j .∈ axes(kf.grid))
                Mvalue =
                    prod(ifelse(di[d] == 0, one(S) * 2 / 3, one(S) / 6) for d = 1:DS; init = one(S)) /
                    2^count(i .== j .== imin .|| i .== j .== imax)
                M[Mindex(i), Mindex(j)] = Mvalue
            end
        end
    end
    M = sparse(M)

    Kernel = VT[
        let
            i = SVector{DS,Int}(Tuple(i0))
            s = zero(VT)
            for di0 in CartesianIndices(ntuple(d -> -1:0, DS))
                di = SVector{DS,Int}(Tuple(di0))
                j = i + di
                if all(imin .<= j .&& j .+ 1 .<= imax)
                    x0 = lincom(imin, first(dom), imax, last(dom), j .+ 0)
                    x1 = lincom(imin, first(dom), imax, last(dom), j .+ 1)

                    # TODO: normalize x to -1...+1 to improve floating
                    # point accuracy, and so that we don't thave to
                    # put the Kernel spacing into `M`

                    y0 = SVector{DS,S}(1 + di[d] for d = 1:DS)
                    y1 = SVector{DS,S}(0 - di[d] for d = 1:DS)
                    basis(x) = prod(lincom(x0[d], y0[d], x1[d], y1[d], x[d]) for d = 1:DS; init = one(S))::S

                    s += hcubature(x -> VT(basis(x) * cat(x)), x0, x1)[1]::VT
                end
            end
            s / prod(dx)
        end for i0 in CartesianIndices(axes(kf.grid))
    ]
    Kernel::Array{VT,DS}

    # Kernel = reshape(M \ reshape(Kernel, :), size(Kernel))
    @assert DT == 1
    Kernel = reinterpret(T, Kernel)
    Kernel = reshape(M \ reshape(Kernel, :), size(Kernel))
    Kernel = reinterpret(VT, Kernel)

    return KernelFunction{DS,S,DT,T}(cat.name, dom, cod, Kernel)
end

function Categories.evaluate(kf::KernelFunction{DS,S,DT,T}, x::SVector{DS,S}) where {DS,S,DT,T}
    isempty(kf.grid) && return zero(SVector{DT,T})::SVector{DT,T}
    length(kf.grid) == 1 && return kf.grid[begin]::SVector{DT,T}

    # x = first <=> i = firstindex
    # x = last  <=> i = lastindex
    dom = domain(kf)
    VS = eltype(dom)
    imin = SVector{DS,Int}(first.(axes(kf.grid)))
    imax = SVector{DS,Int}(last.(axes(kf.grid)))
    ix = lincom(first(dom), VS(imin), last(dom), VS(imax), x)
    i = SVector{DS,Int}(clamp(floor(Int, ix[d]), imin[d]:(imax[d]-1)) for d = 1:DS)
    q = (ix - i)::VS

    fx = zero(SVector{DT,T})
    for di0 = CartesianIndex(ntuple(d -> 0, DS)):CartesianIndex(ntuple(d -> 1, DS))
        di = SVector{DS,Int}(Tuple(di0))
        w = one(S)
        for d = 1:DS
            w *= di[d] == 0 ? 1 - q[d] : q[d]
        end
        fx += SVector{DT,T}(w * kf.grid[CartesianIndex(Tuple(i + di))])
    end
    return fx::SVector{DT,T}
end
Categories.evaluate(kf::KernelFunction{DS,S}, x::SVector{DS}) where {DS,S} = evaluate(kf, SVector{DS,S}(x))

function Categories.integrate(kf::KernelFunction)
    dom = domain(kf)
    cod = codomain(kf)
    VS = eltype(dom)
    DS = ndims(VS)
    S = eltype(VS)
    T = eltype(cod)
    imin = SVector{DS,Int}(first(axes(kf.grid)[d]) for d = 1:DS)
    imax = SVector{DS,Int}(last(axes(kf.grid)[d]) for d = 1:DS)
    h = prod((last(dom) - first(dom)) ./ (size(kf.grid) .- 1))::S
    s = zero(T)
    for i0 = CartesianIndex(Tuple(imin)):CartesianIndex(Tuple(imax))
        i = SVector{DS,Int}(Tuple(i0))
        w = prod(i[d] == imin[d] || i[d] == imax[d] ? one(S) / 2 : one(S) for d = 1:DS)
        s += T(w * kf.grid[i0])
    end
    return T(h * s)
end

end

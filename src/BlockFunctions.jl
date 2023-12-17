module BlockFunctions

using StaticArrays

using ..Boxes
using ..Categories
using ..Defs

struct BlockFunction{DS,S,DT,T} <: Category{Box{DS,S},Box{DT,T}}
    name::AbstractString
    domain::Box{DS,S}
    codomain::Box{DT,T}
    block_size::SVector{DS,Int}
    blocks::AbstractArray{Category{Box{DS,S},Box{DT,T}},DS}
    function BlockFunction{DS,S,DT,T}(
        name::AbstractString,
        domain::Box{DS,S},
        codomain::Box{DT,T},
        block_size::SVector{DS,Int},
        blocks::AbstractArray{<:Category{Box{DS,S},Box{DT,T}},DS},
    ) where {DS,S,DT,T}
        isempty(blocks) && throw(ArgumentError("BlockFunction: blocks cannot be empty so that linear functions can be represented"))
        # This assumes that the blocks are implemented as some kind of `AbstractArray{<:Any,DS}`.
        # Should we instroduce `AbstractGridFunction`?
        all(SVector{DS,Int}(size(b.grid)) == block_size for b in blocks) ||
            throw(ArgumentError("BlockFunction: all blocks must have size `block_size`"))
        # TODO: Check domains and codomains of all blocks
        return new{DS,S,DT,T}(name, domain, codomain, block_size, blocks)
    end
end
function BlockFunction(
    name::AbstractString,
    domain::Box{DS,S},
    codomain::Box{DT,T},
    # Should we use `NTuple{DS,OrdinalRange}` for this?
    block_size::SVector{DS,Int},
    # Should we use `Category` instead of `AbstractArray` here?
    blocks::AbstractArray{<:Category{Box{DS,S},Box{DT,T}},DS},
) where {DS,S,DT,T}
    return BlockFunction{DS,S,DT,T}(name, domain, codomain, block_size, blocks)
end
export BlockFunction

function block_domain(dom::Box{DS,S}, block_size::SVector{DS,Int}, i::SVector{DS,Int}) where {DS,S}
    imin = SVector{DS,Int}(1 for d in 1:DS)
    imax = block_size
    xmin = lincom(imin, first(dom), imax .+ 1, last(dom), i)
    xmax = lincom(imin, first(dom), imax .+ 1, last(dom), i .+ 1)
    dom′ = Box{DS,S}(xmin, xmax)
    return dom′
end

function BlockFunction{DS,S,DT,T}(
    make_block,
    name::AbstractString,
    domain::Box{DS,S},
    codomain::Box{DT,T},
    block_size::SVector{DS,Int},
    blocks_size::SVector{DS,Int},
) where {DS,S,DT,T}
    blocks = [
        let
            name′ = "$name $(Tuple(i))"
            dom′ = block_domain(domain, block_size, SVector(Tuple(i)))
            make_block(name′, dom′, codomain, blocks_size)
        end for i in CartesianIndices(Tuple(block_size))
    ]
    return BlockFunction(name, domain, Box{DT,T}(), block_size, blocks)
end
function BlockFunction(
    make_block,
    name::AbstractString,
    domain::Box{DS,S},
    codomain::Box{DT,T},
    block_size::SVector{DS,Int},
    blocks_size::SVector{DS,Int},
) where {DS,S,DT,T}
    return BlockFunction{DS,S,DT,T}(make_block, name, domain, codomain, block_size, blocks_size)
end
function BlockFunction{DS,S,DT,T}(
    name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, block_size::SVector{DS,Int}, blocks_size::SVector{DS,Int}
) where {DS,S,DT,T}
    make_block(name, dom, cod, npoints) = GridFunction{DS,S,DT,T}(name, dom, cod, npoints)
    return BlockFunction{DS,S,DT,T}(make_block, name, domain, codomain, block_size, blocks_size)
end
function BlockFunction(
    name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, block_size::SVector{DS,Int}, blocks_size::SVector{DS,Int}
) where {DS,S,DT,T}
    return BlockFunction{DS,S,DT,T}(name, domain, codomain, block_size, blocks_size)
end

# Metadata
Categories.name(bf::BlockFunction) = bf.name
Categories.domain(bf::BlockFunction) = bf.domain
Categories.codomain(bf::BlockFunction) = bf.codomain

# Equality
function Base.:(==)(bf1::BlockFunction, bf2::BlockFunction)
    domain(bf1) == domain(bf2) || return false
    codomain(bf1) == codomain(bf2) || return false
    bf1.block_size == bf2.block_size || return false
    return bf1.blocks == bf2.blocks
end

# Collection
Base.isempty(x::BlockFunction) = (isempty(x.blocks) || all(isempty.(x.blocks)))::Bool
Base.length(x::BlockFunction) = sum(length(b) for b in x.blocks)::Int
function Base.getindex(x::BlockFunction{DS,S,DT,T}, i) where {DS,S,DT,T}
    imin = SVector{DS,Int}(1 for d in 1:DS)
    imax = x.block_size .* size(x.blocks)
    ilen = imax - imin .+ 1
    istr = DS == 0 ? SVector{DS,Int}() : SVector{DS,Int}(one(S), view(cumprod(ilen), 1:(DS - 1))...)
    ivec = mod.(fld.(i - 1, istr), ilen) + imin
    ivec::SVector{DS,<:Integer}
    @assert (DS == 0 ? zero(S) : sum((ivec - imin) .* istr)) + 1 == i
    return x[ivec]::SVector{DT,T}
end
function Base.getindex(x::BlockFunction{DS,S,DT,T}, i::SVector{DS,<:Integer}) where {DS,S,DT,T}
    return x.blocks[CartesianIndex(Tuple(fld1.(i, x.block_size)))][mod1.(i, x.block_size)]::SVector{DT,T}
end
function Base.map(f, x::BlockFunction)
    blocks′ = map(b -> map(f, b), x.blocks)
    VT = eltype(eltype(blocks′))
    @assert VT <: SVector
    cod = Box(reduce(minmax, map(extrema, blocks′))...)
    return BlockFunction(x.name, x.domain, cod, blocks′)
end
function Base.map(f, x::BlockFunction{DS}, y::BlockFunction{DS}) where {DS}
    @assert domain(x) == domain(y)
    @assert axes(x.block) == axes(y.block)
    blocks′ = map((bx, by) -> map(f, bx, by), x.blocks, y.blocks)
    VT = eltype(eltype(blocks′))
    @assert VT <: SVector
    cod = Box(reduce(minmax, map(extrema, blocks′))...)
    return BlockFunction1D(x.name, x.domain, cod, blocks′)
end

# Category
function Categories.make_identity(bf::BlockFunction)
    blocks = map(make_identity, bf.blocks)
    return BlockFunction(bf.name, bf.domain, bf.domain, bf.block_size, blocks)
end
function Base.:∘(bf2::BlockFunction{<:Any,<:Any,DT}, bf1::BlockFunction{DT}) where {DT}
    @assert domain(bf2) == codomain(bf1)
    return map(bf2, bf1)
end

# Vector space
Base.zero(x::BlockFunction) = map(zero, x)
Base.:+(x::BlockFunction) = map(+, x)
Base.:-(x::BlockFunction) = map(-, x)
Base.:+(x::BlockFunction, y::BlockFunction) = map(+, x, y)
Base.:-(x::BlockFunction, y::BlockFunction) = map(-, x, y)
Base.:*(a, x::BlockFunction) = map(w -> a * w, x)
Base.:*(x::BlockFunction, a) = map(w -> w * a, x)
Base.:\(a, x::BlockFunction) = map(w -> a \ w, x)
Base.:/(x::BlockFunction, a) = map(w -> w / a, x)

function Categories.project(bf::BlockFunction, cat::Category)
    @assert domain(bf) == domain(cat)
    blocks = map(b -> project(b, cat), bf.blocks)
    return BlockFunction{DS,S,DT,T}(cat.name, dom, cod, blocks)
end

function Categories.evaluate(bf::BlockFunction{DS,S,DT,T}, x::SVector{DS,S}) where {DS,S,DT,T}
    isempty(bf.blocks) && return zero(SVector{DT,T})::SVector{DT,T}

    # x = first <=> i = firstindex
    # x = last  <=> i = lastindex
    dom = domain(bf)
    VS = eltype(dom)
    imin = SVector{DS,Int}(first.(axes(bf.blocks)))
    imax = SVector{DS,Int}(last.(axes(bf.blocks)))
    ix = lincom(first(dom), VS(imin), last(dom), VS(imax .+ 1), x)
    i = SVector{DS,Int}(clamp(floor(Int, ix[d]), imin[d]:imax[d]) for d in 1:DS)

    return evaluate(bf.blocks[CartesianIndex(Tuple(i))], x)::SVector{DT,T}
end
Categories.evaluate(bf::BlockFunction{DS,S}, x::SVector{DS}) where {DS,S} = evaluate(bf, SVector{DS,S}(x))

Categories.integrate(bf::BlockFunction) = sum(integrate, bf.blocks)

end

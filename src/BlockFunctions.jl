module BlockFunctions

using StaticArrays

using ..Boxes
using ..Categories
using ..Defs

struct BlockFunction{DS,S,DT,T} <: Category{Box{DS,S},Box{DT,T}}
    name::AbstractString
    domain::Box{DS,S}
    codomain::Box{DT,T}
    blocks::AbstractArray{Category{Box{DS,S},Box{DT,T}},DS}
    block_lengths_cumsum::Array{Int,DS}
    function BlockFunction{DS,S,DT,T}(
        name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, blocks::AbstractArray{<:Category{Box{DS,S},Box{DT,T}},DS}
    ) where {DS,S,DT,T}
        isempty(blocks) && throw(ArgumentError("BlockFunction: blocks cannot be empty so that linear functions can be represented"))
        block_lengths_cumsum = accumulate(+, length.(blocks))
        # TODO: Check domains and codomains of all blocks
        return new{DS,S,DT,T}(name, domain, codomain, blocks, block_lengths_cumsum)
    end
end
function BlockFunction(
    name::AbstractString,
    domain::Box{DS,S},
    codomain::Box{DT,T},
    # Should we use `Category` instead of `AbstractArray` here?
    blocks::AbstractArray{<:Category{Box{DS,S},Box{DT,T}},DS},
) where {DS,S,DT,T}
    return BlockFunction{DS,S,DT,T}(name, domain, codomain, blocks)
end
export BlockFunction

function block_domain(dom::Box{DS,S}, num_blocks::SVector{DS,Int}, i::SVector{DS,Int}) where {DS,S}
    imin = SVector{DS,Int}(1 for d in 1:DS)
    imax = num_blocks
    xmin = lincom(imin, first(dom), imax .+ 1, last(dom), i)
    xmax = lincom(imin, first(dom), imax .+ 1, last(dom), i .+ 1)
    dom′ = Box{DS,S}(xmin, xmax)
    return dom′::Box{DS,S}
end

function BlockFunction{DS,S,DT,T}(
    make_block, name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, num_blocks::SVector{DS,Int}
) where {DS,S,DT,T}
    blocks = [
        let
            name′ = "$name $(Tuple(i))"
            dom′ = block_domain(domain, num_blocks, SVector{DS,Int}(Tuple(i)))
            make_block(name′, dom′, codomain)::Category{Box{DS,S},Box{DT,T}}
        end for i in CartesianIndices(Tuple(num_blocks))
    ]
    return BlockFunction(name, domain, Box{DT,T}(), blocks)
end
function BlockFunction(
    make_block, name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, num_blocks::SVector{DS,Int}
) where {DS,S,DT,T}
    return BlockFunction{DS,S,DT,T}(make_block, name, domain, codomain, num_blocks)
end

# Metadata
Categories.name(bf::BlockFunction) = bf.name
Categories.domain(bf::BlockFunction) = bf.domain
Categories.codomain(bf::BlockFunction) = bf.codomain

# Equality
function Base.:(==)(bf1::BlockFunction, bf2::BlockFunction)
    domain(bf1) == domain(bf2) || return false
    codomain(bf1) == codomain(bf2) || return false
    return bf1.blocks == bf2.blocks
end

# Collection
Base.isempty(x::BlockFunction) = all(isempty.(x.blocks))::Bool
Base.length(x::BlockFunction) = sum(length(b) for b in x.blocks)::Int
function Base.getindex(x::BlockFunction{DS,S,DT,T}, i::Int) where {DS,S,DT,T}
    @assert i >= 1
    bi = searchsortedfirst(reshape(x.block_lengths_cumsum, :), i)
    @assert bi <= length(x.blocks)
    i0 = bi == 1 ? 0 : x.block_lengths_cumsum[bi - 1]
    return getindex(x.blocks[bi], i - i0)::SVector{DT,T}
end
Base.getindex(x::BlockFunction, i::Integer) = getindex(x, Int(i))
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
    return BlockFunction(bf.name, bf.domain, bf.domain, blocks)
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
    bx = lincom(first(dom), VS(imin), last(dom), VS(imax .+ 1), x)
    bi = SVector{DS,Int}(clamp(floor(Int, bx[d]), imin[d]:imax[d]) for d in 1:DS)

    return evaluate(bf.blocks[CartesianIndex(Tuple(bi))], x)::SVector{DT,T}
end
Categories.evaluate(bf::BlockFunction{DS,S}, x::SVector{DS}) where {DS,S} = evaluate(bf, SVector{DS,S}(x))

Categories.integrate(bf::BlockFunction) = sum(integrate, bf.blocks)

end

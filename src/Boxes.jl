module Boxes

using StaticArrays

using ..Domains

"""
    struct Box{D,T} <: Domain{SVector{D,T}}

A multi-dimensional [`Domains.Domain`](@ref), i.e. a box or
hyper-rectangle.  The number of dimensions `D` must be a non-negative
integer. `T` should be a real-valued scalar type such as
`Float64`. The domain's element type is `SVector{D,T}`.

Zero-dimensional and one-dimensional domains are supported.
Zero-dimensional domains contain only a single point (the origin).
One-dimensional domains correspond to intervals (see
[`Interval`](#GraviPet.Intervals.Interval)).
"""
struct Box{D,T} <: Domain{SVector{D,T}}
    first::SVector{D,T}
    last::SVector{D,T}
    function Box{D,T}(first::SVector{D,T}, last::SVector{D,T}) where {D,T<:Real}
        all(first .<= last) || throw(ArgumentError("Intervals must be non-empty"))
        return new{D,T}(first, last)
    end
end

"""
    Box{D,T}(first::SVector{D}, last::SVector{D})
    Box(first::SVector{D}, last::SVector{D})

Create a box domain with the specified lower and upper bounds. The
bounds are given as vectors, and (some of) their elements can be
infinity. The constraint `first < last` must hold in each dimension.
"""
Box{D,T}(first::SVector{D,<:Real}, last::SVector{D,<:Real}) where {D,T} = Box{D,T}(SVector{D,T}(first), SVector{D,T}(last))
Box(first::SVector{D,T}, last::SVector{D,T}) where {D,T<:Real} = Box{D,T}(first, last)
Box(first::SVector{D,<:Real}, last::SVector{D,<:Real}) where {D} = Box(promote(first, last)...)
export Box

"""
    Box{D,T}()

Create a zero-valued box domain, i.e. a box where both lower and upper
bounds are zero. This is also a convenient way to provide a "skeleton"
for a domain without having to specify any bounds.
"""
Box{D,T}() where {D,T} = Box{D,T}(zero(SVector{D,T}), zero(SVector{D,T}))

# Metadata
Base.ndims(::Box{D}) where {D} = D
Base.first(box::Box) = box.first
Base.last(box::Box) = box.last

# Equality
Base.:(==)(box1::Box, box2::Box) = box1.first == box2.first && box1.last == box2.last

# Domain operations
function Domains.expanded(box::Box{D,T}, delta_lo::SVector{D,T}, delta_hi::SVector{D,T}) where {D,T}
    return Box{D,T}(box.first - delta_lo, box.last + delta_hi)
end

# Set relations
Base.issubset(box1::Box{0}, box2::Box{0}) = true
Base.issubset(box1::Box{D}, box2::Box{D}) where {D} = all((box1.first .≥ box2.first) .& (box1.last .≤ box2.last))
Base.isdisjoint(box1::Box{0}, box2::Box{0}) = false
Base.isdisjoint(box1::Box{D}, box2::Box{D}) where {D} = all((box1.last .< box2.first) .| (box2.last .< box1.first))

# Membership test
Base.in(x, box::Box) = all(first(box) .≤ x .≤ last(box))

end

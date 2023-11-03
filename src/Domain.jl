"""
    abstract type Domain{T}

A computational domain of type `T`. The *domain* of a function is the
type of its argument(s). For example, `Domain{Float64}` would describe
a one-dimensional domain, and `Domain{SVector{3,Float64}}` would
describe a three-dimensional domain.

This type can also describe a codomain, i.e. the result type of a
function. A codomain of `Domain{Float64}` describes a scalar function,
and a codomain of `Domain{SVector{3,Float64}}` describes a
vector-valued function.

See also [`Interval`](@ref) and [`Box`](@ref) for concrete domain
types.
"""
abstract type Domain{T} end
export Domain

# Metadata

"""
    eltype(dom::Domain)::Type

Return the element type of a [`Domain`](@ref).
"""
Base.eltype(::Domain{T}) where {T} = T

"""
    ndims(dom::Domain)::Int

Return the number of dimensions of a [`Domain`](@ref). Scalar domains
are one-dimensional. Zero-dimensional domains consist of a single
point only.
"""
Base.ndims(::Domain) = throw(MethodError(ndims, (dom,)))

"""
    first(dom::Domain{T})::T

Return the lower bound of the domain. (At the moment domains are
assumed to be cuboid, i.e. to form a box or hyper-rectangle.) Domain
bounds can be infinity.

See also [`last`](@ref).
"""
Base.first(dom::Domain) = throw(MethodError(first, (dom,)))

"""
    last(dom::Domain{T})::T

Return the upper bound of the domain. (At the moment domains are
assumed to be cuboid, i.e. to form a box or hyper-rectangle.) Domain
bounds can be infinity.

See also [`first`](@ref).
"""
Base.last(dom::Domain) = throw(MethodError(last, (dom,)))

# Equality

"""
    ==(dom1::Domain, dom2::Domain)::Bool

Equality for domains. Domains are equal if the have the same number of
dimensions (see [`ndims`](@ref)), and the same lower and upper bounds
(see [`first`](@ref) and [`last`](@ref)). The element types might be
different if they can be compared for equality, e.g. `Float64` and
`BigFloat`.

Incompatible domains (living in different number of dimensions) are
treated as unequal.
"""
Base.:(==)(dom1::Domain, dom2::Domain) = false # throw(MethodError(==, (dom1, dom2)))

# Set relations

"""
    issubset(dom1::Domain, dom2::Domain)::Bool

Test whether `dom1` is a subset of `dom2`. Only compatible domains
(which live in the same number of dimensions) can be compared.
"""
Base.issubset(dom1::Domain, dom2::Domain) = throw(MethodError(issubset, (dom1, dom2)))

# Membership test

"""
    in(elem, dom::Domain)::Bool

Test whether `elem` is inside domain `dom`.
"""
Base.in(elem, dom::Domain) = throw(MethodError(in, (elem, dom)))

################################################################################

"""
    struct Interval{T} <: Domain{T}

A one-dimensional [`Domain`](@ref), i.e. an interval. `T` should be a
real-valued scalar type such as `Float64`.
"""
struct Interval{T} <: Domain{T}
    first::T
    last::T
    function Interval{T}(first::T, last::T) where {T<:Real}
        first <= last || throw(ArgumentError("Intervals must be non-empty"))
        return new{T}(first, last)
    end
end

"""
    Interval{T}(first::Real, last::Real)
    Interval(first::Real, last::Real)

Create an interval with the specified bounds. The bounds can be
infinity. The constraint `first < last` must hold.
"""
Interval{T}(first::Real, last::Real) where {T} = Interval{T}(T(first), T(last))
Interval(first::T, last::T) where {T<:Real} = Interval{T}(first, last)
Interval(first::Real, last::Real) = Interval(promote(first, last)...)
export Interval

"""
    Interval{T}()

Create a zero-valued interval, i.e. an interval where both lower and
upper bounds are zero. This is also a convenient way to provide a
"skeleton" for a domain without having to specify any bounds.
"""
Interval{T}() where {T} = Interval{T}(zero(T), zero(T))

# Metadata
Base.ndims(::Interval) = 1
Base.first(iv::Interval) = iv.first
Base.last(iv::Interval) = iv.last

# Equality
Base.:(==)(iv1::Interval, iv2::Interval) = iv1.first == iv2.first && iv1.last == iv2.last

# Set relations
Base.issubset(iv1::Interval, iv2::Interval) = iv1.first >= iv2.first && iv1.last <= iv2.last
Base.isdisjoint(iv1::Interval, iv2::Interval) = iv1.last < iv2.first || iv2.last < iv1.first

# Membership test
Base.in(x, iv::Interval) = first(iv) <= x <= last(iv)

################################################################################

"""
    struct Box{D,T} <: Domain{SVector{D,T}}

A multi-dimensional [`Domain`](@ref), i.e. a box or hyper-rectangle.
The number of dimensions `D` must be a non-negative integer. `T`
should be a real-valued scalar type such as `Float64`. The domain's
element type is `SVector{D,T}`.

Zero-dimensional and one-dimensional domains are supported.
Zero-dimensional domains contain only a single point (the origin).
One-dimensional domains correspond to intervals (see
[`Interval`](@ref)).
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

# Set relations
Base.issubset(box1::Box, box2::Box) = all((box1.first .>= box2.first) .& (box1.last .<= box2.last))
Base.isdisjoint(box1::Box, box2::Box) = all((box1.last .< box2.first) .| (box2.last .< box1.first))

# Membership test
Base.in(x, box::Box) = all(first(box) .<= x .<= last(box))

module Domains

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

See also [`Interval`](#GraviPet.Intervals.Interval) and
[`Box`](#GraviPet.Boxes.Box) for concrete domain types.
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

# Domain operations

"""
    expanded(dom::Domain, delta)::Domain
    expanded(dom::Domain, delta_lo, delta_hi)::Domain

Expand the domain `dom`. The lower and upper domain bounds are moved
by either `delta`, or `delta_lo` and `delta_hi`, respectively.
Positive deltas enlarge the domain, negative values shrink it. Domains
must be non-empty.

`expanded(dom, delta)` is equivalent to `expanded(dom, delta, delta)`.

See also [`shifted`](@ref).
"""
expanded(dom::Domain{T}, delta_lo::T, delta_hi::T) where {T} = throw(MethodError(expanded, (dom, delta_lo, delta_hi)))
expanded(dom::Domain{T}, delta_lo, delta_hi) where {T} = expanded(dom, T(delta_lo), T(delta_hi))
expanded(dom::Domain{T}, delta::T) where {T} = expanded(dom, delta, delta)
expanded(dom::Domain{T}, delta) where {T} = expanded(dom, T(delta))
export expanded

"""
    shifted(dom::Domain, shift)::Domain

Shift the domain `dom` by moving both domain bounds in the same
direction. This does not change the size of the domain.

`shfited(dom, shift)` is equivalent to `expanded(dom, -shift,
+shift)`.

See also [`expanded`](@ref).
"""
shifted(dom::Domain, shift) = expanded(dom, -shift, +shift)
export shifted

# Set relations

"""
    issubset(dom1::Domain, dom2::Domain)::Bool

Test whether `dom1` is a subset of `dom2`. Only compatible domains
(which live in the same number of dimensions) can be compared.
"""
Base.issubset(dom1::Domain, dom2::Domain) = throw(MethodError(issubset, (dom1, dom2)))

"""
    isdisjoint(dom1::Domain, dom2::Domain)::Bool

Test whether `dom1` is disjoint from `dom2`. Only compatible domains
(which live in the same number of dimensions) can be compared.
"""
Base.isdisjoint(dom1::Domain, dom2::Domain) = throw(MethodError(isdisjoint, (dom1, dom2)))

# Membership test

"""
    in(elem, dom::Domain)::Bool

Test whether `elem` is inside domain `dom`.
"""
Base.in(elem, dom::Domain) = throw(MethodError(in, (elem, dom)))

end

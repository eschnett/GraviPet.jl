module Intervals

using ..Domains

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

# Domain operations
Domains.expanded(iv::Interval{T}, delta_lo::T, delta_hi::T) where {T} = Interval{T}(iv.first - delta_lo, iv.last + delta_hi)

# Set relations
Base.issubset(iv1::Interval, iv2::Interval) = iv1.first >= iv2.first && iv1.last <= iv2.last
Base.isdisjoint(iv1::Interval, iv2::Interval) = iv1.last < iv2.first || iv2.last < iv1.first

# Membership test
Base.in(x, iv::Interval) = first(iv) ≤ x ≤ last(iv)

end

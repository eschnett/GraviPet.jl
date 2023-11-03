abstract type Domain{T} end
export Domain

# Metadata
Base.eltype(::Domain{T}) where {T} = T
Base.ndims(::Domain) = throw(MethodError(ndims, (dom,)))
Base.first(dom::Domain) = throw(MethodError(first, (dom,)))
Base.last(dom::Domain) = throw(MethodError(last, (dom,)))

# Equality
Base.:(==)(dom1::Domain, dom2::Domain) = false # throw(MethodError(==, (dom1, dom2)))

# Set relations
Base.issubset(dom1::Domain, dom2::Domain) = throw(MethodError(issubset, (dom1, dom2)))

# Membership test
Base.in(x, dom::Domain) = throw(MethodError(in, (x, dom)))

################################################################################

struct Interval{T} <: Domain{T}
    first::T
    last::T
    function Interval{T}(first::T, last::T) where {T<:Real}
        first <= last || throw(ArgumentError("Intervals must be non-empty"))
        return new{T}(first, last)
    end
end
Interval{T}(first::Real, last::Real) where {T} = Interval{T}(T(first), T(last))
Interval(first::T, last::T) where {T<:Real} = Interval{T}(first, last)
Interval(first::Real, last::Real) = Interval(promote(first, last)...)
export Interval

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

struct Box{D,T} <: Domain{SVector{D,T}}
    first::SVector{D,T}
    last::SVector{D,T}
    function Box{D,T}(first::SVector{D,T}, last::SVector{D,T}) where {D,T<:Real}
        all(first .<= last) || throw(ArgumentError("Intervals must be non-empty"))
        return new{D,T}(first, last)
    end
end
Box{D,T}(first::SVector{D,<:Real}, last::SVector{D,<:Real}) where {D,T} = Box{D,T}(SVector{D,T}(first), SVector{D,T}(last))
Box(first::SVector{D,T}, last::SVector{D,T}) where {D,T<:Real} = Box{D,T}(first, last)
Box(first::SVector{D,<:Real}, last::SVector{D,<:Real}) where {D} = Box(promote(first, last)...)
export Box

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

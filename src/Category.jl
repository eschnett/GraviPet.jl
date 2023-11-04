abstract type Category{Dom,Cod} end
export Category, name, domain, codomain, project, evaluate

# Metadata
name(cat::Category) = throw(MethodError(name, (cat,)))
domain(cat::Category) = throw(MethodError(domain, (cat,)))
codomain(cat::Category) = throw(MethodError(codomain, (cat,)))

# Equality
Base.:(==)(cat1::Category{Dom,Cod}, cat2::Category{Dom,Cod}) where {Dom,Cod} = false # throw(MethodError(==, (cat1, cat2)))

# Collection
Base.isempty(x::Category) = throw(MethodError(isempty, (x,)))
Base.length(x::Category) = throw(MethodError(length, (x,)))
Base.getindex(x::Category, i) = throw(MethodError(getindex, (x, i)))
Base.map(f, x::Category) = throw(MethodError(map, (f, x)))
Base.map(f, x::Category, y::Category) = throw(MethodError(map, (f, x, y)))

# Category
Base.identity(cat::Category) = throw(MethodError(identity, (cat,)))
Base.:∘(cat2::Category, cat1::Category) = throw(MethodError(∘, (cat2, cat1)))

# Vector space
Base.eltype(x::Category) = eltype(codomain(x))
Base.zero(x::Category) = throw(MethodError(zero, (x,)))
Base.:+(x::Category) = throw(MethodError(+, (x,)))
Base.:-(x::Category) = throw(MethodError(-, (y,)))
Base.:+(x::Category, y::Category) = throw(MethodError(-, (x, y)))
Base.:-(x::Category, y::Category) = throw(MethodError(-, (x, y)))
Base.:*(a, x::Category) = throw(MethodError(*, (a, x)))
Base.:*(x::Category, a) = throw(MethodError(*, (x, a)))
Base.:\(a, x::Category) = throw(MethodError(\, (a, x)))
Base.:/(x::Category, a) = throw(MethodError(/, (x, a)))

"""
    project(dst::T, src::Category)::T where {T<:Category}

Project the function `src` onto the representation given by `dst`.
This can be used e.g. to sample Julia functions on a given grid, or to
resample functiont on a different grid. `dst` can be a skeleton
representation.
"""
project(dst::Category, src::Category) = throw(MethodError(project, (dst, src)))

"""
    evaluate(cat::Category, x)::eltype(codomain(cat))

Evaluate the function `cat` at point `x`. This can also be written as
a function call `cat(x)`.
"""
evaluate(cat::Category, x) = throw(MethodError(evaluate, (cat, x)))

"""
    (cat::Category)(x)::eltype(codomain(cat))

Evaluate the function `cat` at point `x`. This is written as a
function call `cat(x)`. See also [`evaluate`](@ref).
"""
(cat::Category)(x) = evaluate(cat, x)

################################################################################

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
name(gf::GridFunction1D) = gf.name
domain(gf::GridFunction1D) = gf.domain
codomain(gf::GridFunction1D) = gf.codomain

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
function Base.identity(gf::GridFunction1D)
    dom = gf.domain
    T = eltype(dom)
    @assert length(gf.grid) >= 2
    grid = if isempty(gf.grid)
        T[]
    elseif length(gf.grid) == 1
        T[first(dom)]
    else
        range(first(dom), last(dom); length=length(gf.grid))
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
function project(gf::GridFunction1D, cat::Category)
    @assert domain(gf) == domain(cat)
    grid = map(cat, gf.grid)
    return GridFunction1D(name(cat), domain(cat), codomain(cat), grid)
end

function evaluate(gf::GridFunction1D{<:Any,T}, x::T) where {T}
    isempty(gf.grid) && return zero(T)::T
    length(gf.grid) == 1 && return gf.grid[begin]::T

    # x = first <=> i = firstindex
    # x = last  <=> i = lastindex
    dom = domain(gf)
    ix =
        firstindex(gf.grid) * ((x - last(dom)) / (first(dom) - last(dom))) +
        lastindex(gf.grid) * ((x - first(dom)) / (last(dom) - first(dom)))
    i = clamp(floor(Int, ix), firstindex(gf.grid):(lastindex(gf.grid) - 1))
    di = ix - i

    fx = (1 - di) * gf.grid[i] + di * gf.grid[i + 1]
    return fx::T
end
evaluate(gf::GridFunction1D{<:Any,T}, x) where {T} = evaluate(gf, T(x))

################################################################################

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
    GridFunction{DS,S,DT,T}(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_shape::SVector{DS,Int})
    GridFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_shape::SVector{DS,Int})

Create a zero-valued grid function, i.e. a grid function that is zero
everywhere. The `codomain` can be a zero-valued domain created via
[`Box{DT,T}()`](@ref).

`grid_shape` specifies the resolution of the grid function. This makes
it convenient to use zero-valued grid functions as "templates" or
"skeletons" when creating other grid functions.

This zero-valued grid function is represented efficiently and does
*not* store one element per grid point.
"""
function GridFunction{DS,S,DT,T}(
    name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_shape::SVector{DS,Int}
) where {DS,S,DT,T}
    return GridFunction(name, domain, Box{DT,T}(), Zeros{T}(Tuple(grid_shape)))
end
function GridFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_shape::SVector{DS,Int}) where {DS,S,DT,T}
    return GridFunction{DS,S,DT,T}(name, domain, codomain, grid_shape)
end

# Metadata
name(gf::GridFunction) = gf.name
domain(gf::GridFunction) = gf.domain
codomain(gf::GridFunction) = gf.codomain

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
    grid′ = map(f, gf.grid)
    VT = eltype(grid′)
    VT::SVector
    cod = Box(extrema(grid′)...)
    return GridFunction(x.name, x.domain, cod, grid′)
end
function Base.map(f, x::GridFunction{DS}, y::GridFunction{DS}) where {DS}
    @assert domain(x) == domain(y)
    @assert axes(x.grid) == axes(y.grid)
    grid′ = map(f, x.grid, y.grid)
    VT = eltype(grid′)
    VT::SVector
    cod = Box(extrema(grid′)...)
    return GridFunction1D(x.name, x.domain, cod, grid′)
end

# Category
function Base.identity(gf::GridFunction)
    dom = gf.domain
    VT = eltype(dom)::Type{<:SVector}
    DT = length(VT)::Int
    xmin = first(dom)
    xmax = last(dom)
    imin = SVector{DT,Int}(first.(axes(gf.grid)))
    imax = SVector{DT,Int}(last.(axes(gf.grid)))
    grid = VT[
        xmin .* (SVector{DT,Int}(Tuple(i)) - imax) ./ (imin - imax) + xmax .* (SVector{DT,Int}(Tuple(i)) - imin) ./ (imax - imin)
        for i in CartesianIndex(Tuple(imin)):CartesianIndex(Tuple(imax))
    ]
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

# We sample instead of projecting...
function project(gf::GridFunction, cat::Category)
    @assert domain(gf) == domain(cat)
    grid = map(cat, gf.grid)
    return GridFunction(name(cat), domain(cat), codomain(cat), grid)
end

function evaluate(gf::GridFunction{DS,S,DT,T}, x::SVector{DS,S}) where {DS,S,DT,T}
    isempty(gf.grid) && return zero(SVector{DT,T})::SVector{DT,T}
    length(gf.grid) == 1 && return gf.grid[begin]::SVector{DT,T}

    # x = first <=> i = firstindex
    # x = last  <=> i = lastindex
    dom = domain(gf)
    ix =
        SVector(first.(axes(gf.grid))) .* ((x - last(dom)) ./ (first(dom) - last(dom))) +
        SVector(last.(axes(gf.grid))) .* ((x - first(dom)) ./ (last(dom) - first(dom)))
    i = SVector{DS,Int}(clamp(floor(Int, ix[d]), first.(axes(gf.grid))[d]:(last.(axes(gf.grid))[d] - 1)) for d in 1:DS)
    q = (ix - i)::SVector{DS,S}

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
evaluate(gf::GridFunction{DS,S}, x::SVector{DS}) where {DS,S} = evaluate(gf, SVector{DS,S}(x))

################################################################################

struct JuliaFunction{DS,S,DT,T} <: Category{Box{DS,S},Box{DT,T}}
    name::AbstractString
    domain::Box{DS,S}
    codomain::Box{DT,T}
    fun::Any
end
export JuliaFunction

"""
    JuliaFunction{DS,S,DT,T}(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T})
    JuliaFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T})

Create a zero-valued Julia function, i.e. a Julia function that is
zero everywhere. The `codomain` can be a zero-valued domain created
via [`Box{DT,T}()`](@ref).
"""
function JuliaFunction{DS,S,DT,T}(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}) where {DS,S,DT,T}
    return JuliaFunction(name, domain, Box{DT,T}(), Returns(zero(SVector{DT,T})))
end
function JuliaFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}) where {DS,S,DT,T}
    return JuliaFunction{DS,S,DT,T}(name, domain, codomain)
end

function project(jf::JuliaFunction, cat::Category)
    @assert domain(jf) == domain(cat)
    return JuliaFunction(name(cat), domain(cat), codomain(cat), cat)
end

# Metadata
name(jf::JuliaFunction) = jf.name
domain(jf::JuliaFunction) = jf.domain
codomain(jf::JuliaFunction) = jf.codomain

# Equality
function Base.:(==)(jf1::JuliaFunction, jf2::JuliaFunction)
    domain(jf1) == domain(jf2) || return false
    codomain(jf1) == codomain(jf2) || return false
    return jf1.grid == jf2.grid
end

# Collection is not available
Base.isempty(x::JuliaFunction) = throw(MethodError(isempty, (x,)))
Base.length(x::JuliaFunction) = throw(MethodError(length, (x,)))
Base.getindex(x::JuliaFunction, i) = throw(MethodError(getindex, (x, i)))
Base.map(f, x::JuliaFunction) = throw(MethodError(map, (f, x)))
Base.map(f, x::JuliaFunction, y::JuliaFunction) = throw(MethodError(map, (f, x, y)))

# Category
Base.identity(jf::JuliaFunction) = JuliaFunction(jf.name, jf.dom, jf.dom, identity)
function Base.:∘(jf2::JuliaFunction{<:Any,<:Any,DT}, jf1::JuliaFunction{DT}) where {DT}
    @assert domain(jf2) == codomain(jf1)
    return JuliaFunction(name(jf1), domain(jf1), codomain(jf2), jf2.fun ∘ jf1.fun)
end

# Vector space
function Base.zero(x::JuliaFunction)
    cod = codomain(x)
    cod′ = Box(zero(first(cod)), zero(last(cod)))
    return JuliaFunction(name(x), domain(x), cod′, Returns(zero(eltype(cod′))))
end
function Base.:+(x::JuliaFunction)
    cod = codomain(x)
    cod′ = Box(first(cod), last(cod))
    return JuliaFunction(name(x), domain(x), cod′, w -> +x.fun(w))
end
function Base.:-(x::JuliaFunction)
    cod = codomain(x)
    cod′ = Box(-last(cod), -first(cod))
    return JuliaFunction(name(x), domain(x), cod′, w -> -x.fun(w))
end
function Base.:+(x::JuliaFunction, y::JuliaFunction)
    @assert domain(x) == domain(y)
    @assert ndims(codomain(x)) == ndims(codomain(y))
    codx = codomain(x)
    cody = codomain(y)
    cod′ = Box(min.(first(codx), first(cody)), max.(last(codx), last(cody)))
    return JuliaFunction(name(x), domain(x), cod′, w -> x.fun(w) + y.fun(w))
end
function Base.:-(x::JuliaFunction, y::JuliaFunction)
    @assert domain(x) == domain(y)
    @assert ndims(codomain(x)) == ndims(codomain(y))
    codx = codomain(x)
    cody = codomain(y)
    cod′ = Box(min.(first(codx), -last(cody)), max.(last(codx), -first(cody)))
    return JuliaFunction(name(x), domain(x), cod′, w -> x.fun(w) - y.fun(w))
end
function Base.:*(a, x::JuliaFunction)
    cod = codomain(x)
    if a < 0
        cod′ = Box(a * last(cody), a * first(codx))
    elseif a > 0
        cod′ = Box(a * first(codx), a * last(cody))
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    return JuliaFunction(name(x), domain(x), cod′, w -> a * x.fun(w))
end
function Base.:*(x::JuliaFunction, a)
    cod = codomain(x)
    if a < 0
        cod′ = Box(last(cody) * a, first(codx) * a)
    elseif a > 0
        cod′ = Box(first(codx) * a, last(cody) * a)
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    return JuliaFunction(name(x), domain(x), cod′, w -> x.fun(w) * a)
end
function Base.:\(a, x::JuliaFunction)
    cod = codomain(x)
    if a < 0
        cod′ = Box(a \ last(cody), a \ first(codx))
    elseif a > 0
        cod′ = Box(a \ first(codx), a \ last(cody))
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    return JuliaFunction(name(x), domain(x), cod′, w -> a \ x.fun(w))
end
function Base.:/(x::JuliaFunction, a)
    cod = codomain(x)
    if a < 0
        cod′ = Box(last(cody) / a, first(codx) / a)
    elseif a > 0
        cod′ = Box(first(codx) / a, last(cody) / a)
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    return JuliaFunction(name(x), domain(x), cod′, w -> x.fun(w) / a)
end

evaluate(jf::JuliaFunction{DS,S,DT,T}, x::SVector{DS,S}) where {DS,S,DT,T} = jf.fun(x)::SVector{DT,T}
evaluate(jf::JuliaFunction{DS,S,DT,T}, x::SVector{DS}) where {DS,S,DT,T} = evaluate(jf, SVector{DS,S}(x))
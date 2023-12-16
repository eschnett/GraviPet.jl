abstract type Category{Dom,Cod} end
export Category, name, domain, codomain, make_identity, project, evaluate

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
make_identity(cat::Category) = throw(MethodError(make_identity, (cat,)))
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
export project

"""
    evaluate(cat::Category, x)::eltype(codomain(cat))

Evaluate the function `cat` at point `x`. This can also be written as
a function call `cat(x)`.
"""
evaluate(cat::Category, x) = throw(MethodError(evaluate, (cat, x)))
export evaluate

"""
    (cat::Category)(x)::eltype(codomain(cat))

Evaluate the function `cat` at point `x`. This is written as a
function call `cat(x)`. See also [`evaluate`](@ref).
"""
(cat::Category)(x) = evaluate(cat, x)

"""
    integrate(cat::Category)::eltype(codomain(cat))

Integrate the function `cat` over its domain.
"""
integrate(cat::Category) = throw(MethodError(integrate, (cat,)))
export integrate

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
function make_identity(gf::GridFunction1D)
    dom = gf.domain
    T = eltype(dom)
    @assert length(gf.grid) >= 2
    grid = if isempty(gf.grid)
        T[]
    elseif length(gf.grid) == 1
        T[first(dom)]
    else
        range(T(first(dom)), T(last(dom)); length=length(gf.grid))
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
function project(gf::GridFunction1D{S}, cat::Category) where {S}
    @assert domain(gf) == domain(cat)
    id = make_identity(gf)
    return map(cat, id)::GridFunction1D{S}
end

function evaluate(gf::GridFunction1D{<:Any,T}, x::T) where {T}
    isempty(gf.grid) && return zero(T)::T
    length(gf.grid) == 1 && return gf.grid[begin]::T

    dom = domain(gf)
    S = eltype(dom)
    ix = lincom(first(dom), S(firstindex(gf.grid)), last(dom), S(lastindex(gf.grid)), x)
    i = clamp(floor(Int, ix), firstindex(gf.grid):(lastindex(gf.grid) - 1))
    di = ix - i

    fx = lincom(S(0), gf.grid[i], S(1), gf.grid[i + 1], di)
    return fx::T
end
evaluate(gf::GridFunction1D{<:Any,T}, x) where {T} = evaluate(gf, T(x))

function integrate(gf::GridFunction1D)
    dom = domain(gf)
    cod = codomain(gf)
    S = eltype(dom)
    T = eltype(cod)
    imin = first(axes(gf.grid)[1])
    imax = last(axes(gf.grid)[1])
    h = (last(dom) - first(dom)) / (length(gf.grid) - 1)
    s = zero(T)
    for i in imin:imax
        w = i == imin || i == imax ? one(S) / 2 : one(S)
        s += T(w * gf.grid[i])
    end
    return T(h * s)
end

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
function make_identity(gf::GridFunction)
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

function project(gf::GridFunction, cat::Category)
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

function evaluate(gf::GridFunction{DS,S,DT,T}, x::SVector{DS,S}) where {DS,S,DT,T}
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
evaluate(gf::GridFunction{DS,S}, x::SVector{DS}) where {DS,S} = evaluate(gf, SVector{DS,S}(x))

function integrate(gf::GridFunction)
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

################################################################################

# TODO: Add an implementation type `X` at the end and declare `fun::X`
# to improve efficiency?
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
make_identity(jf::JuliaFunction) = JuliaFunction(jf.name, jf.dom, jf.dom, identity)
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

function integrate(jf::JuliaFunction)
    dom = domain(jf)
    cod = codomain(jf)
    I, E = hcubature(jf.fun, first(dom), last(dom))
    return I::eltype(cod)
end

################################################################################

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
    xmin = GraviPet.lincom(imin, first(dom), imax .+ 1, last(dom), i)
    xmax = GraviPet.lincom(imin, first(dom), imax .+ 1, last(dom), i .+ 1)
    dom′ = Box{DS,S}(xmin, xmax)
    return dom′
end

function BlockFunction{DS,S,DT,T}(
    name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, block_size::SVector{DS,Int}, blocks_size::SVector{DS,Int}
) where {DS,S,DT,T}
    blocks = [
        let
            name′ = "$name $(Tuple(i))"
            dom′ = block_domain(domain, block_size, SVector(Tuple(i)))
            GridFunction{DS,S,DT,T}(name′, dom′, codomain, block_size)
        end for i in CartesianIndices(Tuple(block_size))
    ]
    return BlockFunction(name, domain, Box{DT,T}(), block_size, blocks)
end
function BlockFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, blocks_size::SVector{DS,Int}) where {DS,S,DT,T}
    return BlockFunction{DS,S,DT,T}(name, domain, codomain, blocks_size)
end

# Metadata
name(bf::BlockFunction) = bf.name
domain(bf::BlockFunction) = bf.domain
codomain(bf::BlockFunction) = bf.codomain

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
function make_identity(bf::BlockFunction)
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

function project(bf::BlockFunction, cat::Category)
    @assert domain(bf) == domain(cat)
    blocks = map(b -> project(b, cat), bf.blocks)
    return BlockFunction{DS,S,DT,T}(cat.name, dom, cod, blocks)
end

function evaluate(bf::BlockFunction{DS,S,DT,T}, x::SVector{DS,S}) where {DS,S,DT,T}
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
evaluate(bf::BlockFunction{DS,S}, x::SVector{DS}) where {DS,S} = evaluate(bf, SVector{DS,S}(x))

integrate(bf::BlockFunction) = sum(integrate, bf.blocks)

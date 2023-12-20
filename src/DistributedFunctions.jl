module DistributedFunctions

using Distributed
using StaticArrays

using ..Boxes
using ..Categories
using ..Funs

struct DistributedFunction{DS,S,DT,T} <: Category{Box{DS,S},Box{DT,T}}
    name::AbstractString
    domain::Box{DS,S}
    codomain::Box{DT,T}
    future::Future           # returns `Category{Box{DS,S},Box{DT,T}}`
end
export DistributedFunction

function DistributedFunction{DS,S,DT,T}(make_cat, name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}) where {DS,S,DT,T}
    future = @spawnat :any make_cat(name, domain, codomain)
    return DistributedFunction(name, domain, Box{DT,T}(), future)
end
function DistributedFunction(make_cat, name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}) where {DS,S,DT,T}
    return DistributedFunction{DS,S,DT,T}(make_cat, name, domain, codomain)
end

Base.wait(df::DistributedFunction) = wait(df.future)
process(df::DistributedFunction) = df.future.where
fetch_future(df::DistributedFunction{DS,S,DT,T}) where {DS,S,DT,T} = fetch(df.future)::Category{Box{DS,S},Box{DT,T}}

function project(df::DistributedFunction, cat::Category)
    @assert domain(df) == domain(cat)
    future = @spawnat process(df) project(fetch_future(df), cat)
    return DistributedFunction(name(cat), domain(cat), codomain(cat), future)
end

# Metadata
Categories.name(df::DistributedFunction) = df.name
Categories.domain(df::DistributedFunction) = df.domain
Categories.codomain(df::DistributedFunction) = df.codomain

# Equality
function Base.:(==)(df1::DistributedFunction, df2::DistributedFunction)
    domain(df1) == domain(df2) || return false
    codomain(df1) == codomain(df2) || return false
    return fetch_future(df1) == fetch_future(df2)
end

# Collection is not available
Base.isempty(x::DistributedFunction) = isempty(fetch_future(x))
Base.length(x::DistributedFunction) = length(fetch_future(x))
Base.getindex(x::DistributedFunction, i) = getindex(fetch_future(x), i)
function Base.map(f::Fun{SVector{DT,T},SVector{DR,R}}, x::DistributedFunction{DS,S,DT,T}) where {DS,S,DT,T,DR,R}
    future = @spawnat process(x) map(f, fetch_future(x))
    # Use a pessimistic codomain
    codomain = Box{DR,R}(SVector{DR,R}(typemin(R) for d in 1:DR), SVector{DR,R}(typemax(R) for d in 1:DR))
    return DistributedFunction{DS,S,DR,R}(x.name, x.domain, codomain, future)
end
function Base.map(f::Category{Box{DT,T},Box{DR,R}}, x::DistributedFunction{DS,S,DT,T}) where {DS,S,DT,T,DR,R}
    @assert f.domain ⊇ x.codomain
    future = @spawnat process(x) map(f, fetch_future(x))
    return DistributedFunction{DS,S,DR,R}(x.name, x.domain, f.codomain, future)
end
function Base.map(
    f::Fun{Tuple{SVector{DT,T},SVector{DU,U}},SVector{DR,R}}, x::DistributedFunction{DS,S,DT,T}, y::DistributedFunction{DS,S,DU,U}
) where {DS,S,DT,T,DU,U,DR,R}
    # We arbitrarily choose the process of the first argument
    future = @spawnat process(x) map(f, fetch_future(x), fetch_future(y))
    # Use a pessimistic codomain
    codomain = Box{DR,R}(SVector{DR,R}(typemin(R) for d in 1:DR), SVector{DR,R}(typemax(R) for d in 1:DR))
    return DistributedFunction{DS,S,DR,R}(x.name, x.domain, codomain, future)
end

# Category
function Categories.make_identity(df::DistributedFunction)
    future = @spawnat process(df) make_identity(fetch_future(df))
    return DistributedFunction(name(df), domain(df), domain(df), future)
end
function Base.:∘(df2::DistributedFunction{<:Any,<:Any,DT}, df1::DistributedFunction{DT}) where {DT}
    # We arbitrarily choose the process of the first function
    @assert domain(df2) == codomain(df1)
    future = @spawnat process(df1) fetch_future(df2) ∘ fetch_future(df1)
    return DistributedFunction(name(df1), domain(df1), codomain(df2), future)
end

# Vector space
function Base.zero(x::DistributedFunction)
    cod = codomain(x)
    cod′ = Box(zero(first(cod)), zero(last(cod)))
    future = @spawnat process(x) zero(fetch_future(x))
    return DistributedFunction(name(x), domain(x), cod′, future)
end
function Base.:+(x::DistributedFunction)
    cod = codomain(x)
    cod′ = Box(first(cod), last(cod))
    future = @spawnat process(x) +(fetch_future(x))
    return DistributedFunction(name(x), domain(x), cod′, future)
end
function Base.:-(x::DistributedFunction)
    cod = codomain(x)
    cod′ = Box(-last(cod), -first(cod))
    future = @spawnat process(x) -(fetch_future(x))
    return DistributedFunction(name(x), domain(x), cod′, future)
end
function Base.:+(x::DistributedFunction, y::DistributedFunction)
    @assert domain(x) == domain(y)
    @assert ndims(codomain(x)) == ndims(codomain(y))
    codx = codomain(x)
    cody = codomain(y)
    cod′ = Box(min.(first(codx), first(cody)), max.(last(codx), last(cody)))
    future = @spawnat process(x) fetch_future(x) + fetch_future(y)
    return DistributedFunction(name(x), domain(x), cod′, future)
end
function Base.:-(x::DistributedFunction, y::DistributedFunction)
    @assert domain(x) == domain(y)
    @assert ndims(codomain(x)) == ndims(codomain(y))
    codx = codomain(x)
    cody = codomain(y)
    cod′ = Box(min.(first(codx), -last(cody)), max.(last(codx), -first(cody)))
    future = @spawnat process(x) fetch_future(x) - fetch_future(y)
    return DistributedFunction(name(x), domain(x), cod′, future)
end
function Base.:*(a, x::DistributedFunction)
    cod = codomain(x)
    if a < 0
        cod′ = Box(a * last(cody), a * first(codx))
    elseif a > 0
        cod′ = Box(a * first(codx), a * last(cody))
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    future = @spawnat process(x) a * fetch_future(x)
    return DistributedFunction(name(x), domain(x), cod′, future)
end
function Base.:*(x::DistributedFunction, a)
    cod = codomain(x)
    if a < 0
        cod′ = Box(last(cody) * a, first(codx) * a)
    elseif a > 0
        cod′ = Box(first(codx) * a, last(cody) * a)
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    future = @spawnat process(x) fetch_future(x) * a
    return DistributedFunction(name(x), domain(x), cod′, future)
end
function Base.:\(a, x::DistributedFunction)
    cod = codomain(x)
    if a < 0
        cod′ = Box(a \ last(cody), a \ first(codx))
    elseif a > 0
        cod′ = Box(a \ first(codx), a \ last(cody))
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    future = @spawnat process(x) a \ fetch_future(x)
    return DistributedFunction(name(x), domain(x), cod′, future)
end
function Base.:/(x::DistributedFunction, a)
    cod = codomain(x)
    if a < 0
        cod′ = Box(last(cody) / a, first(codx) / a)
    elseif a > 0
        cod′ = Box(first(codx) / a, last(cody) / a)
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    future = @spawnat process(x) fetch_future(x) / a
    return DistributedFunction(name(x), domain(x), cod′, future)
end

Categories.evaluate(df::DistributedFunction{DS,S,DT,T}, x::SVector{DS,S}) where {DS,S,DT,T} = fetch_future(df)(x)::SVector{DT,T}
Categories.evaluate(df::DistributedFunction{DS,S,DT,T}, x::SVector{DS}) where {DS,S,DT,T} = evaluate(df, SVector{DS,S}(x))

Categories.integrate(df::DistributedFunction) = integrate(fetch_future(df))

end

module ThreadedFunctions

using Base.Threads
using StaticArrays

using ..Boxes
using ..Categories
using ..Funs

struct ThreadedFunction{DS,S,DT,T} <: Category{Box{DS,S},Box{DT,T}}
    name::AbstractString
    domain::Box{DS,S}
    codomain::Box{DT,T}
    task::Task               # returns `Category{Box{DS,S},Box{DT,T}}`
end
export ThreadedFunction

function ThreadedFunction{DS,S,DT,T}(make_cat, name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}) where {DS,S,DT,T}
    task = @spawn make_cat(name, domain, codomain)
    return ThreadedFunction(name, domain, Box{DT,T}(), task)
end
function ThreadedFunction(make_cat, name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}) where {DS,S,DT,T}
    return ThreadedFunction{DS,S,DT,T}(make_cat, name, domain, codomain)
end

Base.wait(tf::ThreadedFunction) = wait(tf.task)
fetch_task(tf::ThreadedFunction{DS,S,DT,T}) where {DS,S,DT,T} = fetch(tf.task)::Category{Box{DS,S},Box{DT,T}}

function project(tf::ThreadedFunction, cat::Category)
    @assert domain(tf) == domain(cat)
    task = @spawn project(fetch_task(tf), cat)
    return ThreadedFunction(name(cat), domain(cat), codomain(cat), task)
end

# Metadata
Categories.name(tf::ThreadedFunction) = tf.name
Categories.domain(tf::ThreadedFunction) = tf.domain
Categories.codomain(tf::ThreadedFunction) = tf.codomain

# Equality
function Base.:(==)(tf1::ThreadedFunction, tf2::ThreadedFunction)
    domain(tf1) == domain(tf2) || return false
    codomain(tf1) == codomain(tf2) || return false
    return fetch_task(tf1) == fetch_task(tf2)
end

# Collection is not available
Base.isempty(x::ThreadedFunction) = isempty(fetch_task(x))
Base.length(x::ThreadedFunction) = length(fetch_task(x))
Base.getindex(x::ThreadedFunction, i) = getindex(fetch_task(x), i)
function Base.map(f::Fun{SVector{DT,T},SVector{DR,R}}, x::ThreadedFunction{DS,S,DT,T}) where {DS,S,DT,T,DR,R}
    task = @spawn map(f, fetch_task(x))
    # Use a pessimistic codomain
    codomain = Box{DR,R}(SVector{DR,R}(typemin(R) for d in 1:DR), SVector{DR,R}(typemax(R) for d in 1:DR))
    return ThreadedFunction{DS,S,DR,R}(x.name, x.domain, codomain, task)
end
function Base.map(f::Category{Box{DT,T},Box{DR,R}}, x::ThreadedFunction{DS,S,DT,T}) where {DS,S,DT,T,DR,R}
    @assert f.domain ⊇ x.codomain
    task = @spawn map(f, fetch_task(x))
    return ThreadedFunction{DS,S,DR,R}(x.name, x.domain, f.codomain, task)
end
function Base.map(
    f::Fun{Tuple{SVector{DT,T},SVector{DU,U}},SVector{DR,R}}, x::ThreadedFunction{DS,S,DT,T}, y::ThreadedFunction{DS,S,DU,U}
) where {DS,S,DT,T,DU,U,DR,R}
    task = @spawn map(f, fetch_task(x), fetch_task(y))
    # Use a pessimistic codomain
    codomain = Box{DR,R}(SVector{DR,R}(typemin(R) for d in 1:DR), SVector{DR,R}(typemax(R) for d in 1:DR))
    return ThreadedFunction{DS,S,DR,R}(x.name, x.domain, codomain, task)
end

# Category
function Categories.make_identity(tf::ThreadedFunction)
    task = @spawn make_identity(fetch_task(tf))
    return ThreadedFunction(name(tf), domain(tf), domain(tf), task)
end
function Base.:∘(tf2::ThreadedFunction{<:Any,<:Any,DT}, tf1::ThreadedFunction{DT}) where {DT}
    @assert domain(tf2) == codomain(tf1)
    task = @spawn fetch_task(tf2) ∘ fetch_task(tf1)
    return ThreadedFunction(name(tf1), domain(tf1), codomain(tf2), task)
end

# Vector space
function Base.zero(x::ThreadedFunction)
    cod = codomain(x)
    cod′ = Box(zero(first(cod)), zero(last(cod)))
    task = @spawn zero(fetch_task(x))
    return ThreadedFunction(name(x), domain(x), cod′, task)
end
function Base.:+(x::ThreadedFunction)
    cod = codomain(x)
    cod′ = Box(first(cod), last(cod))
    task = @spawn +(fetch_task(x))
    return ThreadedFunction(name(x), domain(x), cod′, task)
end
function Base.:-(x::ThreadedFunction)
    cod = codomain(x)
    cod′ = Box(-last(cod), -first(cod))
    task = @spawn -(fetch_task(x))
    return ThreadedFunction(name(x), domain(x), cod′, task)
end
function Base.:+(x::ThreadedFunction, y::ThreadedFunction)
    @assert domain(x) == domain(y)
    @assert ndims(codomain(x)) == ndims(codomain(y))
    codx = codomain(x)
    cody = codomain(y)
    cod′ = Box(min.(first(codx), first(cody)), max.(last(codx), last(cody)))
    task = @spawn fetch_task(x) + fetch_task(y)
    return ThreadedFunction(name(x), domain(x), cod′, task)
end
function Base.:-(x::ThreadedFunction, y::ThreadedFunction)
    @assert domain(x) == domain(y)
    @assert ndims(codomain(x)) == ndims(codomain(y))
    codx = codomain(x)
    cody = codomain(y)
    cod′ = Box(min.(first(codx), -last(cody)), max.(last(codx), -first(cody)))
    task = @spawn fetch_task(x) - fetch_task(y)
    return ThreadedFunction(name(x), domain(x), cod′, task)
end
function Base.:*(a, x::ThreadedFunction)
    cod = codomain(x)
    if a < 0
        cod′ = Box(a * last(cody), a * first(codx))
    elseif a > 0
        cod′ = Box(a * first(codx), a * last(cody))
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    task = @spawn a * fetch_task(x)
    return ThreadedFunction(name(x), domain(x), cod′, task)
end
function Base.:*(x::ThreadedFunction, a)
    cod = codomain(x)
    if a < 0
        cod′ = Box(last(cody) * a, first(codx) * a)
    elseif a > 0
        cod′ = Box(first(codx) * a, last(cody) * a)
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    task = @spawn fetch_task(x) * a
    return ThreadedFunction(name(x), domain(x), cod′, task)
end
function Base.:\(a, x::ThreadedFunction)
    cod = codomain(x)
    if a < 0
        cod′ = Box(a \ last(cody), a \ first(codx))
    elseif a > 0
        cod′ = Box(a \ first(codx), a \ last(cody))
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    task = @spawn a \ fetch_task(x)
    return ThreadedFunction(name(x), domain(x), cod′, task)
end
function Base.:/(x::ThreadedFunction, a)
    cod = codomain(x)
    if a < 0
        cod′ = Box(last(cody) / a, first(codx) / a)
    elseif a > 0
        cod′ = Box(first(codx) / a, last(cody) / a)
    else
        cod′ = Box(zero(first(codx)), zero(last(codx)))
    end
    task = @spawn fetch_task(x) / a
    return ThreadedFunction(name(x), domain(x), cod′, task)
end

Categories.evaluate(tf::ThreadedFunction{DS,S,DT,T}, x::SVector{DS,S}) where {DS,S,DT,T} = fetch_task(tf)(x)::SVector{DT,T}
Categories.evaluate(tf::ThreadedFunction{DS,S,DT,T}, x::SVector{DS}) where {DS,S,DT,T} = evaluate(tf, SVector{DS,S}(x))

Categories.integrate(tf::ThreadedFunction) = integrate(fetch_task(tf))

end

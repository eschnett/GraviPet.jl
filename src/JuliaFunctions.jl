module JuliaFunctions

using HCubature
using StaticArrays

using ..Boxes
using ..Categories

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
Categories.name(jf::JuliaFunction) = jf.name
Categories.domain(jf::JuliaFunction) = jf.domain
Categories.codomain(jf::JuliaFunction) = jf.codomain

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
Categories.make_identity(jf::JuliaFunction) = JuliaFunction(jf.name, jf.dom, jf.dom, identity)
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

Categories.evaluate(jf::JuliaFunction{DS,S,DT,T}, x::SVector{DS,S}) where {DS,S,DT,T} = jf.fun(x)::SVector{DT,T}
Categories.evaluate(jf::JuliaFunction{DS,S,DT,T}, x::SVector{DS}) where {DS,S,DT,T} = evaluate(jf, SVector{DS,S}(x))

function Categories.integrate(jf::JuliaFunction)
    dom = domain(jf)
    cod = codomain(jf)
    I, E = hcubature(jf.fun, first(dom), last(dom))
    return I::eltype(cod)
end

end

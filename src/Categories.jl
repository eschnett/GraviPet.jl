module Categories

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

end

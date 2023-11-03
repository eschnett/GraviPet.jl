var documenterSearchIndex = {"docs":
[{"location":"#GraviPet.jl","page":"GraviPet.jl","title":"GraviPet.jl","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet]","category":"page"},{"location":"#GraviPet.Box","page":"GraviPet.jl","title":"GraviPet.Box","text":"struct Box{D,T} <: Domain{SVector{D,T}}\n\nA multi-dimensional Domain, i.e. a box or hyper-rectangle. The number of dimensions D must be a non-negative integer. T should be a real-valued scalar type such as Float64. The domain's element type is SVector{D,T}.\n\nZero-dimensional and one-dimensional domains are supported. Zero-dimensional domains contain only a single point (the origin). One-dimensional domains correspond to intervals (see Interval).\n\n\n\n\n\n","category":"type"},{"location":"#GraviPet.Box-Union{Tuple{T}, Tuple{D}, Tuple{StaticArraysCore.SVector{D, var\"#s6\"} where var\"#s6\"<:Real, StaticArraysCore.SVector{D, var\"#s7\"} where var\"#s7\"<:Real}} where {D, T}","page":"GraviPet.jl","title":"GraviPet.Box","text":"Box{D,T}(first::SVector{D}, last::SVector{D})\nBox(first::SVector{D}, last::SVector{D})\n\nCreate a box domain with the specified lower and upper bounds. The bounds are given as vectors, and (some of) their elements can be infinity. The constraint first < last must hold in each dimension.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Domain","page":"GraviPet.jl","title":"GraviPet.Domain","text":"abstract type Domain{T}\n\nA computational domain of type T. The domain of a function is the type of its argument(s). For example, Domain{Float64} would describe a one-dimensional domain, and Domain{SVector{3,Float64}} would describe a three-dimensional domain.\n\nThis type can also describe a codomain, i.e. the result type of a function. A codomain of Domain{Float64} describes a scalar function, and a codomain of Domain{SVector{3,Float64}} describes a vector-valued function.\n\nSee also Interval and Box for concrete domain types.\n\n\n\n\n\n","category":"type"},{"location":"#GraviPet.Interval","page":"GraviPet.jl","title":"GraviPet.Interval","text":"struct Interval{T} <: Domain{T}\n\nA one-dimensional Domain, i.e. an interval. T should be a real-valued scalar type such as Float64.\n\n\n\n\n\n","category":"type"},{"location":"#GraviPet.Interval-Union{Tuple{T}, Tuple{Real, Real}} where T","page":"GraviPet.jl","title":"GraviPet.Interval","text":"Interval{T}(first::Real, last::Real)\nInterval(first::Real, last::Real)\n\nCreate an interval with the specified bounds. The bounds can be infinity. The constraint first < last must hold.\n\n\n\n\n\n","category":"method"},{"location":"#Base.:==-Tuple{Domain, Domain}","page":"GraviPet.jl","title":"Base.:==","text":"==(dom1::Domain, dom2::Domain)::Bool\n\nEquality for domains. Domains are equal if the have the same number of dimensions (see ndims), and the same lower and upper bounds (see first and last). The element types might be different if they can be compared for equality, e.g. Float64 and BigFloat.\n\nIncompatible domains (living in different number of dimensions) are treated as unequal.\n\n\n\n\n\n","category":"method"},{"location":"#Base.eltype-Union{Tuple{Domain{T}}, Tuple{T}} where T","page":"GraviPet.jl","title":"Base.eltype","text":"eltype(dom::Domain)::Type\n\nReturn the element type of a Domain.\n\n\n\n\n\n","category":"method"},{"location":"#Base.first-Tuple{Domain}","page":"GraviPet.jl","title":"Base.first","text":"first(dom::Domain{T})::T\n\nReturn the lower bound of the domain. (At the moment domains are assumed to be cuboid, i.e. to form a box or hyper-rectangle.) Domain bounds can be infinity.\n\nSee also last.\n\n\n\n\n\n","category":"method"},{"location":"#Base.in-Tuple{Any, Domain}","page":"GraviPet.jl","title":"Base.in","text":"in(elem, dom::Domain)::Bool\n\nTest whether elem is inside domain dom.\n\n\n\n\n\n","category":"method"},{"location":"#Base.issubset-Tuple{Domain, Domain}","page":"GraviPet.jl","title":"Base.issubset","text":"issubset(dom1::Domain, dom2::Domain)::Bool\n\nTest whether dom1 is a subset of dom2. Only compatible domains (which live in the same number of dimensions) can be compared.\n\n\n\n\n\n","category":"method"},{"location":"#Base.last-Tuple{Domain}","page":"GraviPet.jl","title":"Base.last","text":"last(dom::Domain{T})::T\n\nReturn the upper bound of the domain. (At the moment domains are assumed to be cuboid, i.e. to form a box or hyper-rectangle.) Domain bounds can be infinity.\n\nSee also first.\n\n\n\n\n\n","category":"method"},{"location":"#Base.ndims-Tuple{Domain}","page":"GraviPet.jl","title":"Base.ndims","text":"ndims(dom::Domain)::Int\n\nReturn the number of dimensions of a Domain. Scalar domains are one-dimensional. Zero-dimensional domains consist of a single point only.\n\n\n\n\n\n","category":"method"}]
}

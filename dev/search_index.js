var documenterSearchIndex = {"docs":
[{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"(Image: GraviPet logo)","category":"page"},{"location":"#GraviPet.jl","page":"GraviPet.jl","title":"GraviPet.jl","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet]","category":"page"},{"location":"#GraviPet.Box","page":"GraviPet.jl","title":"GraviPet.Box","text":"struct Box{D,T} <: Domain{SVector{D,T}}\n\nA multi-dimensional Domain, i.e. a box or hyper-rectangle. The number of dimensions D must be a non-negative integer. T should be a real-valued scalar type such as Float64. The domain's element type is SVector{D,T}.\n\nZero-dimensional and one-dimensional domains are supported. Zero-dimensional domains contain only a single point (the origin). One-dimensional domains correspond to intervals (see Interval).\n\n\n\n\n\n","category":"type"},{"location":"#GraviPet.Box-Union{Tuple{T}, Tuple{D}, Tuple{StaticArraysCore.SVector{D, var\"#s6\"} where var\"#s6\"<:Real, StaticArraysCore.SVector{D, var\"#s7\"} where var\"#s7\"<:Real}} where {D, T}","page":"GraviPet.jl","title":"GraviPet.Box","text":"Box{D,T}(first::SVector{D}, last::SVector{D})\nBox(first::SVector{D}, last::SVector{D})\n\nCreate a box domain with the specified lower and upper bounds. The bounds are given as vectors, and (some of) their elements can be infinity. The constraint first < last must hold in each dimension.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Box-Union{Tuple{}, Tuple{T}, Tuple{D}} where {D, T}","page":"GraviPet.jl","title":"GraviPet.Box","text":"Box{D,T}()\n\nCreate a zero-valued box domain, i.e. a box where both lower and upper bounds are zero. This is also a convenient way to provide a \"skeleton\" for a domain without having to specify any bounds.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Category-Tuple{Any}","page":"GraviPet.jl","title":"GraviPet.Category","text":"(cat::Category)(x)::eltype(codomain(cat))\n\nEvaluate the function cat at point x. This is written as a function call cat(x). See also evaluate.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Domain","page":"GraviPet.jl","title":"GraviPet.Domain","text":"abstract type Domain{T}\n\nA computational domain of type T. The domain of a function is the type of its argument(s). For example, Domain{Float64} would describe a one-dimensional domain, and Domain{SVector{3,Float64}} would describe a three-dimensional domain.\n\nThis type can also describe a codomain, i.e. the result type of a function. A codomain of Domain{Float64} describes a scalar function, and a codomain of Domain{SVector{3,Float64}} describes a vector-valued function.\n\nSee also Interval and Box for concrete domain types.\n\n\n\n\n\n","category":"type"},{"location":"#GraviPet.GridFunction-Union{Tuple{T}, Tuple{DT}, Tuple{S}, Tuple{DS}, Tuple{AbstractString, Box{DS, S}, Box{DT, T}, StaticArraysCore.SVector{DS, Int64}}} where {DS, S, DT, T}","page":"GraviPet.jl","title":"GraviPet.GridFunction","text":"GridFunction{DS,S,DT,T}(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_shape::SVector{DS,Int})\nGridFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_shape::SVector{DS,Int})\n\nCreate a zero-valued grid function, i.e. a grid function that is zero everywhere. The codomain can be a zero-valued domain created via Box{DT,T}().\n\ngrid_shape specifies the resolution of the grid function. This makes it convenient to use zero-valued grid functions as \"templates\" or \"skeletons\" when creating other grid functions.\n\nThis zero-valued grid function is represented efficiently and does not store one element per grid point.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.GridFunction1D-Union{Tuple{T}, Tuple{S}, Tuple{AbstractString, Interval{S}, Interval{T}, Int64}} where {S, T}","page":"GraviPet.jl","title":"GraviPet.GridFunction1D","text":"GridFunction1D{S,T}(name::AbstractString, domain::Interval{S}, codomain::Interval{T}, npoints::Int)\nGridFunction1D(name::AbstractString, domain::Interval{S}, codomain::Interval{T}, npoints::Int)\n\nCreate a zero-valued grid function, i.e. a grid function that is zero everywhere. The codomain can be a zero-valued domain created via Interval{T}().\n\nnpoints specifies the resolution of the grid function. This makes it convenient to use zero-valued grid functions as \"templates\" or \"skeletons\" when creating other grid functions.\n\nThis zero-valued grid function is represented efficiently and does not store one element per grid point.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Interval","page":"GraviPet.jl","title":"GraviPet.Interval","text":"struct Interval{T} <: Domain{T}\n\nA one-dimensional Domain, i.e. an interval. T should be a real-valued scalar type such as Float64.\n\n\n\n\n\n","category":"type"},{"location":"#GraviPet.Interval-Union{Tuple{T}, Tuple{Real, Real}} where T","page":"GraviPet.jl","title":"GraviPet.Interval","text":"Interval{T}(first::Real, last::Real)\nInterval(first::Real, last::Real)\n\nCreate an interval with the specified bounds. The bounds can be infinity. The constraint first < last must hold.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Interval-Union{Tuple{}, Tuple{T}} where T","page":"GraviPet.jl","title":"GraviPet.Interval","text":"Interval{T}()\n\nCreate a zero-valued interval, i.e. an interval where both lower and upper bounds are zero. This is also a convenient way to provide a \"skeleton\" for a domain without having to specify any bounds.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.JuliaFunction-Union{Tuple{T}, Tuple{DT}, Tuple{S}, Tuple{DS}, Tuple{AbstractString, Box{DS, S}, Box{DT, T}}} where {DS, S, DT, T}","page":"GraviPet.jl","title":"GraviPet.JuliaFunction","text":"JuliaFunction{DS,S,DT,T}(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T})\nJuliaFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T})\n\nCreate a zero-valued Julia function, i.e. a Julia function that is zero everywhere. The codomain can be a zero-valued domain created via Box{DT,T}().\n\n\n\n\n\n","category":"method"},{"location":"#Base.:==-Tuple{Domain, Domain}","page":"GraviPet.jl","title":"Base.:==","text":"==(dom1::Domain, dom2::Domain)::Bool\n\nEquality for domains. Domains are equal if the have the same number of dimensions (see ndims), and the same lower and upper bounds (see first and last). The element types might be different if they can be compared for equality, e.g. Float64 and BigFloat.\n\nIncompatible domains (living in different number of dimensions) are treated as unequal.\n\n\n\n\n\n","category":"method"},{"location":"#Base.eltype-Union{Tuple{Domain{T}}, Tuple{T}} where T","page":"GraviPet.jl","title":"Base.eltype","text":"eltype(dom::Domain)::Type\n\nReturn the element type of a Domain.\n\n\n\n\n\n","category":"method"},{"location":"#Base.first-Tuple{Domain}","page":"GraviPet.jl","title":"Base.first","text":"first(dom::Domain{T})::T\n\nReturn the lower bound of the domain. (At the moment domains are assumed to be cuboid, i.e. to form a box or hyper-rectangle.) Domain bounds can be infinity.\n\nSee also last.\n\n\n\n\n\n","category":"method"},{"location":"#Base.in-Tuple{Any, Domain}","page":"GraviPet.jl","title":"Base.in","text":"in(elem, dom::Domain)::Bool\n\nTest whether elem is inside domain dom.\n\n\n\n\n\n","category":"method"},{"location":"#Base.isdisjoint-Tuple{Domain, Domain}","page":"GraviPet.jl","title":"Base.isdisjoint","text":"isdisjoint(dom1::Domain, dom2::Domain)::Bool\n\nTest whether dom1 is disjoint from dom2. Only compatible domains (which live in the same number of dimensions) can be compared.\n\n\n\n\n\n","category":"method"},{"location":"#Base.issubset-Tuple{Domain, Domain}","page":"GraviPet.jl","title":"Base.issubset","text":"issubset(dom1::Domain, dom2::Domain)::Bool\n\nTest whether dom1 is a subset of dom2. Only compatible domains (which live in the same number of dimensions) can be compared.\n\n\n\n\n\n","category":"method"},{"location":"#Base.last-Tuple{Domain}","page":"GraviPet.jl","title":"Base.last","text":"last(dom::Domain{T})::T\n\nReturn the upper bound of the domain. (At the moment domains are assumed to be cuboid, i.e. to form a box or hyper-rectangle.) Domain bounds can be infinity.\n\nSee also first.\n\n\n\n\n\n","category":"method"},{"location":"#Base.ndims-Tuple{Domain}","page":"GraviPet.jl","title":"Base.ndims","text":"ndims(dom::Domain)::Int\n\nReturn the number of dimensions of a Domain. Scalar domains are one-dimensional. Zero-dimensional domains consist of a single point only.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.evaluate-Tuple{Category, Any}","page":"GraviPet.jl","title":"GraviPet.evaluate","text":"evaluate(cat::Category, x)::eltype(codomain(cat))\n\nEvaluate the function cat at point x. This can also be written as a function call cat(x).\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.project-Tuple{Category, Category}","page":"GraviPet.jl","title":"GraviPet.project","text":"project(dst::T, src::Category)::T where {T<:Category}\n\nProject the function src onto the representation given by dst. This can be used e.g. to sample Julia functions on a given grid, or to resample functiont on a different grid. dst can be a skeleton representation.\n\n\n\n\n\n","category":"method"}]
}

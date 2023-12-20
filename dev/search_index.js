var documenterSearchIndex = {"docs":
[{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"(Image: GraviPet logo)","category":"page"},{"location":"#GraviPet.jl","page":"GraviPet.jl","title":"GraviPet.jl","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"GraviPet is the General Relativistic Astrophysics Visualization, Initialization, and Postprocessing Efficient Toolkit.","category":"page"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet]\nPrivate = false","category":"page"},{"location":"#Domains","page":"GraviPet.jl","title":"Domains","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.Domains]\nPrivate = false","category":"page"},{"location":"#GraviPet.Domains.Domain","page":"GraviPet.jl","title":"GraviPet.Domains.Domain","text":"abstract type Domain{T}\n\nA computational domain of type T. The domain of a function is the type of its argument(s). For example, Domain{Float64} would describe a one-dimensional domain, and Domain{SVector{3,Float64}} would describe a three-dimensional domain.\n\nThis type can also describe a codomain, i.e. the result type of a function. A codomain of Domain{Float64} describes a scalar function, and a codomain of Domain{SVector{3,Float64}} describes a vector-valued function.\n\nSee also Interval and Box for concrete domain types.\n\n\n\n\n\n","category":"type"},{"location":"#GraviPet.Domains.expanded-Union{Tuple{T}, Tuple{Domain{T}, T, T}} where T","page":"GraviPet.jl","title":"GraviPet.Domains.expanded","text":"expanded(dom::Domain, delta)::Domain\nexpanded(dom::Domain, delta_lo, delta_hi)::Domain\n\nExpand the domain dom. The lower and upper domain bounds are moved by either delta, or delta_lo and delta_hi, respectively. Positive deltas enlarge the domain, negative values shrink it. Domains must be non-empty.\n\nexpanded(dom, delta) is equivalent to expanded(dom, delta, delta).\n\nSee also shifted.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Domains.shifted-Tuple{Domain, Any}","page":"GraviPet.jl","title":"GraviPet.Domains.shifted","text":"shifted(dom::Domain, shift)::Domain\n\nShift the domain dom by moving both domain bounds in the same direction. This does not change the size of the domain.\n\nshfited(dom, shift) is equivalent to expanded(dom, -shift, +shift).\n\nSee also expanded.\n\n\n\n\n\n","category":"method"},{"location":"#Intervals","page":"GraviPet.jl","title":"Intervals","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.Intervals]\nPrivate = false","category":"page"},{"location":"#GraviPet.Intervals.Interval","page":"GraviPet.jl","title":"GraviPet.Intervals.Interval","text":"struct Interval{T} <: Domain{T}\n\nA one-dimensional Domain, i.e. an interval. T should be a real-valued scalar type such as Float64.\n\n\n\n\n\n","category":"type"},{"location":"#GraviPet.Intervals.Interval-Union{Tuple{T}, Tuple{Real, Real}} where T","page":"GraviPet.jl","title":"GraviPet.Intervals.Interval","text":"Interval{T}(first::Real, last::Real)\nInterval(first::Real, last::Real)\n\nCreate an interval with the specified bounds. The bounds can be infinity. The constraint first < last must hold.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Intervals.Interval-Union{Tuple{}, Tuple{T}} where T","page":"GraviPet.jl","title":"GraviPet.Intervals.Interval","text":"Interval{T}()\n\nCreate a zero-valued interval, i.e. an interval where both lower and upper bounds are zero. This is also a convenient way to provide a \"skeleton\" for a domain without having to specify any bounds.\n\n\n\n\n\n","category":"method"},{"location":"#Boxes","page":"GraviPet.jl","title":"Boxes","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.Boxes]\nPrivate = false","category":"page"},{"location":"#GraviPet.Boxes.Box","page":"GraviPet.jl","title":"GraviPet.Boxes.Box","text":"struct Box{D,T} <: Domain{SVector{D,T}}\n\nA multi-dimensional Domains.Domain, i.e. a box or hyper-rectangle.  The number of dimensions D must be a non-negative integer. T should be a real-valued scalar type such as Float64. The domain's element type is SVector{D,T}.\n\nZero-dimensional and one-dimensional domains are supported. Zero-dimensional domains contain only a single point (the origin). One-dimensional domains correspond to intervals (see Interval).\n\n\n\n\n\n","category":"type"},{"location":"#GraviPet.Boxes.Box-Union{Tuple{T}, Tuple{D}, Tuple{StaticArraysCore.SVector{D, <:Real}, StaticArraysCore.SVector{D, <:Real}}} where {D, T}","page":"GraviPet.jl","title":"GraviPet.Boxes.Box","text":"Box{D,T}(first::SVector{D}, last::SVector{D})\nBox(first::SVector{D}, last::SVector{D})\n\nCreate a box domain with the specified lower and upper bounds. The bounds are given as vectors, and (some of) their elements can be infinity. The constraint first < last must hold in each dimension.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Boxes.Box-Union{Tuple{}, Tuple{T}, Tuple{D}} where {D, T}","page":"GraviPet.jl","title":"GraviPet.Boxes.Box","text":"Box{D,T}()\n\nCreate a zero-valued box domain, i.e. a box where both lower and upper bounds are zero. This is also a convenient way to provide a \"skeleton\" for a domain without having to specify any bounds.\n\n\n\n\n\n","category":"method"},{"location":"#Categories","page":"GraviPet.jl","title":"Categories","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.Categories]\nPrivate = false","category":"page"},{"location":"#GraviPet.Categories.Category-Tuple{Any}","page":"GraviPet.jl","title":"GraviPet.Categories.Category","text":"(cat::Category)(x)::eltype(codomain(cat))\n\nEvaluate the function cat at point x. This is written as a function call cat(x). See also evaluate.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Categories.evaluate-Tuple{Category, Any}","page":"GraviPet.jl","title":"GraviPet.Categories.evaluate","text":"evaluate(cat::Category, x)::eltype(codomain(cat))\n\nEvaluate the function cat at point x. This can also be written as a function call cat(x).\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Categories.integrate-Tuple{Category}","page":"GraviPet.jl","title":"GraviPet.Categories.integrate","text":"integrate(cat::Category)::eltype(codomain(cat))\n\nIntegrate the function cat over its domain.\n\n\n\n\n\n","category":"method"},{"location":"#GraviPet.Categories.project-Tuple{Category, Category}","page":"GraviPet.jl","title":"GraviPet.Categories.project","text":"project(dst::T, src::Category)::T where {T<:Category}\n\nProject the function src onto the representation given by dst. This can be used e.g. to sample Julia functions on a given grid, or to resample functiont on a different grid. dst can be a skeleton representation.\n\n\n\n\n\n","category":"method"},{"location":"#Julia-Functions","page":"GraviPet.jl","title":"Julia Functions","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.JuliaFunctions]\nPrivate = false","category":"page"},{"location":"#GraviPet.JuliaFunctions.JuliaFunction-Union{Tuple{T}, Tuple{DT}, Tuple{S}, Tuple{DS}, Tuple{AbstractString, Box{DS, S}, Box{DT, T}}} where {DS, S, DT, T}","page":"GraviPet.jl","title":"GraviPet.JuliaFunctions.JuliaFunction","text":"JuliaFunction{DS,S,DT,T}(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T})\nJuliaFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T})\n\nCreate a zero-valued Julia function, i.e. a Julia function that is zero everywhere. The codomain can be a zero-valued domain created via Box{DT,T}().\n\n\n\n\n\n","category":"method"},{"location":"#D-Grid-Functions","page":"GraviPet.jl","title":"1D Grid Functions","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.GridFunctions1D]\nPrivate = false","category":"page"},{"location":"#GraviPet.GridFunctions1D.GridFunction1D-Union{Tuple{T}, Tuple{S}, Tuple{AbstractString, Interval{S}, Interval{T}, Int64}} where {S, T}","page":"GraviPet.jl","title":"GraviPet.GridFunctions1D.GridFunction1D","text":"GridFunction1D{S,T}(name::AbstractString, domain::Interval{S}, codomain::Interval{T}, npoints::Int)\nGridFunction1D(name::AbstractString, domain::Interval{S}, codomain::Interval{T}, npoints::Int)\n\nCreate a zero-valued grid function, i.e. a grid function that is zero everywhere. The codomain can be a zero-valued domain created via Interval{T}().\n\nnpoints specifies the resolution of the grid function. This makes it convenient to use zero-valued grid functions as \"templates\" or \"skeletons\" when creating other grid functions.\n\nThis zero-valued grid function is represented efficiently and does not store one element per grid point.\n\n\n\n\n\n","category":"method"},{"location":"#Multi-dimensional-Grid-Functions","page":"GraviPet.jl","title":"Multi-dimensional Grid Functions","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.GridFunctions]\nPrivate = false","category":"page"},{"location":"#GraviPet.GridFunctions.GridFunction-Union{Tuple{T}, Tuple{DT}, Tuple{S}, Tuple{DS}, Tuple{AbstractString, Box{DS, S}, Box{DT, T}, StaticArraysCore.SVector{DS, Int64}}} where {DS, S, DT, T}","page":"GraviPet.jl","title":"GraviPet.GridFunctions.GridFunction","text":"GridFunction{DS,S,DT,T}(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_size::SVector{DS,Int})\nGridFunction(name::AbstractString, domain::Box{DS,S}, codomain::Box{DT,T}, grid_size::SVector{DS,Int})\n\nCreate a zero-valued grid function, i.e. a grid function that is zero everywhere. The codomain can be a zero-valued domain created via Box{DT,T}().\n\ngrid_size specifies the resolution of the grid function. This makes it convenient to use zero-valued grid functions as \"templates\" or \"skeletons\" when creating other grid functions.\n\nThis zero-valued grid function is represented efficiently and does not store one element per grid point.\n\n\n\n\n\n","category":"method"},{"location":"#Blocked-Functions-(Domain-Decomposition)","page":"GraviPet.jl","title":"Blocked Functions (Domain Decomposition)","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.BlockFunctions]\nPrivate = false","category":"page"},{"location":"#Threaded-Functions","page":"GraviPet.jl","title":"Threaded Functions","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.ThreadedFunctions]\nPrivate = false","category":"page"},{"location":"#Distributed-Functions","page":"GraviPet.jl","title":"Distributed Functions","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.DistributedFunctions]\nPrivate = false","category":"page"},{"location":"#Miscellaneous","page":"GraviPet.jl","title":"Miscellaneous","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Base.:(==)(::Domain, ::Domain)\nBase.last(::Domain)\nBase.isdisjoint(::Domain, ::Domain)\nBase.eltype(::Domain)\nBase.ndims(::Domain)\nBase.issubset(::Domain, ::Domain)\nBase.in(::Any, ::Domain)\nBase.first(::Domain)","category":"page"},{"location":"#Base.:==-Tuple{Domain, Domain}","page":"GraviPet.jl","title":"Base.:==","text":"==(dom1::Domain, dom2::Domain)::Bool\n\nEquality for domains. Domains are equal if the have the same number of dimensions (see ndims), and the same lower and upper bounds (see first and last). The element types might be different if they can be compared for equality, e.g. Float64 and BigFloat.\n\nIncompatible domains (living in different number of dimensions) are treated as unequal.\n\n\n\n\n\n","category":"method"},{"location":"#Base.last-Tuple{Domain}","page":"GraviPet.jl","title":"Base.last","text":"last(dom::Domain{T})::T\n\nReturn the upper bound of the domain. (At the moment domains are assumed to be cuboid, i.e. to form a box or hyper-rectangle.) Domain bounds can be infinity.\n\nSee also first.\n\n\n\n\n\n","category":"method"},{"location":"#Base.isdisjoint-Tuple{Domain, Domain}","page":"GraviPet.jl","title":"Base.isdisjoint","text":"isdisjoint(dom1::Domain, dom2::Domain)::Bool\n\nTest whether dom1 is disjoint from dom2. Only compatible domains (which live in the same number of dimensions) can be compared.\n\n\n\n\n\n","category":"method"},{"location":"#Base.eltype-Tuple{Domain}","page":"GraviPet.jl","title":"Base.eltype","text":"eltype(dom::Domain)::Type\n\nReturn the element type of a Domain.\n\n\n\n\n\n","category":"method"},{"location":"#Base.ndims-Tuple{Domain}","page":"GraviPet.jl","title":"Base.ndims","text":"ndims(dom::Domain)::Int\n\nReturn the number of dimensions of a Domain. Scalar domains are one-dimensional. Zero-dimensional domains consist of a single point only.\n\n\n\n\n\n","category":"method"},{"location":"#Base.issubset-Tuple{Domain, Domain}","page":"GraviPet.jl","title":"Base.issubset","text":"issubset(dom1::Domain, dom2::Domain)::Bool\n\nTest whether dom1 is a subset of dom2. Only compatible domains (which live in the same number of dimensions) can be compared.\n\n\n\n\n\n","category":"method"},{"location":"#Base.in-Tuple{Any, Domain}","page":"GraviPet.jl","title":"Base.in","text":"in(elem, dom::Domain)::Bool\n\nTest whether elem is inside domain dom.\n\n\n\n\n\n","category":"method"},{"location":"#Base.first-Tuple{Domain}","page":"GraviPet.jl","title":"Base.first","text":"first(dom::Domain{T})::T\n\nReturn the lower bound of the domain. (At the moment domains are assumed to be cuboid, i.e. to form a box or hyper-rectangle.) Domain bounds can be infinity.\n\nSee also last.\n\n\n\n\n\n","category":"method"},{"location":"#Public-Helpers","page":"GraviPet.jl","title":"Public Helpers","text":"","category":"section"},{"location":"","page":"GraviPet.jl","title":"GraviPet.jl","text":"Modules = [GraviPet.Funs]","category":"page"},{"location":"#GraviPet.Funs.Fun","page":"GraviPet.jl","title":"GraviPet.Funs.Fun","text":"struct Fun{Tuple{Ts...},R}\n\nA function that takes inputs with the types Ts and returns a result of type R.\n\n\n\n\n\n","category":"type"}]
}

function test_Domain(::Val{D}, ::Type{T}, ::Type{DomainT}, makeDomainT) where {D,T,DomainT}
    @assert DomainT <: Domain
    dom, values, notvalues = makeDomainT()::Tuple{DomainT,AbstractVector{T},AbstractVector{T}}

    @test dom isa Domain{T}
    @test eltype(dom) == T
    @test ndims(dom) == D

    @test dom == dom

    for x in values
        @test x in dom
    end
    for x in notvalues
        @test !(x in dom)
    end

    # TODO: test issubset, isdisjoint
end

Random.seed!(0)
@testset "Interval T=$T" for T in [Float32, Float64, Double64, BigRat]
    function make_dom()
        xmin = rand((-T(1)):T(0.01):(+T(1)))
        xmax = rand((xmin + 1):T(0.01):(xmin + 10))
        dom = Interval(xmin, xmax)
        values = rand(xmin:T(0.01):xmax, 10)
        notvalues = filter(x -> !(xmin <= x <= xmax), rand(T(xmin - 1):T(0.01):T(xmax + 2), 10))
        return dom, values, notvalues
    end
    for iter in 1:10
        test_Domain(Val(1), T, Interval{T}, make_dom)
    end
end

Random.seed!(0)
@testset "Box D=$D T=$T" for D in 0:4, T in [Float32, Float64, Double64, BigRat]
    function make_dom()
        xmin = random_SVector(
            SVector{D,T}(T(-1) for d in 1:D), SVector{D,T}(T(0.01) for d in 1:D), SVector{D,T}(T(+1) for d in 1:D)
        )
        xmax = random_SVector(xmin .+ 1, SVector{D,T}(T(0.01) for d in 1:D), xmin .+ 10)
        dom = Box(xmin, xmax)
        values = [random_SVector(xmin, SVector{D,T}(T(0.01) for d in 1:D), xmax) for n in 1:10]
        notvalues = filter(
            x -> any((x .< xmin) .| (x .> xmax)),
            [random_SVector(xmin .- 1, SVector{D,T}(T(0.01) for d in 1:D), xmax .+ 2) for n in 1:100],
        )
        return dom, values, notvalues
    end
    for iter in 1:10
        test_Domain(Val(D), SVector{D,T}, Box{D,T}, make_dom)
    end
end

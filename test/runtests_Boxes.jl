Random.seed!(0)
@testset "Box D=$D T=$T" for D in 0:4, T in [Float32, Float64, Double64, BigRat]
    function make_dom()
        xmin = random_SVector(
            SVector{D,T}(T(-1) for d in 1:D), SVector{D,T}(T(1//100) for d in 1:D), SVector{D,T}(T(+1) for d in 1:D)
        )
        xmax = random_SVector(xmin .+ 1, SVector{D,T}(T(1//100) for d in 1:D), xmin .+ 10)
        dom = Box(xmin, xmax)
        values = [random_SVector(xmin, SVector{D,T}(T(1//100) for d in 1:D), xmax) for n in 1:10]
        notvalues = filter(
            x -> any((x .< xmin) .| (x .> xmax)),
            [random_SVector(xmin .- 1, SVector{D,T}(T(1//100) for d in 1:D), xmax .+ 2) for n in 1:100],
        )
        return dom, values, notvalues
    end
    for iter in 1:10
        test_Domain(Val(D), SVector{D,T}, Box{D,T}, make_dom)
    end
end

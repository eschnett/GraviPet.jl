Random.seed!(0)
@testset "Interval T=$T" for T in [Float32, Float64, Double64, BigRat]
    function make_dom()
        xmin = rand((-T(1)):T(1//100):(+T(1)))
        xmax = rand((xmin + 1):T(1//100):(xmin + 10))
        dom = Interval(xmin, xmax)
        values = rand(xmin:T(1//100):xmax, 10)
        notvalues = filter(x -> !(xmin ≤ x ≤ xmax), rand(T(xmin - 1):T(1//100):T(xmax + 2), 10))
        return dom, values, notvalues
    end
    for iter in 1:10
        test_Domain(Val(1), T, Interval{T}, make_dom)
    end
end

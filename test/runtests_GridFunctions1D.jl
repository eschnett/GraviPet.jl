Random.seed!(0)
@testset "GridFunction1D S=$T T=$T" for S in [Float32, Float64, Double64, BigRat], T in [S]
    function make_dom()
        xmin = rand((-S(1)):S(1//100):(+S(1)))
        xmax = rand((xmin + 1):S(1//100):(xmin + 10))
        dom = Interval(xmin, xmax)
        values = rand(xmin:S(1//100):xmax, 10)
        return dom, values
    end
    function make_gf(dom, cod)
        npoints = rand(2:20)
        grid = rand(first(cod):T(1//100):last(cod), npoints)
        gf = GridFunction1D("gf1d", dom, cod, grid)
        return gf
    end
    for iter in 1:10
        test_Category(Interval{T}, Interval{T}, GridFunction1D{T,T}, make_dom, make_dom, make_dom, make_gf)
    end
end

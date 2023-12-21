Random.seed!(0)
@testset "ThreadedFunction{GridFunction} S=$S^$DS T=$T^$DT" for DS in 0:4,
    S in [Float32, Float64, Double64, BigRat],
    DT in 0:4,
    T in [S]

    function make_dom()
        xmin = random_SVector(
            SVector{DS,S}(S(-1) for d in 1:DS), SVector{DS,S}(S(1//100) for d in 1:DS), SVector{DS,S}(S(+1) for d in 1:DS)
        )
        xmax = random_SVector(xmin .+ 1, SVector{DS,S}(S(1//100) for d in 1:DS), xmin .+ 10)
        dom = Box(xmin, xmax)
        values = [random_SVector(xmin, SVector{DS,S}(S(1//100) for d in 1:DS), xmax) for n in 1:10]
        return dom, values
    end
    function make_cod()
        ymin = random_SVector(
            SVector{DT,T}(T(-1) for d in 1:DT), SVector{DT,T}(T(1//100) for d in 1:DT), SVector{DT,T}(T(+1) for d in 1:DT)
        )
        ymax = random_SVector(ymin .+ 1, SVector{DT,T}(T(1//100) for d in 1:DT), ymin .+ 10)
        cod = Box(ymin, ymax)
        values = [random_SVector(ymin, SVector{DT,T}(T(1//100) for d in 1:DT), ymax) for n in 1:10]
        return cod, values
    end
    function make_gf(dom, cod)
        npoints = rand(2:10, DS)
        grid = [
            random_SVector(first(cod), SVector{DT,T}(1//100 for d in 1:DT), last(cod)) for i in CartesianIndices(Tuple(npoints))
        ]
        grid::AbstractArray{SVector{DT,T},DS}
        gf = GridFunction("gf", dom, cod, grid)
        return gf
    end
    function make_tf(dom, cod)
        make_gf′(name, dom, cod) = make_gf(dom, cod)
        tf = ThreadedFunction(make_gf′, "tf", dom, cod)
        return tf
    end
    for iter in 1:10
        test_Category(Box{DS,S}, Box{DT,T}, ThreadedFunction{DS,S,DT,T}, make_dom, make_cod, nothing, make_tf)
    end
end

Random.seed!(0)
@testset "BlockFunction{ThreadedFunction{GridFunction}} S=$S^$DS T=$T^$DT" for DS in 0:4,
    S in [Float32, Float64, Double64, BigRat],
    DT in 0:4,
    T in [S]

    function make_dom()
        xmin = random_SVector(
            SVector{DS,S}(S(-1) for d in 1:DS), SVector{DS,S}(S(1//100) for d in 1:DS), SVector{DS,S}(S(+1) for d in 1:DS)
        )
        xmax = random_SVector(xmin .+ 1, SVector{DS,S}(S(1//100) for d in 1:DS), xmin .+ 10)
        dom = Box(xmin, xmax)
        values = [random_SVector(xmin, SVector{DS,S}(S(1//100) for d in 1:DS), xmax) for n in 1:10]
        return dom, values
    end
    function make_cod()
        ymin = random_SVector(
            SVector{DT,T}(T(-1) for d in 1:DT), SVector{DT,T}(T(1//100) for d in 1:DT), SVector{DT,T}(T(+1) for d in 1:DT)
        )
        ymax = random_SVector(ymin .+ 1, SVector{DT,T}(T(1//100) for d in 1:DT), ymin .+ 10)
        cod = Box(ymin, ymax)
        values = [random_SVector(ymin, SVector{DT,T}(T(1//100) for d in 1:DT), ymax) for n in 1:10]
        return cod, values
    end
    function make_gf(dom, cod)
        npoints = rand(2:5, DS)
        grid = [
            random_SVector(first(cod), SVector{DT,T}(1//100 for d in 1:DT), last(cod)) for i in CartesianIndices(Tuple(npoints))
        ]
        grid::AbstractArray{SVector{DT,T},DS}
        gf = GridFunction("gf", dom, cod, grid)
        return gf
    end
    function make_tf(dom, cod)
        make_gf′(name, dom, cod) = make_gf(dom, cod)
        tf = ThreadedFunction(make_gf′, "tf", dom, cod)
        return tf
    end
    function make_bf(dom, cod)
        nblocks = rand(1:5, DS)
        make_tf′(name, dom, cod) = make_tf(dom, cod)
        bf = BlockFunction(make_tf′, "bf", dom, cod, SVector{DS,Int}(nblocks))
        return bf
    end
    for iter in 1:10
        test_Category(Box{DS,S}, Box{DT,T}, BlockFunction{DS,S,DT,T}, make_dom, make_cod, nothing, make_bf)
    end
end

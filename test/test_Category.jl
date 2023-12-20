function test_Category(
    ::Type{Dom}, ::Type{Cod}, ::Type{CategoryDomCod}, makeDom, makeCod, makeDomCod, makeCategoryDomCod
) where {Dom,Cod,CategoryDomCod}
    @assert Dom <: Domain
    @assert Cod <: Domain
    @assert CategoryDomCod <: Category

    dom, xvalues = makeDom()
    cod, yvalues = makeCod()
    cat = makeCategoryDomCod(dom, cod)::CategoryDomCod
    dom::Dom
    cod::Cod
    xvalues::AbstractVector{eltype(dom)}
    yvalues::AbstractVector{eltype(cod)}

    @test cat isa Category{Dom,Cod}
    @test domain(cat) isa Dom
    @test codomain(cat) isa Cod

    VS = eltype(dom)
    DS = ndims(dom)
    S = eltype(VS)
    atolS = S <: Rational ? zero(S) : 100 * eps(one(S))
    dx = VS <: AbstractArray ? VS(atolS for d in 1:DS) : atolS

    VT = eltype(cod)
    DT = ndims(cod)
    T = eltype(VT)
    atolT = T <: Rational ? zero(T) : 100 * eps(real(one(T)))
    dy = VT <: AbstractArray ? VT(atolT for d in 1:DT) : atolT
    for x in xvalues
        @test x in expanded(dom, dx)
    end

    @test cat == cat

    len = length(cat)
    @test isempty(cat) == (len == 0)
    if !isempty(cat)
        for n in 1:10
            i = rand(1:len)
            y = cat[i]
            @test y ∈ expanded(cod, dy)
        end
    end

    for x in xvalues
        y = cat(x)
        @test y isa eltype(cod)
        @test y ∈ expanded(cod, dy)
    end

    id = make_identity(cat)
    @test domain(id) == domain(cat)
    @test codomain(id) == domain(id)
    for x in xvalues
        y = id(x)
        @test y isa eltype(dom)
        @test y ≈ x atol = atolT
    end

    # TODO: test collection better
    # TODO: test function composition
    # TODO: test vector spaces
    # TODO: test `project`
end

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

Random.seed!(0)
@testset "GridFunction S=$S^$DS T=$T^$DT" for DS in 0:4, S in [Float32, Float64, Double64, BigRat], DT in 0:4, T in [S]
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
    for iter in 1:10
        test_Category(Box{DS,S}, Box{DT,T}, GridFunction{DS,S,DT,T}, make_dom, make_cod, nothing, make_gf)
    end
end

Random.seed!(0)
@testset "BlockFunction{GridFunction} S=$S^$DS T=$T^$DT" for DS in 0:4,
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
    function make_bf(dom, cod)
        nblocks = rand(1:5, DS)
        make_gf′(name, dom, cod) = make_gf(dom, cod)
        bf = BlockFunction(make_gf′, "bf", dom, cod, SVector{DS,Int}(nblocks))
        return bf
    end
    for iter in 1:10
        test_Category(Box{DS,S}, Box{DT,T}, BlockFunction{DS,S,DT,T}, make_dom, make_cod, nothing, make_bf)
    end
end

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

Random.seed!(0)
@testset "DistributedFunction{GridFunction} S=$S^$DS T=$T^$DT" for DS in 0:4,
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
    function make_df(dom, cod)
        make_gf′(name, dom, cod) = make_gf(dom, cod)
        df = DistributedFunction(make_gf′, "df", dom, cod)
        return df
    end
    for iter in 1:10
        test_Category(Box{DS,S}, Box{DT,T}, DistributedFunction{DS,S,DT,T}, make_dom, make_cod, nothing, make_df)
    end
end

Random.seed!(0)
@testset "BlockFunction{DistributedFunction{GridFunction}} S=$S^$DS T=$T^$DT" for DS in 0:4,
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
    function make_df(dom, cod)
        make_gf′(name, dom, cod) = make_gf(dom, cod)
        df = DistributedFunction(make_gf′, "df", dom, cod)
        return df
    end
    function make_bf(dom, cod)
        nblocks = rand(1:5, DS)
        make_df′(name, dom, cod) = make_df(dom, cod)
        bf = BlockFunction(make_df′, "bf", dom, cod, SVector{DS,Int}(nblocks))
        return bf
    end
    for iter in 1:10
        test_Category(Box{DS,S}, Box{DT,T}, BlockFunction{DS,S,DT,T}, make_dom, make_cod, nothing, make_bf)
    end
end

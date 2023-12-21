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

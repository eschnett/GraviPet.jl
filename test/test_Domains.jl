function test_Domain(::Val{D}, ::Type{T}, ::Type{DomainT}, makeDomainT) where {D,T,DomainT}
    @assert DomainT <: Domain
    dom, values, notvalues = makeDomainT()::Tuple{DomainT,AbstractVector{T},AbstractVector{T}}

    @test dom isa Domain{T}
    @test eltype(dom) == T
    @test ndims(dom) == D

    @test dom == dom

    @test dom ⊆ dom
    @test !isdisjoint(dom, dom)

    for x in values
        @test x ∈ dom
    end
    for x in notvalues
        @test x ∉ dom
    end

    # TODO: test expanded, shifted
    # TODO: test issubset, isdisjoint better
end

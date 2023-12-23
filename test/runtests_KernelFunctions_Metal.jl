@static if VERSION >= v"1.8"
    # Metal requires at least Julia 1.8
    if Metal.functional()
        Random.seed!(0)
        # `Float64` is not supported on Metal, `Double32` does not work
        @testset "KernelFunction [Metal] S=$S^$DS T=$T^$DT" for DS in 0:4, S in [Float16, Float32], DT in 0:4, T in [S]
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
            function make_kf(dom, cod)
                npoints = rand(2:10, DS)
                grid = [
                    random_SVector(first(cod), SVector{DT,T}(1//100 for d in 1:DT), last(cod)) for
                    i in CartesianIndices(Tuple(npoints))
                ]
                grid::AbstractArray{SVector{DT,T},DS}
                grid = MtlArray(grid)
                kf = KernelFunction("kf", dom, cod, grid)
                return kf
            end
            for iter in 1:10
                test_Category(Box{DS,S}, Box{DT,T}, KernelFunction{DS,S,DT,T}, make_dom, make_cod, nothing, make_kf)
            end
        end
    end
end

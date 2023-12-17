# Need to use Julia 1.9

@static if (VERSION.major, VERSION.minor) != (1, 9)
    @warn "*** This code might require using Julia 1.9 ***"
end

using KernelAbstractions
using LinearAlgebra
# using Random
using StaticArrays

@kernel function init_kernel!(u::AbstractArray{T,D}, @Const(xmin::SVector{D,S}), @Const(dx::SVector{D,S})) where {D,S,T}
    i0 = @index(Global, NTuple)
    i = SVector{D}(i0)
    # ndr = @ndrange()
    # grs = @groupsize()
    # gli = @index(Global, NTuple)
    # gri = @index(Group, NTuple)
    # loi = @index(Local, NTuple)
    # if i == j == 100
    #     @print("ndrange:      ", ndr, "\n")
    #     @print("groupsize:    ", grs, "\n")
    #     @print("global index: ", gli, "\n")
    #     @print("group index:  ", gri, "\n")
    #     @print("local index:  ", loi, "\n")
    # end
    x = xmin + (i .- 1) .* dx
    ci = CartesianIndex(Tuple(i))
    u[ci] = prod(sinpi(a) for a in x)
    if false
    end
end

function init!(u::AbstractArray{T,D}, xmin::SVector{D,S}, dx::SVector{D,S}) where {D,S,T}
    D::Integer
    backend = KernelAbstractions.get_backend(u)
    kernel! = init_kernel!(backend)
    kernel!(u, xmin, dx; ndrange=size(u))
    return nothing
end

@kernel function diffx_kernel!(
    du::AbstractArray{T,D}, @Const(u::AbstractArray{T,D}), @Const(xmin::SVector{D,S}), @Const(dx::SVector{D,S})
) where {D,S,T}
    i0 = @index(Global, NTuple)
    i = SVector{D}(i0)
    imin = SVector{D}(ntuple(d -> 1, D))
    imax = SVector{D}(@ndrange())

    dir = 1

    di = SVector{D}(ntuple(d -> d == dir, D))
    ci = CartesianIndex(Tuple(i))
    cim = CartesianIndex(Tuple(i - di))
    cip = CartesianIndex(Tuple(i + di))

    if i[dir] == imin[dir]
        du[ci] = (u[cip] - u[ci]) / dx[dir]
    elseif i[dir] == imax[dir]
        du[ci] = (u[ci] - u[cim]) / dx[dir]
    else
        du[ci] = (u[cip] - u[cim]) / 2dx[dir]
    end
end

function diffx!(du::AbstractArray{T,D}, u::AbstractArray{T,D}, xmin::SVector{D,S}, dx::SVector{D,S}) where {D,S,T}
    D::Integer
    @assert axes(du) == axes(u)
    backend = KernelAbstractions.get_backend(du)
    kernel! = diffx_kernel!(backend)
    kernel!(du, u, xmin, dx; ndrange=size(du))
    return nothing
end

@kernel function sum2_kernel!(s::AbstractArray{T,0}, @Const(u::AbstractArray{T,D})) where {D,T}
    i = @index(Global, Linear)
    KernelAbstractions.@atomic :acquire_release s[] += u[i]^2
end

function sum2!(s::AbstractArray{T,0}, u::AbstractArray{T,D}) where {D,T}
    D::Integer
    backend = KernelAbstractions.get_backend(u)
    kernel! = sum2_kernel!(backend)
    kernel!(s, u; ndrange=size(u))
    return nothing
end

# using CUDA
# using CUDA.CUDAKernels
# const backend = CUDABackend()
# CUDA.allowscalar(false)

using Metal

function main()
    backend = CPU()

    # Metal.versioninfo()
    # display(Metal.devices())
    # Metal.device!(Metal.devices()[2])
    # backend = MetalBackend()

    D = 3
    T = Float32
    # n = SVector{D}(256, 256, 256)
    n = SVector{D}(512, 512, 512)
    # n = SVector{D}(1024, 1024, 1024)
    xmin = T(0) * n
    dx = T(1) ./ (n .- 1)

    u = KernelAbstractions.zeros(backend, T, n...)
    init!(u, xmin, dx)

    du = KernelAbstractions.zeros(backend, T, n...)
    diffx!(du, u, xmin, dx)

    #TODO s = KernelAbstractions.zeros(backend, T)
    #TODO sum2!(s, u)

    KernelAbstractions.synchronize(backend)

    u = Array(u)
    s = zeros(T)
    s[] = sum(u .^ 2)

    s[] = sqrt(s[] / length(u))

    @show s[]

    return nothing
end

main()

nothing

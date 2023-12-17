module ExampleInitialize

using GraviPet
using LinearAlgebra
using StaticArrays

################################################################################

renamed(name::AbstractString, gf::GridFunction) = typeof(gf)(name, gf.domain, gf.codomain, gf.grid)

################################################################################

# # Set up domain and codomain
# xmin = SVector(0.0, 0.0)
# xmax = SVector(1.0, 1.0)
# dom = Box(xmin, xmax)
# yinf = SVector(Inf)
# cod = Box(-yinf, +yinf)
# 
# # Define a Julia function
# f(x, y) = sinpi(x) * cospi(y)
# jf = JuliaFunction("wave", dom, cod, xs -> SVector(f(xs...)))
# 
# I0 = nothing
# for res in 0:0   #TODO 0:4
# 
#     # Choose a discretization
#     npoints = 10 * 2^res + 1
#     print("npoints: $npoints... ")
#     gf0 = GridFunction("skeleton", dom, dom, SVector(npoints, npoints))
# 
#     # Project the Julia function onto the discretization
#     gf = project(gf0, jf)
# 
#     # Calculate the discretization error
#     jf1 = project(jf, gf)
#     # @show integrate(jf)
#     @show integrate(gf)
#     @show integrate(jf1)
# 
#     THIS DOES NOT WORK BECAUSE HCUBATURE CANNOT HANDLE THE KINKS
# 
#     Δjf = jf1 - jf
#     I = integrate(Δjf)
#     print("   ‖error‖₂: $(norm(I))")
# 
#     if I0 !== nothing
#         print("   factor: $(norm(I0) / norm(I))")
#     end
#     global I0 = I
# 
#     println()
# end

################################################################################

function main()

    # Set up domain and codomain
    xmin = SVector(0.0)
    xmax = SVector(1.0)
    dom = Box(xmin, xmax)
    yinf = SVector(Inf)
    cod = Box(-yinf, +yinf)

    # Define a Julia function
    f(x) = 1.0 # sinpi(x)
    jf = JuliaFunction("wave", dom, cod, xs -> SVector(f(xs...)))

    npoints = 11
    gf0 = GridFunction("skeleton", dom, dom, SVector(npoints))

    # Project the Julia function onto the discretization
    gf = project(gf0, jf)
    @show map(xs -> xs[1], gf.grid)

    coords = renamed("coordinates", make_identity(gf0))
    gf1 = map(xs -> SVector(f(xs...)), coords)
    @show map(xs -> xs[1], gf1.grid)

    return nothing
end

end

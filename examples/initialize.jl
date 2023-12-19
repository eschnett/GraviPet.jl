module ExampleInitialize

using CairoMakie
using Dates
using GraviPet
using LinearAlgebra
using MappedArrays
using SixelTerm
using StaticArrays

# Both GraviPet and CairoMakie export `Box`
const Box = GraviPet.Box

################################################################################

renamed(name::AbstractString, gf::GridFunction) = typeof(gf)(name, gf.domain, gf.codomain, gf.grid)

function provenance()
    user = get(ENV, "USER", "anonymous")
    host = gethostname()
    date = today()
    # return "GraviPet\n$user@$host\n$date"
    return "GraviPet\n$user\n$date"
end

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

    ########################################

    # Define a Julia function
    # This function takes an `SVector` as input and produces an `SVector` as output.
    f(x) = SVector(sinpi(4 * x[1]))
    jf = JuliaFunction("wave", dom, cod, f)

    npoints = 21
    gf0 = GridFunction("skeleton", dom, dom, SVector(npoints))

    # Project the Julia function onto the discretization
    gf = project(gf0, jf)

    coords = renamed("coordinates", make_identity(gf0))
    gf1 = map(xs -> SVector(f(xs...)), coords)

    ########################################

    # Prepare data
    xs = mappedarray(xy -> xy[1], coords.grid)
    data = mappedarray(val -> val[1], gf.grid)
    data1 = mappedarray(val -> val[1], gf1.grid)

    # Visualize it!
    fig = Figure(; fontsize=30, size=(1280, 960))
    ax = Axis(fig[1, 1]; title=gf.name)
    obj1 = scatterlines!(xs, data1; color=:blue, marker=:circle)
    obj = scatterlines!(xs, data; color=:red, marker=:rect)

    Label(fig[1, 2][1, 1], provenance(); justification=:left, padding=(10, 10, 10, 10))
    CairoMakie.Box(fig[1, 2][1, 1]; color=(:black, 0.15), strokewidth=0)
    Legend(fig[1, 2][2, 1], [obj1, obj], ["sampled", "projected"])

    # colsize!(fig.layout, 1, Aspect(1, 1.0))
    rowsize!(fig.layout, 1, Aspect(1, 0.5))

    display(fig)

    ########################################

    return nothing
end

end

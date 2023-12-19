module ExampleVisualize

using CairoMakie
using Dates
using GraviPet
using MappedArrays
using SixelTerm
using StaticArrays

# Both GraviPet and CairoMakie export `Box`
const Box = GraviPet.Box

################################################################################

renamed(name::AbstractString, gf::GridFunction) = GridFunction(name, gf.domain, gf.codomain, gf.grid)

function provenance()
    user = get(ENV, "USER", "anonymous")
    host = gethostname()
    date = today()
    # return "GraviPet\n$user@$host\n$date"
    return "GraviPet\n$user\n$date"
end

################################################################################

function main()

    # Set up domain
    xmin = SVector(0.0, 0.0)
    xmax = SVector(1.0, 1.0)
    dom = Box(xmin, xmax)

    # Choose a discretization
    gf0 = GridFunction("skeleton", dom, dom, SVector(101, 101))
    coords = renamed("coordinates", make_identity(gf0))

    # Define a function
    # This function takes an `SVector` as input and produces an `SVector` as output.
    f(x) = SVector(sinpi(4 * hypot(x[1], x[2])))

    # Set up a grid function
    gf = renamed("wave", map(f, coords))

    # Prepare data
    xs = mappedarray(xy -> xy[1], view(coords.grid, :, 1))
    ys = mappedarray(xy -> xy[2], view(coords.grid, 1, :))
    data = mappedarray(val -> val[1], gf.grid)

    # Visualize it!
    fig = Figure(; fontsize=30, size=(1280, 960))
    ax = Axis(fig[1, 1]; title=gf.name)
    obj = contourf!(xs, ys, data; colormap=:plasma)

    Label(fig[1, 2][1, 1], provenance(); justification=:left, padding=(10, 10, 10, 10))
    CairoMakie.Box(fig[1, 2][1, 1]; color=(:black, 0.15), strokewidth=0)
    cb = Colorbar(fig[1, 2][2, 1], obj; label="amplitude", width=20)
    # Remove superfluous white space on the right of the colorbar
    cb.alignmode = Mixed(; right=0)

    # Set aspect ratio
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    # rowsize!(fig.layout, 1, Aspect(1, 1.0))

    display(fig)

    return nothing
end

end

using Plots

gc(y, c, α) = c/√(1-y^α)
dgc(y, c, α) = (c*α*y^(α-1))/(2√(1-y^α)*(1-y^α))

"""
  Connects the x->x^2 and the gc part of the g function so that g have C^2 regularity
"""
function g_spline_interpolation(c, α)
    y1, y2 = 1/2, 3/4; z1, z2 = 1/4, gc(3/4, c, α)
    # Compute independant part of the spline
    t(y) = (y-y1)/(y2-y1);
    a = 2*y1*(y2-y1) - (z2-z1)
    b = -dgc(y2,c,α)*(y2-y1) + (z2-z1)
    # Assemble linear part + derivative part
    y -> (1-t(y))*z1 + t(y)*z2 + t(y)*(1-t(y))*( (1-t(y))*a + t(y)*b )
end

function g_cα(y, c, α)
    # Define 3 parts of the function
    g1 = y->y'y; g2 = g_spline_interpolation(c,α); g3 = y->gc(y,c,α)
    (0 ≤ y < 1/2)   && (return g1(y))
    (1/2 ≤ y < 3/4) && (return g2(y))
    (3/4 ≤ y < 1)   && (return g3(y))
    error("The regularization g function is defined on [0,1). Did you devide by √Ecut ?")
    nothing
end

function plot_g_cα(c, α)
    x_axis = LinRange(0, 0.99, 100)
    p = plot(x_axis, g_cα.(x_axis, c, α), label=:none, xlims=[0,1])
    vline!([1/2,3/4], label="C^2 spline interpolation", linestyle=:dash)
    vline!([1], label=:none, linestyle=:dash)
    plot!(x->x^2, label="x->x^2", linestyle=:dash, color=:black)
    plot!(legendfontsize=12, legend=:topleft)
    plot!(size=(750,500))
    p
end

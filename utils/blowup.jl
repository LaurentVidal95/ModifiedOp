using Optim

"""
New blowup structure that encodes the wanted regularity.
Only used for this paper.
"""
struct VariableBlowupCHV{F}
    p::Real             # Targeted regularity
    blowup_function::F
end

# Contruct a blowup given targeted regularity and interpolation interval
VariableBlowupCHV(p; interval=[0.85, 0.90]) =
    VariableBlowupCHV(p, blowup_function(p; interval))
(blowup::VariableBlowupCHV)(x) = blowup.blowup_function(x)

"""
Smooth interpolation but which in turn can give results slightly under the x-> x^2 curve
Negligeable for a choice of blow-up under 3.
"""
function smooth_interpolation(g3, interval)
    a, b = interval
    @assert b>a
    f(x) = (x==0) ? 0 : exp(-1/x)
    step(x) = f((x-a)/(b-a)) / (f((x-a)/(b-a)) + f(1-(x-a)/(b-a)))
    x -> (1-step(x))*x^2 + step(x)*g3(x)
end

"""
Compute the blowup_function at x giving the 3 parts of the function
"""
function blowup_function(x, g1, g2, g3, interval)
    x1, x2 = interval
    # Define 3 parts of the function
    (0 ≤ x < x1)   && (return g1(x))
    (x1 ≤ x < x2)  && (return g2(x))
    (x2 ≤ x < 1)   && (return g3(x))
    error("The blow-up function is defined on [0,1). Did you devide by √Ecut ?")
    nothing
end
function blowup_function(g3, interval)
    g2 = smooth_interpolation(g3, interval)
    g1 = x->x'x
    x->blowup_function(x, g1, g2, g3, interval)
end
function blowup_function(blowup_rate::T; interval) where {T<:Real}
    g3 = optimal_g3(blowup_rate, interval)
    blowup_function(g3, interval)
end

"""
Blow up part of the blow-up function. Good choice is a=0.48 and ε=-1.
ε gives the blow up rate of gm.
"""
Ca(a, ε) = (3/2)*(a^2)*(1-a)^(ε)
ha(a, ε) = y-> Ca(a, ε)/( (1-y)^(ε) )
function optimal_g3(ε, interval)
    x_axis = LinRange(0,0.99, 1000)
    function f(a)
        !(0≤a[1]<1) && (return 1e6)
        g3 = ha(only(a), ε)
        Gm = blowup_function(g3, interval)
        out = Gm.(x_axis) - x_axis .^2
        # Penalize Gm that are bellow the x->x^2 curve
        norm(out,1) .+ norm(out[out .< 0],1)*1e4
    end
    a_opti = only(optimize(f, [0.3], LBFGS()).minimizer)
    ha(a_opti, ε)
end

# """
# Only C2 pol interpolation.
# -> Above x^2 garantied but coefficients can blow-up with the blow-up rate.
# """
# function C2_pol_interpolation(g3, interval)
#     # ForwardDiff second derivative of g3
#     dg3(x) = ForwardDiff.derivative(g3, x)
#     d2g3(x) = ForwardDiff.derivative(dg3, x)

#     # Points of interpolation
#     a, b = interval
#     x1, x2 = a^2, g3(b)
#     y1, y2 = 2*a, dg3(b)
#     z1, z2 = 2, d2g3(b)

#     # Solve interpolation linear system
#     A = [1 a a^2  a^3   a^4    a^5;
#          1 b b^2  b^3   b^4    b^5;
#          0 1 2*a  3*a^2 4*a^3  5*a^4;
#          0 1 2*b  3*b^2 4*b^3  5*b^4;
#          0 0 2    6*a   12*a^2 20*a^3;
#          0 0 2    6*b   12*b^2 20*b^3]
#     B = [x1, x2, y1, y2, z1, z2]                      
#     C = A\B

#     # Assemble polynomial
#     x -> C'*((x * ones(6)) .^(0:5))
# end

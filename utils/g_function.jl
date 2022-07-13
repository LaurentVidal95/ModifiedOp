function C2_pol_interpolation(ha; interval=[0.5, 0.75])
    # ForwardDiff second derivative of ha
    dha(x) = ForwardDiff.derivative(ha, x)
    d2ha(x) = ForwardDiff.derivative(dha, x)

    # Points of interpolation
    a, b = interval
    x1, x2 = a^2, ha(b)
    y1, y2 = 2*a, dha(b)
    z1, z2 = 2, d2ha(b)

    # Solve interpolation linear system
    A = [1 a a^2  a^3   a^4    a^5;
         1 b b^2  b^3   b^4    b^5;
         0 1 2*a  3*a^2 4*a^3  5*a^4;
         0 1 2*b  3*b^2 4*b^3  5*b^4;
         0 0 2    6*a   12*a^2 20*a^3;
         0 0 2    6*b   12*b^2 20*b^3]
    B = [x1, x2, y1, y2, z1, z2]                      
    C = A\B

    # Assemble polynomial
    x -> C'*((x * ones(6)) .^(0:5))
end

function gm(x, g1, g2, g3; interval=[0.5, 0.75])
    x1, x2 = interval
    # Define 3 parts of the function
    (0 ≤ x < x1)   && (return g1(x))
    (x1 ≤ x < x2)  && (return g2(x))
    (x2 ≤ x < 1)   && (return g3(x))
    error("The regularization g function is defined on [0,1). Did you devide by √Ecut ?")
    nothing
end
function gm(g3; interval=[0.5, 0.75])
    g2 = C2_pol_interpolation(g3; interval)
    g1 = x->x'x
    x->gm(x, g1, g2, g3; interval)
end

"""
Some tests parameters:
For [0.85, 0.90]:
    C = 0.51
    (blow_up_rate==3//2) && (C=0.175)
    (blow_up_rate==5//2) && (C=0.06)
"""
function gm(blow_up_rate::T; interval=[0.5, 0.75]) where {T<:Real}
    C = 0.51
    (blow_up_rate==3//2) && (C=0.175)
    (blow_up_rate==5//2) && (C=0.06)
    g3 = ha(C, blow_up_rate)
    gm(g3; interval)
end

"""
Blow up part of gm. Good choice is a=0.48 and ε=-1.
ε gives the blow up rate of gm.
"""
Ca(a, ε) = (3/2)*(a^2)*(1-a)^(ε)
ha(a, ε) = y-> Ca(a, ε)/( (1-y)^(ε) )

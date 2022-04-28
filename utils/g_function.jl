using ForwardDiff

"""
  Connects the x->x^2 and the ha part of the g function so that g have C^2 regularity
"""
function gm_C2_spline_interpolation(ha)
    y1, y2 = 1/2, 3/4; z1, z2 = 1/4, ha(3/4)
    # Compute independant part of the spline
    t(y) = (y-y1)/(y2-y1);
    t1 = 2*y1*(y2-y1) - (z2-z1)
    t2 = -ForwardDiff.derivative(ha, y2)*(y2-y1) + (z2-z1)
    # Assemble linear part + derivative part
    y -> (1-t(y))*z1 + t(y)*z2 + t(y)*(1-t(y))*( (1-t(y))*t1 + t(y)*t2 )
end

function gm(y, ha)
    # Define 3 parts of the function
    g1 = y->y'y; g2 = gm_C2_spline_interpolation(ha); g3 = y->ha(y)
    (0 ≤ y < 1/2)   && (return g1(y))
    (1/2 ≤ y < 3/4) && (return g2(y))
    (3/4 ≤ y < 1)   && (return g3(y))
    error("The regularization g function is defined on [0,1). Did you devide by √Ecut ?")
    nothing
end

"""
Blow up part of gm. Good choice is a=0.48 and ε=-1.
ε gives the blow up rate of gm.
"""
Ca(a, ε) = (3/2)*(a^2)*(1-a)^(1/2 - ε)
ha(a, ε) = y-> Ca(a, ε)/( (1-y)^(1/2 - ε ))

"""
Old ha with blow up |⋅|^{-1/2}
"""
# Good choice is a = 1/(3.5) and ε = 1
# ha_1(y, a, ε) = a/√(1-y^ε)
# ha_1(y; a=1/(3.5), ε=1) = ha_1(y, a, ε)
# dha(y, a, ε) = (a*ε*y^(ε-1))/(2√(1-y^ε)*(1-y^ε))

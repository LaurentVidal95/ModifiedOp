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
(blowup::VariableBlowupCHV)(x, Ecut) = blowup.blowup_function(x, Ecut)

"""
Smooth interpolation but which in turn can give results slightly under the x-> x^2 curve
Negligeable for a choice of blow-up under 3.
"""
function smooth_interpolation(blowup_part, interval)
    a, b = interval
    @assert b>a
    f(x) = (x==0) ? 0 : exp(-1/x)
    step(x) = f((x-a)/(b-a)) / (f((x-a)/(b-a)) + f(1-(x-a)/(b-a)))
    x -> (1-step(x))*x^2 + step(x)*blowup_part(x)
end

"""
Compute the blowup_function at x giving the 3 parts of the function
"""
function blowup_function(y, Ecut, blowup_part, interval)
    # precomputations
    E_kin = y^2/2
    x = y / √(2Ecut)
    x1, x2 = interval
    # Define the 3 parts of the function
    (0 ≤ x < x1)   && (return 1)
    (x1 ≤ x < x2)  && (return (Ecut/E_kin)*smooth_interpolation(blowup_part, interval)(x))
    (x2 ≤ x < 1)   && (return (Ecut/E_kin)*blowup_part(x))
    Inf # Handle the case |k+G|^2 = 2Ecut 
end
function blowup_function(blowup_rate::T; interval) where {T<:Real}
    blowup_part = optimal_blowup_part(blowup_rate, interval)
    (x, Ecut) -> blowup_function(x, Ecut, blowup_part, interval)
end

"""
Blow up part of the blow-up function. Good choice is a=0.48 and ε=-1.
ε gives the blow up rate of gm.
"""
Ca(a, ε) = (3/2)*(a^2)*(1-a)^(ε)
ha(a, ε) = y-> Ca(a, ε)/( (1-y)^(ε) )
function optimal_blowup_part(ε, interval)
    x_axis = LinRange(0,0.99, 1000)
    function f(a)
        !(0≤a[1]<1) && (return 1e6)
        blowup_part = ha(only(a), ε)
        Gm(y) = blowup_function(y, 1/2, blowup_part, interval)
        out = (Gm.(x_axis) .- 1) .* ( x_axis .^2 )
        # Penalize Gm that are bellow the x->x^2 curve
        norm(out,1) .+ norm(out[out .< 0],1)*1e4
    end
    a_opti = only(optimize(f, [0.3], LBFGS()).minimizer)
    ha(a_opti, ε)
end

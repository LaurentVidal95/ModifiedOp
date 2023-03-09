using Optim

"""
Defines a blow-up function 𝒢 with custom blow-up parameter p. It is defined as such.
Let ``I = [x₁, x₂]⊂(1/2,1)`` be an given interval. Define ``𝒽(x) = C(a,p)(1-x)^{-p}``
with ``C(a,p) = (3/2)a²(1-a)ᵖ``. The function ``𝒢`` is such that:

  •  ``𝒢(x) = x²``  on ``[0,\frac{1}{2})``
  •  ``𝒢(x) = 𝒽(x)`` on ``(x₂,1)``
  • The value of 𝒢 on I is obtained by a the smooth interpolation ``itp(x) = (1-step(x))x² + step(x)𝒽(x)``
    with `` step(x) = f((x-x₁)/(x₂-x₁)) / (f((x-x₁)/(x₂-x₁)) + f(1-(x-x₁)/(x₂-x₁)))``, ``f(x) = e^{-1/x}``.

The parameter a, hence the constant ``C(a,p)`` is chosen so that the interpolation ``itp`` is
as close as possible to the x→x² curve on I.
"""

struct VariableBlowupCHV{F}
    p::Real             # Targeted regularity
    blowup_function::F
end

# Contruct a blowup given targeted regularity and interpolation interval
DefaultInterval=[0.85, 0.90]
VariableBlowupCHV(p; interval=DefaultInterval) =
    VariableBlowupCHV(p, 𝒢(p; interval))
(blowup::VariableBlowupCHV)(x, Ecut) = blowup.blowup_function(x, Ecut)

"""
Smooth interpolation but which in turn can give results slightly under the x-> x^2 curve
Negligeable for a choice of blow-up under 3.
"""
function smooth_interpolation(blowup_part, interval)
    x₁, x₂ = interval
    @assert x₂>x₁
    f(x) = (x==0) ? 0 : exp(-1/x)
    step(x) = f((x-x₁)/(x₂-x₁)) / (f((x-x₁)/(x₂-x₁)) + f(1-(x-x₁)/(x₂-x₁)))
    x -> (1-step(x))*x^2 + step(x)*blowup_part(x)
end

"""
Compute the blowup_function at x giving the 3 parts of the function
"""
function 𝒢(y, Ecut, 𝒽, interval)
    # precomputations
    E_kin = y^2/2
    x = y / √(2Ecut)
    x₁, x₂ = interval
    # Define the 3 parts of the function
    (0 ≤ x < x₁)   && (return 1)
    (x₁ ≤ x < x₂)  && (return (Ecut/E_kin)*smooth_interpolation(𝒽, interval)(x))
    (x₂ ≤ x < 1)   && (return (Ecut/E_kin)*𝒽(x))
    Inf # Handle the case |k+G|^2 = 2Ecut 
end
function 𝒢(blowup_rate::T; interval) where {T<:Real}
    𝒽 = optimal_𝒽(blowup_rate, interval)
    (x, Ecut) -> 𝒢(x, Ecut, 𝒽, interval)
end

"""
Blow up part 𝒽 of the blow-up function 𝒢, that ensures that 𝒢 is
as close as possible to the x→x² curve.
"""
Ca(a, p) = (3/2)*(a^2)*(1-a)^(p)
ha(a, p) = y-> Ca(a, p)/( (1-y)^(p) )
function optimal_𝒽(p, interval)
    x_axis = LinRange(0,0.99, 1000)
    function f(a)
        !(0≤a[1]<1) && (return 1e6)
        𝒽 = ha(only(a), p)
        G(y) = 𝒢(y, 1/2, 𝒽, interval)
        out = (G.(x_axis) .- 1) .* ( x_axis .^2 )
        # Penalize blow-up functions that fall bellow the x->x^2 curve
        norm(out,1) .+ norm(out[out .< 0],1)*1e4
    end
    a_opti = only(optimize(f, [0.3], LBFGS()).minimizer)
    ha(a_opti, p)
end

using DiffEqBase
# DDE Example Problems
# examples with constant delays
export prob_dde_1delay, prob_dde_1delay_notinplace, prob_dde_1delay_scalar_notinplace,
       prob_dde_2delays, prob_dde_2delays_notinplace, prob_dde_2delays_scalar_notinplace,
       prob_dde_1delay_long, prob_dde_1delay_long_notinplace,
       prob_dde_1delay_long_scalar_notinplace, prob_dde_2delays_long,
       prob_dde_2delays_long_notinplace, prob_dde_2delays_long_scalar_notinplace,
       prob_dde_mackey, prob_dde_wheldon, prob_dde_qs,
# examples with vanishing time dependent delays
       prob_ddde_neves1, prob_dde_neves_thompson,
# examples with state dependent delays
       prob_dde_paul1, prob_dde_paul2, prob_dde_mahaffy1, prob_dde_mahaffy2,
# examples with vanishing state dependent delays
       prob_neves2, prob_dde_gatica


# DDE examples with analytical solution

## Single constant delay

### In-place function

function f_1delay(du, u, h, p, t)
    du[1] = - h(p, t - oneunit(t))[1] / oneunit(t)
end

function f_1delay_analytic(u₀, p, t)
    z = t/oneunit(t)

    if z < 0
        return zero(u₀)
    elseif z < 1
        return u₀
    elseif z ≤ 10
        if z < 2
            c = @evalpoly(z, 2, -1)
        elseif z < 3
            c = @evalpoly(z, 4, -3, 1//2)
        elseif z < 4
            c = @evalpoly(z, 17//2, -15//2, 2, -1//6)
        elseif z < 5
            c = @evalpoly(z, 115//6, -109//6, 6, -5//6, 1//24)
        elseif z < 6
            c = @evalpoly(z, 1085//24, -1061//24, 197//12, -35//12, 1//4, -1//120)
        elseif z < 7
            c = @evalpoly(z, 13201//120, -13081//120, 521//12, -107//12, 1, -7//120, 1//720)
        elseif z < 8
            c = @evalpoly(z, 39371//144, -39227//144, 27227//240, -3685//144, 487//144,
                          -21//80, 1//90, -1//5040)
        elseif z < 9
            c = @evalpoly(z, 1158379//1680, -1156699//1680, 212753//720, -51193//720,
                          1511//144, -701//720, 1//18, -1//560, 1//40320)
        else
            c = @evalpoly(z, 23615939//13440, -23602499//13440, 7761511//10080,
                          -279533//1440, 89269//2880, -1873//576, 323//1440, -11//1120,
                          1//4032, -1//362880)
        end

        return c .* u₀
    else
        error("This analytical solution is only valid on (-∞,10]")
    end
end

ff_1delay = DDEFunction(f_1delay,analytic=f_1delay_analytic)

build_prob_dde_1delay(u₀, ::T=u₀) where {T} =
    DDEProblem(ff_1delay, [u₀], (p, t)->[zero(u₀)], (zero(T), T(10)),
               constant_lags = [oneunit(T)])

"""
    prob_dde_1delay

Model problem of finding a solution ``u(t)`` in the time span ``t \\in [0,10]`` to the delay
differential equation

```math
\\frac{du}{dt} = -u(t-1)
```

with history function

```math
u(t) = \\begin{cases}
0 & \\text{if } t < 0, \\
u₀ & \\text{if } t = 0,
\\end{cases}
```

for ``t \\leq 0`` with ``u₀ = 1``. Note that the problem is discontinuous at ``t = 0`` for
all ``u₀ \\neq 0``.

An analytical solution of this problem is provided for ``t \\in (-\\infty,10]``.
"""
prob_dde_1delay = build_prob_dde_1delay(1.0)

### Not in-place function

function f_1delay_notinplace(u, h, p, t)
    - h(p, t - oneunit(t)) ./ oneunit(t)
end

ff_1delay_notinplace = DDEFunction(f_1delay_notinplace,analytic=f_1delay_analytic)

#### Vectorized history function

build_prob_dde_1delay_notinplace(u₀, ::T=u₀) where {T} =
    DDEProblem(ff_1delay_notinplace, [u₀], (p, t)->[zero(u₀)], (zero(T), T(10)),
               constant_lags = [oneunit(T)])

"""
    prob_dde_1delay_notinplace

Same as [`prob_dde_1delay`](@ref), but purposefully implemented with a not in-place
function.
"""
prob_dde_1delay_notinplace = build_prob_dde_1delay_notinplace(1.0)

#### Scalar history function

build_prob_dde_1delay_scalar_notinplace(u₀, ::T=u₀) where {T} =
    DDEProblem(ff_1delay_notinplace, u₀, (p, t) -> zero(u₀), (zero(T), T(10)),
               constant_lags = [oneunit(T)])

"""
    prob_dde_1delay_scalar_notinplace

Same as [`prob_dde_1delay_notinplace`](@ref), but purposefully implemented with a scalar
history function.
"""
prob_dde_1delay_scalar_notinplace = build_prob_dde_1delay_scalar_notinplace(1.0)

## Two constant delays

### In-place function

function f_2delays(du, u, h, p, t::T) where T
    du[1] = (- h(p, t - T(1//3))[1] - h(p, t - T(1//5))[1]) / oneunit(t)
end

function f_2delays_analytic(u₀, p, t)
    z = t/oneunit(t)

    if z < 0
        return zero(u₀)
    elseif z < 1/5
        return u₀
    elseif z ≤ 1
        if z < 1/3
            c = @evalpoly(z, 6//5, -1)
        elseif z < 2/5
            c = @evalpoly(z, 23//15, -2)
        elseif z < 8/15
            c = @evalpoly(z, 121//75, -12//5, 1//2)
        elseif z < 3/5
            c = @evalpoly(z, 427//225, -52//15, 3//2)
        elseif z < 2/3
            c = @evalpoly(z, 4351//2250, -547//150, 9//5, -1//6)
        elseif z < 11/15
            c = @evalpoly(z, 539//250, -647//150, 23//10, -1//6)
        elseif z < 4/5
            c = @evalpoly(z, 7942//3375, -128//25, 17//5, -2//3)
        elseif z < 13/15
            c = @evalpoly(z, 39998//16875, -1952//375, 89//25, -4//5, 1//24)
        elseif z < 14/15
            c = @evalpoly(z, 10109//3750, -1583//250, 243//50, -13//10, 1//24)
        else
            c = @evalpoly(z, 171449//60750, -139199//20250, 2579//450, -173//90, 5//24)
        end

        return c .* u₀
    else
        error("This analytical solution is only valid on (-∞,1]")
    end
end

ff_2delays = DDEFunction(f_2delays,analytic=f_2delays_analytic)

build_prob_dde_2delays(u₀, ::T=u₀) where {T} =
    DDEProblem(ff_2delays, [u₀], (p, t) -> [zero(u₀)], (zero(T), oneunit(T)),
    constant_lags = [T(1//3), T(1//5)])

"""
    prob_dde_2delays

Model problem of finding a solution ``u(t)`` in the time span ``t \\in [0,1]``, where
``t`` is of type `T`, to the delay differential equation

```math
\\frac{du}{dt} = -u(t-1/3)-u(t-1/5)
```

with history function

```math
u(t) = \\begin{cases}
0 & \\text{if } t < 0 \\
u₀ & \\text{if } t = 0
\\end{cases}
```

for ``t \\leq 0`` with ``u₀ = 1``. Note that the problem is discontinuous at ``t = 0`` for
all ``u₀ \\neq 0``.

An analytical solution of this problem is provided for ``t \\in (-\\infty,1]``.
"""
prob_dde_2delays = build_prob_dde_2delays(1.0)

### Not in-place function

function f_2delays_notinplace(u, h, p, t::T) where T
    (- h(p, t - T(1//3)) .- h(p, t - T(1//5))) ./ oneunit(t)
end

ff_2delays_notinplace = DDEFunction(f_2delays_notinplace,analytic=f_2delays_analytic)

#### Vectorized history function

build_prob_dde_2delays_notinplace(u₀, ::T=u₀) where {T} =
    DDEProblem(ff_2delays_notinplace, [u₀], (p, t) -> [zero(u₀)], (zero(T), oneunit(T)),
               constant_lags =  [T(1//3), T(1//5)])

"""
    prob_dde_2delays_notinplace

Same as [`prob_dde_2delays`](@ref), but purposefully implemented with a not in-place
function.
"""
prob_dde_2delays_notinplace = build_prob_dde_2delays_notinplace(1.0)

#### Scalar history function

build_prob_dde_2delays_scalar_notinplace(u₀, ::T=u₀) where {T} =
    DDEProblem(ff_2delays_notinplace, u₀, (p, t) -> zero(u₀), (zero(T), oneunit(T)),
               constant_lags = [T(1//3), T(1//5)])

"""
    prob_dde_2delays_scalar_notinplace

Same as [`prob_dde_2delays_notinplace`](@ref), but purposefully implemented with a
scalar history function.
"""
prob_dde_2delays_scalar_notinplace = build_prob_dde_2delays_scalar_notinplace(1.0)

# DDE examples without analytical solution

## Single constant delay

### In-place function

function f_1delay_long(du, u, h, p, t::T) where T
    du[1] = (- h(p, t - T(1//5))[1] + u[1]) / oneunit(t)
end

build_prob_dde_1delay_long(u₀, ::T=u₀) where {T} =
    DDEProblem(f_1delay_long, [u₀], (p, t) -> [zero(u₀)], (zero(T), T(100)),
    constant_lags = [T(1//5)])

"""
    prob_dde_1delay_long

Model problem of finding a solution ``u(t)`` in the time span ``t \\in [0,100]`` to the
delay differential equation

```math
\\frac{du}{dt} = -u(t-0.2) + u(t)
```

with history function

```math
u(t) = \\begin{cases}
0 & \\text{if } t < 0,\\
u₀ & \\text{if } t = 0,
\\end{cases}
```

for ``t \\leq 0`` with ``u₀ = 1``. Note that the problem is discontinuous at ``t = 0`` for
all ``u₀ \\neq 0``.
"""
prob_dde_1delay_long = build_prob_dde_1delay_long(1.0)

### Not in-place function

function f_1delay_long_notinplace(u, h, p, t::T) where T
    (- h(p, t - T(1//5)) .+ u ) ./ oneunit(t)
end

build_prob_dde_1delay_long_notinplace(u₀, ::T=u₀) where {T} =
    DDEProblem(f_1delay_long_notinplace, [u₀], (p, t) -> [zero(u₀)], (zero(T), T(100)),
               constant_lags = [T(1//5)])

"""
    prob_dde_1delay_long_notinplace

Same as [`prob_dde_1delay_long`](@ref), but purposefully implemented with a not
in-place function.
"""
prob_dde_1delay_long_notinplace = build_prob_dde_1delay_long_notinplace(1.0)

build_prob_dde_1delay_long_scalar_notinplace(u₀, ::T=u₀) where {T} =
    DDEProblem(f_1delay_long_notinplace, u₀, (p, t) -> zero(u₀), (zero(T), T(100)),
               constant_lags = [T(1//5)])

"""
    prob_dde_1delay_long_scalar_notinplace

Same as [`prob_dde_1delay_long_notinplace`](@ref), but purposefully implemented with a
scalar history function.
"""
prob_dde_1delay_long_scalar_notinplace = build_prob_dde_1delay_long_scalar_notinplace(1.0)

## Two constant delays

### In-place function

function f_2delays_long(du, u, h, p, t::T) where T
    du[1] = (- h(p, t - T(1//3))[1] - h(p, t - T(1//5))[1]) / oneunit(t)
end

build_prob_dde_2delays_long(u₀, ::T=u₀) where {T} =
    DDEProblem(f_2delays_long, [u₀], (p, t) -> [zero(u₀)], (zero(T), T(100)),
               constant_lags = [T(1//3), T(1//5)])

"""
    prob_dde_2delays_long

Model problem of finding a solution ``u(t)`` in the time span ``t \\in [0,100]`` to the
delay differential equation

```math
\\frac{du}{dt} = -u(t-1/3) - u(1-1/5)
```

with history function

```math
u(t) = \\begin{cases}
0 & \\text{if } t < 0,\\
u₀ & \\text{if } t = 0,
\\end{cases}
```

for ``t < 0`` with ``u₀ = 1``. Note that the problem is discontinuous at ``t = 0`` for all
``u₀ \\neq 0``.
"""
prob_dde_2delays_long = build_prob_dde_2delays_long(1.0)

### Not in-place function

function f_2delays_long_notinplace(u, h, p, t::T) where T
    (- h(p, t - T(1//3)) .- h(p, t - T(1//5))) ./ oneunit(t)
end

#### Vectorized history function

build_prob_dde_2delays_long_notinplace(u₀, ::T=u₀) where {T} =
     DDEProblem(f_2delays_long_notinplace, [u₀], (p, t) -> [zero(u₀)], (zero(T), T(100)),
                constant_lags = [T(1//3), T(1//5)])

"""
    prob_dde_2delays_long_notinplace

Same as [`prob_dde_2delays_long`](@ref), but purposefully implemented with a not
in-place function.
"""
prob_dde_2delays_long_notinplace = build_prob_dde_2delays_long_notinplace(1.0)

#### Scalar history function

build_prob_dde_2delays_long_scalar_notinplace(u₀, ::T=u₀) where {T} =
     DDEProblem(f_2delays_long_notinplace, u₀, (p, t) -> zero(u₀), (zero(T), T(100)),
                constant_lags = [T(1//3), T(1//5)])

"""
    prob_dde_2delays_long_scalar_notinplace

Same as [`prob_dde_2delays_long_notinplace`](@ref), but purposefully implemented with a
scalar history function.
"""
prob_dde_2delays_long_scalar_notinplace = build_prob_dde_2delays_long_scalar_notinplace(1.0)

#=
The following examples are taken from:
W.H. Enright and H. Hayashi, The evaluation of numerical software for delay
differential equations.
=#

function f_dde_mackey(du, u, h, p, t)
    du[1] = 0.2*h(p, t-14)[1]/(1 + h(p, t-14)[1]^10) - 0.1*u[1]
end

"""
    prob_dde_mackey

Model of blood production with constant delay (M. C. Mackey and L. Glass, Oscillation and
chaos in physiological control systems, 1977).
"""
prob_dde_mackey = DDEProblem(f_dde_mackey, [0.5], (p, t) -> [0.5], (0.0, 500.0),
                             constant_lags = [14])

function f_dde_wheldon(du, u, h, p, t)
    du[1] = 1.1/(1 + sqrt(10)*(h(p, t-20)[1])^(5/4)) - 10*u[1]/(1 + 40*u[2])
    du[2] = 100*u[1]/(1 + 40*u[2]) - 2.43*u[2]
end
u0_wheldon = [1.05767027/3; 1.030713491/3]

"""
    prob_dde_wheldon

Model of chronic granulocytic leukemia with constant delay (T. Wheldon, J. Kirk and
H. Finlay, Cyclical granulopoiesis in chronic granulocytic leukemia: A simulation study, 1974).
"""
prob_dde_wheldon = DDEProblem(f_dde_wheldon, u0_wheldon, (p, t) -> u0_wheldon, (0., 100.),
                              constant_lags = [20])

function f_dde_neves1(du, u, h, p, t)
    du[1] = 1 - h(p, exp(1-1/t))[1]
end
# only valid for specific history function
function f_dde_neves1_analytic(u₀, p, t )
    0 < t ≤ 10 && u₀ == [log(0.1)] && return [log(t)]
    error("This analytical solution is only valid on (0, 10] and for history function ϕ(t) = ln(t) for 0 < t ≤ 0.1")
end

ff_dde_neves1 = DDEFunction(f_dde_neves1,analytic=f_dde_neves1_analytic)

"""
    prob_dde_neves1

DDE with vanishing time dependent delay at ``t = 1`` (K. W. Neves, Automatic integratorion
of functional differential equations: An approach, 1975).
"""
prob_dde_neves_1 = DDEProblem(ff_dde_neves1, [log(0.1)], (p, t) -> [log(t)], (0.1, 10.),
                              dependent_lags = [(u,p,t) -> t - exp(1 - 1/t)])

function f_dde_neves_thompson(du, u, h, p, t)
    if h(p, t/2)[1] < 0
        du[1] = 1 - u[1]
    else
        du[1] = -1 - u[1]
    end
end
# only valid for specific history function
function f_dde_neves_thompson_analytic(u₀, p, t )
    if u₀ == [1]
        if 0 ≤ t ≤ 2*log(2)
            return [2*exp(-t) - 1]
        elseif t ≤ 2*log(6)
            return [1 - 6*exp(-t)]
        elseif t ≤ 2*log(66)
            return [66*exp(-t) - 1]
        end
    end
    error("This analytical solution is only valid on [0, 2ln(66)] and for history function ϕ(t) = 1 for t ≤ 0")
end

ff_dde_neves_thompson = DDEFunction(f_dde_neves_thompson,
                                    analytic=f_dde_neves_thompson_analytic)

"""
    prob_dde_neves_thompson

DDE with vanishing time dependent delay at ``t = 0`` (K. W. Neves and S. Thompson, Solution
of systems of functional differential equations with state dependent delays, 1992), given by

```math
u'(t) = \\begin{cases}
-1 - u(t) & \\text{if } u(t/2) \\geq 0 \\
1 - u(t) & \\text{if } u(t/2) < 0
\\end{case}
```

for ``t \\in [0, 2*\\log 66]`` with history function

```math
u(t) = 1
```

for ``t \\leq 0``.
"""
prob_dde_neves_thompson = DDEProblem(ff_dde_neves_thompson, [1.],
                                     (p, t) -> [1.], (0., 2*log(66)),
                                     constant_lags = [(u,p,t) -> t/2])

function f_dde_paul1(du, u, h, p, t)
    du[1] = - 2*h(p, t - 1 - abs(u[1]))[1]*(1 - u[1]^2)
end

"""
    prob_dde_paul1

DDE with state dependent delay (C. A. H. Paul, A test set of functional differential
equations, 1994), given by

```math
u'(t) = - 2u(t - 1 - |u(t)|)(1 - u(t)^2)
```

for ``t \\in [0, 30]`` with history function

```math
u(t) = 1/2
```

for ``t \\leq 0``.
"""
prob_dde_paul1 = DDEProblem(f_dde_paul1, [0.5], (p, t) -> [0.5], (0., 30.),
                            dependent_lags = [(u,p,t) -> 1 + abs(u[1])])

function f_dde_paul2(du, u, h, p, t)
    h1 = h(p, t - u[2])[1]
    du[1] = -2*h1
    du[2] = (abs(h1) - abs(u[1]))/(1 + abs(h1))
end

"""
    prob_dde_paul2

DDE with state dependent delay (C. A. H. Paul, A test set of functional differential equations, 1994).
"""
prob_dde_paul2 = DDEProblem(f_dde_paul2, [1; 0.5], (p, t) -> [1; 0.5], (0., 30.),
                            dependent_lags = [(u,p,t) -> u[2]])

function build_prob_dde_mahaffy(tspan, h, σ₀, T₁, γ, Q, k, a, K, r)
    function f(du, u, h, p, t)
        du[1] = σ₀*h(p, t-T₁)[2] - γ*u[1] - Q
        du[2] = a/(1 + K*u[1]^r) - k*u[2]
        du[3] = 1 - Q*exp(γ*u[3])/(σ₀*h(p, t-T₁-u[3])[2])
    end

    DDEProblem(f, h(nothing, 0), h, tspan,
                constant_lags = [T₁],
                dependent_lags = [(u,p,t) -> T₁ + u[3]])
end

function h_mahaffy1(p, t)
    if t < -6
        return [3.325; 9.5; 120]
    else
        return [3.325; 10; 120]
    end
end
prob_dde_mahaffy1 = build_prob_dde_mahaffy((0., 300.), h_mahaffy1, 0.0031, 6, 0.001, 0.0275,
                                           2.8, 6570, 0.0382, 6.96)

prob_dde_mahaffy2 = build_prob_dde_mahaffy((0., 100.), (p, t) -> [3.5; 10; 50], 0.00372, 3, 0.1,
                                           0.00178, 6.65, 15600, 0.0382, 6.96)

"""
Two models of hematopoiesis with state dependent delay (J. M. Mahaffy, J. Bélair and
M. C. Mackey, Hematopoietic model with moving boundary condition and state dependent delay, 1996).
"""
prob_dde_mahaffy1, prob_dde_mahaffy2

function f_dde_neves2(du, u, h, p, t)
    t2 = exp(1 - u[2])
    du[1] = u[2]
    du[2] = -h(p, t2)[2]*u[2]^2*t2
end
# only valid for specific history function
function f_dde_neves2(::Type{Val{:analytic}}, u₀, p, t )
    0 < t ≤ 5 && u₀ == [log(0.1); 10] && return [log(t); 1/t]
    error("This analytical solution is only valid on (0, 5) and for history function ϕ(t) = [ln(t); 1/t] for 0 < t ≤ 0.1")
end

"""
    prob_dde_neves2

DDE with vanishing state dependent delay at ``t = 1`` (K. W. Neves, Automatic integration
of functional differential equations: An approach, 1975).
"""
prob_dde_neves2 = DDEProblem(f_dde_neves2, [log(0.1); 10], (p, t) -> [log(t); 1/t], (0.1, 5.),
                             dependent_lags = [(u,p,t) -> t - exp(1 - u[2])])

function build_f_dde_gatica(r₁, r₂, α, δ)
    function f_dde_gatica(du, u, h, p, t)
        u₁u₂ = u[1]*u[2]
        r₁u₁u₂ = r₁*u₁u₂
        r₂u₃ = r₂*u[3]
        uₕ = h(p, t - u[4])
        uₕ₁₂ = uₕ[1]*uₕ[2]

        du[1] = -r₁u₁u₂ + r₂u₃
        du[2] = -r₁u₁u₂ + α*r₁*uₕ₁₂
        du[3] = r₁u₁u₂ - r₂u₃
        du[4] = 1 + (3*δ - u₁u₂ - u[3])/(uₕ₁₂ + uₕ[3])*exp(δ*u[4])
    end
end

"""
    prob_dde_gatica

Model of antigen antibody dynamics with fading memory, with vanishing state dependent delay
at ``t = 0`` (J. Gatica and P. Waltman, A threshold model of antigen antibody dynamics with fading memory, 1982).
"""
prob_dde_gatica = DDEProblem(build_f_dde_gatica(0.02, 0.005, 3, 0.01),
                             [5; 0.1; 0; 0], (p, t) -> [5; 0.1; 0; 0],
                             (0., 40.), dependent_lags = [(u,p,t) -> u[4]])

#=
Quorum Sensing model
=#
function build_prob_dde_qs(u₀, tspan, τ, D, γₛ, Kₘ, nₛ, a, αₐ, βₐ, C₁, n₁, γₐ, αC, R, γ₃, Kₑ, αₗ, C₂, n₂, γₗ)
    function hill_growth(x, K, n)
        if x > K
            return 1/(1 + (K / x)^n)
        end

        if K == 0
            return one(x)
        end

        y = (x / K)^n
        return y / (1 + y)
    end

    S₀ = u₀[1]

    function f_dde_qs(du, u, h, p, t)
        if u[1] < 0
            # eq. 1 and 2 not defined for u[1] < 0
            du[1] = 0
            du[2] = 0
        else
            tmp = u[2] * hill_growth(u[1], Kₘ, nₛ)
            du[1] = D * (S₀ - u[1]) - γₛ * tmp
            du[2] = a * tmp - D * u[2]
        end

        if u[4] < 0
            # eq. 3 not defined for u[4] < 0
            du[3] = 0
            du[4] = max(0, αC * (R - u[4]) * u[3] - γ₃ * u[4])
        else
            du[3] = (αₐ + βₐ * hill_growth(u[4], C₁, n₁)) * u[2] -
                (γₐ + D) * u[3] - αC * (R - u[4]) * u[3] + γ₃ * u[4] -
                Kₑ * u[3] * u[5]
            du[4] = αC * (R - u[4]) * u[3] - γ₃ * u[4]
        end

        Cₜ = h(p, t - τ)[4]
        if Cₜ < 0
            # eq. 5 not defined for Cₜ < 0
            du[5] = 0
        else
            du[5] = αₗ * hill_growth(Cₜ, C₂, n₂) * u[2] - (γₗ + D) * u[5]
        end
    end

    DDEProblem(f_dde_qs, u₀, (p, t) -> u₀, tspan, constant_lags = [τ])
end

"""
    prob_dde_qs

Model of Quorum Sensing (QS) of Pseudomonas putida IsoF in continuous cultures with constant
delay (Buddrus-Schiemann et al., Analysis of N-Acylhomoserine Lactone Dynamics in Continuous
Cultures of Pseudomonas Putida IsoF By Use of ELISA and UHPLC/qTOF-MS-derived Measurements
and Mathematical Models, Analytical and Bioanalytical Chemistry, 2014).
"""
prob_dde_qs = build_prob_dde_qs([1; 8.4e8; 2.5e-9; 7.6e-8; 5e-15], # initial values
                                (0., 45.), # time span
                                2, # delay
                                0.1, 1.3e-12, 0.38, 1.3, 0.66, 2.3e-19, 2.3e-18, 70e-9,
                                2.3, 0.05, 4e4, 5e-7, 0.080, 1.5e-4, 1.1e-8, 70e-9, 2.5,
                                0.005)

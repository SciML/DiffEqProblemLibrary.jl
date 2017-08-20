# DDE examples with analytical solution

## Single constant delay

### In-place function

f_1delay = function (t,u,h,du)
    du[1] = - h(t-1)[1]
end

function (f::typeof(f_1delay))(::Type{Val{:analytic}}, t, u₀)
    if t < 0
        return 0
    elseif t < 1
        return u₀
    elseif t < 2
        return u₀ * (2 - t)
    elseif t < 3
        return u₀ * (8 - 6t + t^2) / 2
    elseif t < 4
        return u₀ * (51 - 45t + 12t^2 - t^3) / 6
    elseif t < 5
        return u₀ * (460 - 436t + 144t^2 - 20t^3 + t^4) / 24
    elseif t < 6
        return u₀ * (5_425 - 5_305t + 1_970t^2 - 350t^3 + 30t^4 - t^5) / 120
    elseif t < 7
        return u₀ * (79_206 - 78_486t + 31_260t^2 - 6_420t^3 + 720t^4 - 42t^5 + t^6) / 720
    elseif t < 8
        return u₀ * (1_377_985 - 1_372_945t + 571_767t^2 - 128_975t^3 + 17_045t^4 -
                     1_323t^5 + 56t^6 - t^7) / 5040
    elseif t < 9
        return u₀ * (27_801_096 - 27_760_776t + 11_914_168t^2 - 2_866_808t^3 + 423_080t^4 -
                     39_256t^5 + 2_240t^6 - 72t^7 + t^8)  / 40_320
    elseif t ≤ 10
        return u₀ * (637_630_353 - 637_267_473t + 279_414_396t^2 - 70_442_316t^3 +
                     11_247_894t^4 - 1_179_990t^5 + 81_396t^6 - 3_564t^7 + 90t^8 - t^9) /
                     362_880
    else
        error("This analytical solution is only valid on (-∞,10]")
    end
end

"""
    prob_dde_1delay(u₀)

Model problem of finding a solution ``u(t)`` in the time span ``t \\in [0,10]`` to the
delay differential equation

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

for ``t \\leq 0``. Hence the problem is discontinuous at ``t = 0`` for all ``u₀ \\neq 0``.

An analytical solution of this problem is provided for ``t \\in (-\\infty,10]``.
"""
prob_dde_1delay(u₀) = DDEProblem(f_1delay, t->[0.0], [u₀], (0.0, 10.0), [1])

### Not in-place function

f_1delay_notinplace = function (t,u,h)
    - h(t-1)
end

(f::typeof(f_1delay_notinplace))(::Type{Val{:analytic}}, t, u0) = f_1delay(Val{:analytic},
                                                                           t, u0)

#### Vectorized history function

"""
    prob_dde_1delay_notinplace(u₀)

Same as [`prob_dde_1delay`](@ref), but purposefully implemented with a not in-place
function.
"""
prob_dde_1delay_notinplace(u₀) =
    DDEProblem(f_1delay_notinplace, t->[0.0], [u₀], (0.0, 10.0), [1])

#### Scalar history function

"""
    prob_dde_1delay_scalar_notinplace(u₀)

Same as [`prob_dde_1delay_notinplace`](@ref), but purposefully implemented with a scalar
history function.
"""
prob_dde_1delay_scalar_notinplace(u₀) =
    DDEProblem(f_1delay_notinplace, t->0.0, u₀, (0.0, 10.0), [1])

## Two constant delays

### In-place function

f_2delays = function (t,u,h,du)
    du[1] = - h(t-1/3)[1] - h(t-1/5)[1]
end

function (f::typeof(f_2delays))(::Type{Val{:analytic}}, t, u₀)
    if t < 0
        return 0
    elseif t < 1/5
        return u₀
    elseif t < 1/3
        return u₀ * (6 - 5t) / 5
    elseif t < 2/5
        return u₀ * (23 - 30t) / 15
    elseif t < 8/15
        return u₀ * (242 - 360t + 75t^2) / 150
    elseif t < 3/5
        return u₀ * (854 - 1_560t + 675t^2) / 450
    elseif t < 2/3
        return u₀ * (4_351 - 8_205t + 4_050t^2 - 375t^3) / 2_250
    elseif t < 11/15
        return u₀ * (1_617 - 3_235t + 1_725t^2 - 125t^3) / 750
    elseif t < 4/5
        return u₀ * (7_942 - 17_280t + 11_475t^2 - 2_250t^3) / 3_375
    elseif t < 13/15
        return u₀ * (319_984 - 702_720t + 480_600t^2 - 108_000t^3 + 5_625t^4) / 135_000
    elseif t < 14/15
        return u₀ * (40_436 - 94_980t + 72_900t^2 - 19_500t^3 + 625t^4) / 15_000
    elseif t ≤ 1
        return u₀ * (685_796 - 1_670_388t + 1_392_660t^2 - 467_100t^3 + 50_625t^4) / 243_000
    else
        error("This analytical solution is only valid on (-∞,1]")
    end
end

"""
    prob_dde_2delays(u₀)

Model problem of finding a solution ``u(t)`` in the time span ``t \\in [0,1]`` to the delay
differential equation

```math
\\frac{du}{dt} = -u(t-1/3)-u(t-1/5)
```

with history function

```math
u(t) = \\begin{cases}
0 & \text{if } t < 0 \\
u₀ & \text{if } t = 0
\end{cases}
```

for ``t \\leq 0``. Hence the problem is discontinuous at ``t = 0`` for all ``u₀ \\neq 0``.

An analytical solution of this problem is provided for ``t \\in (-\\infty,1]``.
"""
prob_dde_2delays(u₀) = DDEProblem(f_2delays, t->[0.0], [u₀], (0.0, 1.0), [1//3, 1//5])

### Not in-place function

f_2delays_notinplace = function (t,u,h)
    - h(t-1/3) - h(t-1/5)
end

(f::typeof(f_2delays_notinplace))(::Type{Val{:analytic}}, t, u0) =
    f_2delays(Val{:analytic}, t, u0)

#### Vectorized history function

"""
    prob_dde_2delays_notinplace(u₀)

Same as [`prob_dde_2delays`](@ref), but purposefully implemented with a not in-place
function.
"""
prob_dde_2delays_notinplace(u₀) =
    DDEProblem(f_2delays_notinplace, t->[0.0], [u₀], (0.0, 1.0), [1//3, 1//5])

#### Scalar history function

"""
    prob_dde_2delays_scalar_notinplace(u₀)

Same as [`prob_dde_2delays_notinplace`](@ref), but purposefully implemented with a scalar
history function.
"""
prob_dde_2delays_scalar_notinplace(u₀) =
    DDEProblem(f_2delays_notinplace, t->0.0, u₀, (0.0, 1.0), [1//3, 1//5])

# DDE examples without analytical solution

## Single constant delay

### In-place function

f_1delay_long = function (t,u,h,du)
    du[1] = - h(t-0.2)[1] + u[1]
end

"""
Model problem of finding a solution ``u(t)`` in the time span ``t \\in [0,100]`` to the
delay differential equation

```math
\\frac{du}{dt} = -u(t-0.2) + u(t)
```

with history function

```math
u(t) = \\begin{cases}
0 & \\text{if } t < 0,\\
1 & \\text{if } t = 0,
\\end{cases}
```

for ``t \\leq 0``. Hence the problem is discontinuous at ``t = 0``.
"""
prob_dde_1delay_long = DDEProblem(f_1delay_long, t->[0.0], [1.0],
                                             (0.0, 100.0),[0.2])

### Not in-place function

f_1delay_long_notinplace = function (t,u,h)
    - h(t-0.2) + u
end

"""
Same as [`prob_dde_1delay_long`](@ref), but purposefully implemented with a not in-place
function.
"""
prob_dde_1delay_long_notinplace =
    DDEProblem(f_1delay_long_notinplace, t->[0.0], [1.0], (0.0, 100.0), [0.2])

"""
Same as [`prob_dde_1delay_long_notinplace`](@ref), but purposefully implemented with a
scalar history function.
"""
prob_dde_1delay_long_scalar_notinplace =
    DDEProblem(f_1delay_long_notinplace, t->0.0, 1.0, (0.0, 100.0), [0.2])

## Two constant delays

### In-place function

f_2delays_long = function (t,u,h,du)
    du[1] = - h(t-1/3)[1] - h(t-1/5)[1]
end

"""
Model problem of finding a solution ``u(t)`` in the time span ``t \\in [0,100]`` to the
delay differential equation

```math
\\frac{du}{dt} = -u(t-1/3) - u(1-1/5)
```

with history function

```math
u(t) = \\begin{cases}
0 & \\text{if } t < 0,\\
1 & \\text{if } t = 0,
\\end{cases}
```

for ``t < 0``. Hence the problem is discontinuous at ``t = 0``.
"""
prob_dde_2delays_long = DDEProblem(f_2delays_long, t->[0.0], [1.0],
                                              (0.0, 100.0), [1//3, 1//5])

### Not in-place function

f_2delays_long_notinplace = function (t,u,h)
    - h(t-1/3) - h(t-1/5)
end

#### Vectorized history function

"""
Same as [`prob_dde_2delays_long`](@ref), but purposefully implemented with a not in-place
function.
"""
prob_dde_2delays_long_notinplace =
    DDEProblem(f_2delays_long_notinplace, t->[0.0], [1.0],
                          (0.0, 100.0), [1//3, 1//5])

#### Scalar history function

"""
Same as [`prob_dde_2delays_long_notinplace`](@ref), but purposefully implemented with a scalar
history function.
"""
prob_dde_2delays_long_scalar_notinplace =
    DDEProblem(f_2delays_long_notinplace, t->0.0, 1.0,
                          (0.0, 100.0), [1//3, 1//5])

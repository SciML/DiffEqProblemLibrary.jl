# Helper function for changing initial value and time type of the following DDEs

"""
    remake_dde_constant_u0_tType(prob::DDEProblem, u₀, tType)

Create a new delay differential problem by replacing the initial state of problem `prob`
with `u0` and setting the type of time `t` to `tType`.

This function makes special assumptions about the structure of `prob` and is intended to be
used for all problems of name `prob_dde_constant_*` in `DDEProblemLibrary`. The functions of
these delay differential equation problems with constant delays are purposefully
implemented such that they work for arbitrary types of state `u` and time `t`, and hence in
particular for number types with units. The element type of `u` is saved as parameter `p` to
ensure that the return type of the history functions is correct (which are functions of `p`
but not of `u`).
"""
function remake_dde_constant_u0_tType(prob::DDEProblem, u₀, tType)
    remake(prob; u0 = u₀, tspan = tType.(prob.tspan), p = eltype(u₀),
        constant_lags = tType.(prob.constant_lags))
end

# History functions

function h_dde_constant_ip(p, t; idxs = 1)
    idxs == 1 || error("history function is only implemented for the first component")
    t ≤ zero(t) || error("history function is only implemented for t ≤ 0")

    p === nothing ? zero(t) : zero(p)
end

function h_dde_constant_oop(p, t)
    t ≤ zero(t) || error("history function is only implemented for t ≤ 0")

    p === nothing ? [zero(t)] : [zero(p)]
end

function h_dde_constant_scalar(p, t)
    t ≤ zero(t) || error("history function is only implemented for t ≤ 0")

    p === nothing ? zero(t) : zero(p)
end

# Single constant delay

## Short time span

### In-place function

@doc raw"""
    prob_dde_constant_1delay_ip

Delay differential equation

```math
u'(t) = -u(t - 1)
```

for ``t \in [0, 1]`` with history function ``\phi(t) = 0`` if ``t < 0`` and ``\phi(0) = 1``.

# Solution

The analytical solution for ``t \in [0, 10]`` can be obtained by the method of steps and
is provided in this implementation.
"""
prob_dde_constant_1delay_ip

function f_dde_constant_1delay_ip!(du, u, h, p, t)
    e = oneunit(t)
    du[1] = -h(p, t - e; idxs = 1) / e

    nothing
end

function fanalytic_dde_constant_1delay(u₀,
        ::Union{typeof(h_dde_constant_ip),
            typeof(h_dde_constant_oop),
            typeof(h_dde_constant_scalar)},
        p, t)
    z = t * inv(oneunit(t))
    0 ≤ z ≤ 10 || error("analytical solution is only implemented for t ∈ [0, 10]")

    if z < 1
        copy(u₀)
    else
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

        c .* u₀
    end
end

const prob_dde_constant_1delay_ip = DDEProblem(
    DDEFunction(f_dde_constant_1delay_ip!;
        analytic = fanalytic_dde_constant_1delay),
    [1.0], h_dde_constant_ip, (0.0, 10.0),
    typeof(1.0);
    constant_lags = [1])

### Out-of-place function

"""
    prob_dde_constant_1delay_oop

Same delay differential equation as [`prob_dde_constant_1delay_ip`](@ref), but purposefully
implemented with an out-of-place function.
"""
prob_dde_constant_1delay_oop

function f_dde_constant_1delay_oop(u, h, p, t)
    e = oneunit(t)
    h(p, t - e) ./ (-e)
end

const prob_dde_constant_1delay_oop = DDEProblem(
    DDEFunction(f_dde_constant_1delay_oop;
        analytic = fanalytic_dde_constant_1delay),
    [1.0], h_dde_constant_oop, (0.0, 10.0),
    typeof(1.0);
    constant_lags = [1])

### Scalar function

"""
    prob_dde_constant_1delay_scalar

Same delay differential equation as [`prob_dde_constant_1delay_ip`](@ref), but purposefully
implemented with a scalar function.
"""
prob_dde_constant_1delay_scalar

function f_dde_constant_1delay_scalar(u, h, p, t)
    e = oneunit(t)
    -h(p, t - e) / e
end

const prob_dde_constant_1delay_scalar = DDEProblem(
    DDEFunction(f_dde_constant_1delay_scalar;
        analytic = fanalytic_dde_constant_1delay),
    1.0, h_dde_constant_scalar, (0.0, 10.0),
    typeof(1.0);
    constant_lags = [1])

## Long time span

### In-place function

@doc raw"""
    prob_dde_constant_1delay_long_ip

Delay differential equation

```math
u'(t) = u(t) - u(t - 1/5)
```

for ``t \in [0, 100]`` with history function ``\phi(t) = 0`` if ``t < 0`` and
``\phi(0) = 1``.
"""
prob_dde_constant_1delay_long_ip

function f_dde_constant_1delay_long_ip!(du, u, h, p, t)
    T = typeof(t)
    du[1] = (u[1] - h(p, t - T(1 / 5); idxs = 1)) / oneunit(t)

    nothing
end

const prob_dde_constant_1delay_long_ip = DDEProblem(f_dde_constant_1delay_long_ip!, [1.0],
    h_dde_constant_ip, (0.0, 100.0),
    typeof(1.0);
    constant_lags = [1 / 5])

### Out-of-place function

"""
    prob_dde_constant_1delay_long_oop

Same delay differential equation as [`prob_dde_constant_1delay_long_ip`](@ref), but
purposefully implemented with an out-of-place function.
"""
prob_dde_constant_1delay_long_oop

function f_dde_constant_1delay_long_oop(u, h, p, t)
    T = typeof(t)
    (u .- h(p, t - T(1 / 5))) ./ oneunit(t)
end

const prob_dde_constant_1delay_long_oop = DDEProblem(f_dde_constant_1delay_long_oop, [1.0],
    h_dde_constant_oop, (0.0, 100.0),
    typeof(1.0);
    constant_lags = [1 / 5])

### Scalar function

"""
    prob_dde_constant_1delay_long_scalar

Same delay differential equation as [`prob_dde_constant_1delay_long_ip`](@ref), but
purposefully implemented with a scalar function.
"""
prob_dde_constant_1delay_long_scalar

function f_dde_constant_1delay_long_scalar(u, h, p, t)
    T = typeof(t)
    (u - h(p, t - T(1 / 5))) / oneunit(t)
end

const prob_dde_constant_1delay_long_scalar = DDEProblem(f_dde_constant_1delay_long_scalar,
    1.0, h_dde_constant_scalar,
    (0.0, 100.0),
    typeof(1.0);
    constant_lags = [1 / 5])

# Two constant delays

## Short time span

### In-place function

@doc raw"""
    prob_dde_constant_2delays_ip

Delay differential equation

```math
u'(t) = -u(t - 1/3) - u(t - 1/5)
```

for ``t \in [0, 1]`` with history function ``\phi(t) = 0`` if ``t < 0`` and ``\phi(0) = 1``.

# Solution

The analytical solution for ``t \in [0, 10]`` can be obtained by the method of steps and
is provided in this implementation.
"""
prob_dde_constant_2delays_ip

function f_dde_constant_2delays_ip!(du, u, h, p, t)
    T = typeof(t)
    du[1] = -(h(p, t - T(1 / 3); idxs = 1) + h(p, t - T(1 / 5); idxs = 1)) / oneunit(t)

    nothing
end

function fanalytic_dde_constant_2delays(u₀,
        ::Union{typeof(h_dde_constant_ip),
            typeof(h_dde_constant_oop),
            typeof(h_dde_constant_scalar)},
        p, t)
    z = t * inv(oneunit(t))
    0 ≤ z ≤ 1 || error("analytical solution is only implemented for t ∈ [0, 1]")

    if z < 1 / 5
        copy(u₀)
    else
        if z < 1 / 3
            c = @evalpoly(z, 6//5, -1)
        elseif z < 2 / 5
            c = @evalpoly(z, 23//15, -2)
        elseif z < 8 / 15
            c = @evalpoly(z, 121//75, -12//5, 1//2)
        elseif z < 3 / 5
            c = @evalpoly(z, 427//225, -52//15, 3//2)
        elseif z < 2 / 3
            c = @evalpoly(z, 4351//2250, -547//150, 9//5, -1//6)
        elseif z < 11 / 15
            c = @evalpoly(z, 539//250, -647//150, 23//10, -1//6)
        elseif z < 4 / 5
            c = @evalpoly(z, 7942//3375, -128//25, 17//5, -2//3)
        elseif z < 13 / 15
            c = @evalpoly(z, 39998//16875, -1952//375, 89//25, -4//5, 1//24)
        elseif z < 14 / 15
            c = @evalpoly(z, 10109//3750, -1583//250, 243//50, -13//10, 1//24)
        else
            c = @evalpoly(z, 171449//60750, -139199//20250, 2579//450, -173//90, 5//24)
        end

        c .* u₀
    end
end

const prob_dde_constant_2delays_ip = DDEProblem(
    DDEFunction(f_dde_constant_2delays_ip!;
        analytic = fanalytic_dde_constant_2delays),
    [1.0], h_dde_constant_ip, (0.0, 1.0),
    typeof(1.0);
    constant_lags = [1 / 3, 1 / 5])

### Out-of-place function

"""
    prob_dde_constant_2delays_oop

Same delay differential equation as [`prob_dde_constant_2delays_ip`](@ref), but purposefully
implemented with an out-of-place function.
"""
prob_dde_constant_2delays_oop

function f_dde_constant_2delays_oop(u, h, p, t)
    T = typeof(t)
    (h(p, t - T(1 / 3)) .+ h(p, t - T(1 / 5))) ./ (-oneunit(t))
end

const prob_dde_constant_2delays_oop = DDEProblem(
    DDEFunction(f_dde_constant_2delays_oop;
        analytic = fanalytic_dde_constant_2delays),
    [1.0], h_dde_constant_oop, (0.0, 1.0),
    typeof(1.0);
    constant_lags = [1 / 3, 1 / 5])

### Scalar function

"""
    prob_dde_constant_2delays_scalar

Same delay differential equation as [`prob_dde_constant_2delays_ip`](@ref), but purposefully
implemented with a scalar function.
"""
prob_dde_constant_2delays_scalar

function f_dde_constant_2delays_scalar(u, h, p, t)
    T = typeof(t)
    -(h(p, t - T(1 / 3)) + h(p, t - T(1 / 5))) / oneunit(t)
end

const prob_dde_constant_2delays_scalar = DDEProblem(
    DDEFunction(f_dde_constant_2delays_scalar;
        analytic = fanalytic_dde_constant_2delays),
    1.0, h_dde_constant_scalar, (0.0, 1.0),
    typeof(1.0);
    constant_lags = [1 / 3, 1 / 5])

## Long time span

### In-place function

@doc raw"""
    prob_dde_constant_2delays_long_ip

Delay differential equation

```math
u'(t) = - u(t - 1/3) - u(t - 1/5)
```

for ``t \in [0, 100]`` with history function ``\phi(t) = 0`` if ``t < 0`` and
``\phi(0) = 1``.
"""
prob_dde_constant_2delays_long_ip

function f_dde_constant_2delays_long_ip!(du, u, h, p, t)
    T = typeof(t)
    du[1] = -(h(p, t - T(1 / 3); idxs = 1) + h(p, t - T(1 / 5); idxs = 1)) / oneunit(t)
    nothing
end

const prob_dde_constant_2delays_long_ip = DDEProblem(
    f_dde_constant_2delays_long_ip!, [1.0],
    h_dde_constant_ip, (0.0, 100.0),
    typeof(1.0);
    constant_lags = [1 / 3, 1 / 5])

### Out-of-place function

"""
    prob_dde_constant_2delays_long_oop

Same delay differential equation as [`prob_dde_constant_2delays_long_ip`](@ref), but
purposefully implemented with an out-of-place function.
"""
prob_dde_constant_2delays_long_oop

function f_dde_constant_2delays_long_oop(u, h, p, t)
    T = typeof(t)
    (h(p, t - T(1 / 3)) .+ h(p, t - T(1 / 5))) ./ (-oneunit(t))
end

const prob_dde_constant_2delays_long_oop = DDEProblem(f_dde_constant_2delays_long_oop,
    [1.0], h_dde_constant_oop,
    (0.0, 100.0),
    typeof(1.0);
    constant_lags = [1 / 3, 1 / 5])

#### Scalar function

"""
    prob_dde_constant_2delays_long_scalar

Same delay differential equation as [`prob_dde_constant_2delays_long_ip`](@ref), but
purposefully implemented with a scalar function.
"""
prob_dde_constant_2delays_long_scalar

function f_dde_constant_2delays_long_scalar(u, h, p, t)
    T = typeof(t)
    -(h(p, t - T(1 / 3)) + h(p, t - T(1 / 5))) / oneunit(t)
end

const prob_dde_constant_2delays_long_scalar = DDEProblem(
    f_dde_constant_2delays_long_scalar,
    1.0, h_dde_constant_scalar,
    (0.0, 100.0),
    typeof(1.0);
    constant_lags = [1 / 3, 1 / 5])

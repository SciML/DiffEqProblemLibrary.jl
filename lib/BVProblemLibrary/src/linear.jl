λ = 1

################### linear_bvp1 ############################
function prob_bvp_linear_1_analytic(u, λ, t)
    a = 1 / sqrt(λ)
    [(exp(-a * t) - exp((t - 2) * a)) / (1 - exp(-2 * a)),
        (-a * exp(-t * a) - a * exp((t - 2) * a)) / (1 - exp(-2 * a))]
end
function prob_bvp_linear_1_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * u[1]
end
function prob_bvp_linear_1_bca!(resid_a, u_a, p)
    resid_a[1] = u_a[1] - 1
end
function prob_bvp_linear_1_bcb!(resid_b, u_b, p)
    resid_b[1] = u_b[1]
end

prob_bvp_linear_1_function = BVPFunction(
    prob_bvp_linear_1_f!, (prob_bvp_linear_1_bca!, prob_bvp_linear_1_bcb!),
    analytic = prob_bvp_linear_1_analytic, bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_linear_1_tspan = (0.0, 1.0)
@doc raw"""
    prob_bvp_linear_1

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}y_1
```

with boundary condition

```math
y_1(0)=1, y_1(1)=0
```

# Solution

The analytical solution for ``t \in [0, 1]`` is

```math
y_1(t) = \frac{\exp(-t/\sqrt{\lambda}) - \exp((t-2)/\sqrt{\lambda})}{1-\exp(-2/\sqrt{\lambda})}
```
```math
y_2(t)=y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=75)
"""
prob_bvp_linear_1 = BVProblem(prob_bvp_linear_1_function,
    [1.0, 0.0],
    prob_bvp_linear_1_tspan,
    λ)

################### linear_bvp2 ############################
function prob_bvp_linear_2_analytic(u, λ, t)
    a = 1 / sqrt(λ)
    [(1 - exp((t - 1) * λ)) / (1 - exp(-a)),
        -a * exp((t - 1) * a) / (1 - exp(-a))]
end
function prob_bvp_linear_2_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * u[2]
end
function prob_bvp_linear_2_bca!(resid_a, u_a, p)
    resid_a[1] = u_a[1] - 1
end
function prob_bvp_linear_2_bcb!(resid_b, u_b, p)
    resid_b[1] = u_b[1]
end
prob_bvp_linear_2_function = BVPFunction(
    prob_bvp_linear_2_f!, (prob_bvp_linear_2_bca!, prob_bvp_linear_2_bcb!),
    analytic = prob_bvp_linear_2_analytic, bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_linear_2_tspan = (0.0, 1.0)
@doc raw"""
    prob_bvp_linear_2

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}y_2
```

with boundary condition

```math
y_1(0)=1, y_1(1)=0
```

# Solution

The analytical solution for ``t \in [0, 1]`` is

```math
y_1(t) = \frac{1-\exp((t-1)/\lambda)}{1-\exp(-1/\lambda)}
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=103)
"""
prob_bvp_linear_2 = BVProblem(prob_bvp_linear_2_function,
    [1.0, 0.0],
    prob_bvp_linear_2_tspan,
    λ)

################### linear_bvp3 ############################
function prob_bvp_linear_3_analytic(u, λ, t)
    [cos(π * t),
        -π * sin(π * t)]
end
function prob_bvp_linear_3_f(t, u1, u2)
    -(2 + cos(π * t)) * u2 + u1 - (1 + λ * π^2) * cos(π * t) -
    (2 + cos(π * t)) * π * sin(π * t)
end
function prob_bvp_linear_3_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_3_f(t, u[1], u[2])
end
function prob_bvp_linear_3_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] + 1
end
function prob_bvp_linear_3_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] + 1
end
prob_bvp_linear_3_function = BVPFunction(
    prob_bvp_linear_3_f!, (prob_bvp_linear_3_bca!, prob_bvp_linear_3_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_3_analytic, twopoint = Val(true))
prob_bvp_linear_3_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_3

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(t, y_1, y_2)
```

where

```math
f(t, y_1, y_2) = -(2+\cos(\pi t))y_2 + y_1 -(1+\lambda \pi^2)\cos(\pi t) - (2+\cos(\pi t))\pi\sin(\pi t)
```

with boundary condition

```math
y_1(-1)=-1, y_1(1)=-1
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = \cos(\pi t)
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=111)
"""
prob_bvp_linear_3 = BVProblem(prob_bvp_linear_3_function,
    [1.0, 0.0],
    prob_bvp_linear_3_tspan,
    λ)

################### linear_bvp4 ############################
function prob_bvp_linear_4_analytic(u, λ, t)
    [exp(t - 1) + exp(-(1 + λ) * (1 + t) / λ),
        exp(t - 1) - (1 + λ) / λ * exp(-(1 + λ) * (1 + t) / λ)]
end
prob_bvp_linear_4_f(u1, u2, λ) = -u2 + (1 + λ) * u1
function prob_bvp_linear_4_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_4_f(u[1], u[2], p)
end
function prob_bvp_linear_4_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1 - exp(-2)
end
function prob_bvp_linear_4_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1 - exp(-2 * (1 + p) / p)
end
prob_bvp_linear_4_function = BVPFunction(
    prob_bvp_linear_4_f!, (prob_bvp_linear_4_bca!, prob_bvp_linear_4_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_4_analytic, twopoint = Val(true))
prob_bvp_linear_4_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_4

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(y_1, y_2)
```

where

```math
f(y_1, y_2)=-y2+(1+\lambda)y1
```

with boundary condition

```math
y_1(-1)=1+\exp(-2), y_1(1)=1+\exp(-2(1+\lambda))
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = \exp(t-1)+\exp(-(1+\lambda)(1+t)/\lambda)
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=171)
"""
prob_bvp_linear_4 = BVProblem(prob_bvp_linear_4_function,
    [1.0, 0.0],
    prob_bvp_linear_4_tspan,
    λ)

################### linear_bvp5 ############################
function prob_bvp_linear_5_analytic(u, λ, t)
    [cos(π * t),
        -π * sin(π * t)]
end
function prob_bvp_linear_5_f(t, u1, u2, λ)
    t * u2 + u1 - (1 + λ * π^2) * cos(π * t) + π * t * sin(π * t)
end
function prob_bvp_linear_5_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_5_f(t, u[1], u[2], p)
end
function prob_bvp_linear_5_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] + 1
end
function prob_bvp_linear_5_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] + 1
end
prob_bvp_linear_5_function = BVPFunction(
    prob_bvp_linear_5_f!, (prob_bvp_linear_5_bca!, prob_bvp_linear_5_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_5_analytic, twopoint = Val(true))
prob_bvp_linear_5_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_5

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(t, y_1, y_2)
```

where

```math
f(t, y_1, y_2)=ty_2+y_1-(1+\lambda\pi^2)\cos(\pi t)+\pi t\sin(\pi t)
```

with boundary condition

```math
y_1(-1)=-1, y_1(1)=-1
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = \cos(\pi t)
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=174)
"""
prob_bvp_linear_5 = BVProblem(prob_bvp_linear_5_function,
    [1.0, 0.0],
    prob_bvp_linear_5_tspan,
    λ)

################### linear_bvp6 ############################
function prob_bvp_linear_6_analytic(u, λ, t)
    [cos(π * t) + erf(t / sqrt(2 * λ)) / erf(1 / sqrt(2 * λ)),
        -π * sin(π * t) +
        1 / sqrt(2 * λ) * 2 / sqrt(π) * exp(-t^2 / (2 * λ)) / erf(1 / sqrt(2 * λ))]
end
prob_bvp_linear_6_f(t, u2, λ) = -t * u2 - λ * π^2 * cos(π * t) - π * t * sin(π * t)
function prob_bvp_linear_6_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_6_f(t, u[2], p)
end
function prob_bvp_linear_6_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] + 2
end
function prob_bvp_linear_6_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1]
end
prob_bvp_linear_6_function = BVPFunction(
    prob_bvp_linear_6_f!, (prob_bvp_linear_6_bca!, prob_bvp_linear_6_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_6_analytic, twopoint = Val(true))
prob_bvp_linear_6_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_6

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(t, y_2)
```

where

```math
f(t, y_2)=ty_2 - \lambda\pi^2\cos(\pi t)-\pi t\sin(\pi t)
```

with boundary condition

```math
y_1(-1)=-2, y_1(1)=0
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = \cos(\pi t) + \erf(t/\sqrt{2\lambda})/\erf(1/\sqrt{2\lambda})
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=176)
"""
prob_bvp_linear_6 = BVProblem(prob_bvp_linear_6_function,
    [1.0, 0.0],
    prob_bvp_linear_6_tspan,
    λ)

################### linear_bvp7 ############################
function prob_bvp_linear_7_analytic(u, λ, t)
    [
        cos(π * t) + t +
        (t * erf(t / sqrt(2 * λ)) + sqrt(2 * λ / π) * exp(-t^2 / (2 * λ))) /
        (erf(1 / sqrt(2 * λ)) + sqrt(2 * λ / π) * exp(-1 / (2 * λ))),
        -π * sin(π * t) + 1 +
        (erf(t / sqrt(2 * λ)) + 2 * t / sqrt(π) * exp(-t^2 / (2 * λ)) / sqrt(2 * λ) -
         t / λ * sqrt(2 * λ / π) * exp(-t^2 / (2 * λ))) /
        (erf(1 / sqrt(2 * λ)) + sqrt(2 * λ / π) * exp(-1 / (2 * λ)))]
end
function prob_bvp_linear_7_f(t, u1, u2, λ)
    -t * u2 + u1 - (1 + λ * π^2) * cos(π * t) - π * t * sin(π * t)
end
function prob_bvp_linear_7_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_7_f(t, u[1], u[2], p)
end
function prob_bvp_linear_7_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] + 1
end
function prob_bvp_linear_7_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1
end
prob_bvp_linear_7_function = BVPFunction(
    prob_bvp_linear_7_f!, (prob_bvp_linear_7_bca!, prob_bvp_linear_7_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_7_analytic, twopoint = Val(true))
prob_bvp_linear_7_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_7

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(t, y_1, y_2)
```

where

```math
f(t, y_1, y_2)=ty_2+y_1-(1+\lambda\pi^2)\cos(\pi t)+\pi t\sin(\pi t)
```

with boundary condition

```math
y_1(-1)=-1, y_1(1)=1
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = \cos(\pi t) + t + \frac{t\erf(t/\sqrt{2\lambda}) + \sqrt{2\lambda/\pi}\exp(-t^2/2\lambda)}{}
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=178)
"""
prob_bvp_linear_7 = BVProblem(prob_bvp_linear_7_function,
    [1.0, 0.0],
    prob_bvp_linear_7_tspan,
    λ)

################### linear_bvp8 ############################
function prob_bvp_linear_8_analytic(u, λ, t)
    a = 1 / λ
    [(2 - exp(-a) - exp(-t * a)) / (1 - exp(-a)),
        a * exp(-t * a) / (1 - exp(-a))]
end
function prob_bvp_linear_8_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = -1 / p * u[2]
end
function prob_bvp_linear_8_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1
end
function prob_bvp_linear_8_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 2
end
prob_bvp_linear_8_function = BVPFunction(
    prob_bvp_linear_8_f!, (prob_bvp_linear_8_bca!, prob_bvp_linear_8_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_8_analytic, twopoint = Val(true))
prob_bvp_linear_8_tspan = (0.0, 1.0)
@doc raw"""
    prob_bvp_linear_8

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = -\frac{1}{\lambda}y_2
```

with boundary condition

```math
y_1(0)=1, y_1(1)=2
```

# Solution

The analytical solution for ``t \in [0, 1]`` is

```math
y_1(t) = (2-\exp(-1/\lambda)-\exp(-t/\lambda))/(1-\exp(-1/\lambda))
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=181)
"""
prob_bvp_linear_8 = BVProblem(prob_bvp_linear_8_function,
    [1.0, 0.0],
    prob_bvp_linear_8_tspan,
    λ)

################### linear_bvp9 ############################
function prob_bvp_linear_9_analytic(u, λ, t)
    [1 / (λ + t^2),
        -2 * t / (λ + t^2)^2]
end
prob_bvp_linear_9_f(t, u1, u2, λ) = -4 * t * u2 - 2 * u1
function prob_bvp_linear_9_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / (p + t^2) * prob_bvp_linear_9_f(t, u[1], u[2], p)
end
function prob_bvp_linear_9_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1 / (1 + p)
end
function prob_bvp_linear_9_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1 / (1 + p)
end
prob_bvp_linear_9_function = BVPFunction(
    prob_bvp_linear_9_f!, (prob_bvp_linear_9_bca!, prob_bvp_linear_9_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_9_analytic, twopoint = Val(true))
prob_bvp_linear_9_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_9

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda+t^2}f(t, y_1, y_2)
```

where

```math
f(t, y_1, y_2)=-4ty_2 - 2y_1
```

with boundary condition

```math
y_1(-1)=1/(1+\lambda), y_1(1)=1/(1+\lambda)
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = 1/(\lambda+t^2)
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=184)
"""
prob_bvp_linear_9 = BVProblem(prob_bvp_linear_9_function,
    [1.0, 0.0],
    prob_bvp_linear_9_tspan,
    λ)

################### linear_bvp10 ############################
function prob_bvp_linear_10_analytic(u, λ, t)
    a = 1 / sqrt(2 * λ)
    [1 + erf(t * a) / erf(a),
        a / erf(a) * 2 / sqrt(π) * exp(-t^2 / (2 * λ))]
end
prob_bvp_linear_10_f(t, u2, λ) = -t * u2
function prob_bvp_linear_10_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_10_f(t, u[2], p)
end
function prob_bvp_linear_10_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
end
function prob_bvp_linear_10_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 2
end
prob_bvp_linear_10_function = BVPFunction(
    prob_bvp_linear_10_f!, (prob_bvp_linear_10_bca!, prob_bvp_linear_10_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_10_analytic, twopoint = Val(true))
prob_bvp_linear_10_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_10

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(t, y_2)
```

where

```math
f(t, y_2)=-ty_2
```

with boundary condition

```math
y_1(-1)=0, y_1(1)=2
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = 1+\erf(t/\sqrt{2\lambda})/\erf(1/\sqrt{2\lambda})
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=186)
"""
prob_bvp_linear_10 = BVProblem(prob_bvp_linear_10_function,
    [1.0, 0.0],
    prob_bvp_linear_10_tspan,
    λ)

################### linear_bvp11 ############################
function prob_bvp_linear_11_analytic(u, λ, t)
    [cos(π * t),
        -π * sin(π * t)]
end
prob_bvp_linear_11_f(t, u1, λ) = u1 - (1 + λ * π^2) * cos(π * t)
function prob_bvp_linear_11_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_11_f(t, u[1], p)
end
function prob_bvp_linear_11_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] + 1
end
function prob_bvp_linear_11_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] + 1
end
prob_bvp_linear_11_function = BVPFunction(
    prob_bvp_linear_11_f!, (prob_bvp_linear_11_bca!, prob_bvp_linear_11_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_11_analytic, twopoint = Val(true))
prob_bvp_linear_11_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_11

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(t, y_1)
```

where

```math
f(t, y_1)=y_1-(1+\lambda\pi^2)\cos(\pi t)
```

with boundary condition

```math
y_1(-1)=0, y_1(1)=2
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = \cos(\pi t)
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=188)
"""
prob_bvp_linear_11 = BVProblem(prob_bvp_linear_11_function,
    [1.0, 0.0],
    prob_bvp_linear_11_tspan,
    λ)

################### linear_bvp12 ############################
function prob_bvp_linear_12_analytic(u, λ, t)
    a = 1 / sqrt(λ)
    [cos(π * t) + (exp((t + 1) * a) - exp((-t - 1) * a)) / (exp(2 * a) - exp(-2 * a)),
        -π * sin(π * t) +
        (a * exp((t + 1) * a) + a * exp((-t - 1) * a)) / (exp(2 * a) - exp(-2 * a))]
end
prob_bvp_linear_12_f(t, u1, λ) = u1 - (1 + λ * π^2) * cos(π * t)
function prob_bvp_linear_12_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_12_f(t, u[1], p)
end
function prob_bvp_linear_12_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] + 1
end
function prob_bvp_linear_12_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1]
end
prob_bvp_linear_12_function = BVPFunction(
    prob_bvp_linear_12_f!, (prob_bvp_linear_12_bca!, prob_bvp_linear_12_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_12_analytic, twopoint = Val(true))
prob_bvp_linear_12_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_12

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(t, y_1)
```

where

```math
f(t, y_1)=y_1-(1+\lambda\pi^2)\cos(\pi t)
```

with boundary condition

```math
y_1(-1)=-1, y_1(1)=0
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = \cos(\pi t)+\frac{\exp((t+1)/\sqrt{\lambda})-\exp((-t-1))/\sqrt{\lambda}}{\exp(2/\sqrt{\lambda})-\exp(-2/\sqrt{\lambda})}
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=209)
"""
prob_bvp_linear_12 = BVProblem(prob_bvp_linear_12_function,
    [1.0, 0.0],
    prob_bvp_linear_12_tspan,
    λ)

################### linear_bvp13 ############################
function prob_bvp_linear_13_analytic(u, λ, t)
    a = 1 / sqrt(λ)
    [cos(π * t) + exp((-t - 1) * a),
        -π * sin(π * t) - a * exp((-t - 1) * a)]
end
prob_bvp_linear_13_f(t, u1, λ) = u1 - (1 + λ * π^2) * cos(π * t)
function prob_bvp_linear_13_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_13_f(t, u[1], p)
end
function prob_bvp_linear_13_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
end
function prob_bvp_linear_13_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] + 1 - exp(-2 / sqrt(p))
end
prob_bvp_linear_13_function = BVPFunction(
    prob_bvp_linear_13_f!, (prob_bvp_linear_13_bca!, prob_bvp_linear_13_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_13_analytic, twopoint = Val(true))
prob_bvp_linear_13_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_13

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(t, y_1)
```

where

```math
f(t, y_1)=y_1-(1+\lambda\pi^2)\cos(\pi t)
```

with boundary condition

```math
y_1(-1)=0, y_1(1)=-1
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = \cos(\pi t)+\exp(-(t+1)/\sqrt{\lambda})
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=212)
"""
prob_bvp_linear_13 = BVProblem(prob_bvp_linear_13_function,
    [1.0, 0.0],
    prob_bvp_linear_13_tspan,
    λ)

################### linear_bvp14 ############################
function prob_bvp_linear_14_analytic(u, λ, t)
    a = 1 / sqrt(λ)
    [cos(π * t) + exp((t - 1) * a) + exp((-t - 1) * a),
        -π * sin(π * t) + a * exp((t - 1) * a) - a * exp((-t - 1) * a)]
end
prob_bvp_linear_14_f(t, u1, λ) = u1 - (1 + λ * π^2) * cos(π * t)
function prob_bvp_linear_14_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_14_f(t, u[1], p)
end
function prob_bvp_linear_14_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - exp(-2 / sqrt(p))
end
function prob_bvp_linear_14_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - exp(-2 / sqrt(p))
end
prob_bvp_linear_14_function = BVPFunction(
    prob_bvp_linear_14_f!, (prob_bvp_linear_14_bca!, prob_bvp_linear_14_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_14_analytic, twopoint = Val(true))
prob_bvp_linear_14_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_14

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(t, y_1)
```

where

```math
f(t, y_1)=y_1-(1+\lambda\pi^2)\cos(\pi t)
```

with boundary condition

```math
y_1(-1)=\exp(-2/\sqrt{\lambda}, y_1(1)=\exp(-2/\sqrt{\lambda})
```

# Solution

The analytical solution for ``t \in [-1, 1]`` is

```math
y_1(t) = \cos(\pi t)+\exp((t-1)/\sqrt{\lambda})+\exp(-(t+1)/\sqrt{\lambda})
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=215)
"""
prob_bvp_linear_14 = BVProblem(prob_bvp_linear_14_function,
    [1.0, 0.0],
    prob_bvp_linear_14_tspan,
    λ)

################### linear_bvp15 ############################
# No analytical solution
prob_bvp_linear_15_f(t, u1, λ) = t * u1
function prob_bvp_linear_15_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_linear_15_f(t, u[1], p)
end
function prob_bvp_linear_15_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1
end
function prob_bvp_linear_15_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1
end
prob_bvp_linear_15_function = BVPFunction(
    prob_bvp_linear_15_f!, (prob_bvp_linear_15_bca!, prob_bvp_linear_15_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_linear_15_tspan = (-1.0, 1.0)
@doc raw"""
    prob_bvp_linear_15

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(t, y_1)
```

where

```math
f(t, y_1)=ty_1
```

with boundary condition

```math
y_1(-1)=1, y_1(1)=1
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=218)
"""
prob_bvp_linear_15 = BVProblem(prob_bvp_linear_15_function,
    [1.0, 0.0],
    prob_bvp_linear_15_tspan,
    λ)

################### linear_bvp16 ############################
function prob_bvp_linear_16_analytic(u, λ, t)
    a = π / (2 * λ)
    [sin(a * t),
        a * cos(a * t)]
end
function prob_bvp_linear_16_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = -1 / p^2 * π^2 * u[1] / 4
end
function prob_bvp_linear_16_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
end
function prob_bvp_linear_16_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - sin(π / (2 * p))
end
prob_bvp_linear_16_function = BVPFunction(
    prob_bvp_linear_16_f!, (prob_bvp_linear_16_bca!, prob_bvp_linear_16_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_16_analytic, twopoint = Val(true))
prob_bvp_linear_16_tspan = (0.0, 1.0)
@doc raw"""
    prob_bvp_linear_16

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda^2}f(y_1)
```

where

```math
f(t, y_1)=-π^2y_1/4
```

with boundary condition

```math
y_1(0)=0, y_1(1)=\sin(\pi/(2*\lambda))
```

# Solution

The analytical solution for ``t \in [0, 1]`` is

```math
y_1(t) = \sin(\pi t/2\lambda) when 1/\lambda is odd
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=221)
"""
prob_bvp_linear_16 = BVProblem(prob_bvp_linear_16_function,
    [1.0, 0.0],
    prob_bvp_linear_16_tspan,
    λ)

################### linear_bvp17 ############################
function prob_bvp_linear_17_analytic(u, λ, t)
    [t / sqrt(λ + t^2),
        1 / sqrt(λ + t^2) - t^2 / (λ + t^2)^(3 / 2)]
end
function prob_bvp_linear_17_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = -3 * p * u[1] / (p + t^2)^2
end
function prob_bvp_linear_17_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] + 0.1 / sqrt(p + 0.01)
end
function prob_bvp_linear_17_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 0.1 / sqrt(p + 0.01)
end
prob_bvp_linear_17_function = BVPFunction(
    prob_bvp_linear_17_f!, (prob_bvp_linear_17_bca!, prob_bvp_linear_17_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_17_analytic, twopoint = Val(true))
prob_bvp_linear_17_tspan = (-0.1, 0.1)
@doc raw"""
    prob_bvp_linear_17

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = f(y_1)
```

where

```math
f(t, y_1)=-3\lambda y_1/(\lambda+t^2)^2
```

with boundary condition

```math
y_1(-0.1)=-0.1/\sqrt{\lambda+0.01}, y_1(0.1)=0.1/\sqrt{\lambda+0.01}
```

# Solution

The analytical solution for ``t \in [-0.1, 0.1]`` is

```math
y_1(t) = t/\sqrt{\lambda+t^2}
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=224)
"""
prob_bvp_linear_17 = BVProblem(prob_bvp_linear_17_function,
    [1.0, 0.0],
    prob_bvp_linear_17_tspan,
    λ)

################### linear_bvp18 ############################
function prob_bvp_linear_18_analytic(u, λ, t)
    a = -1 / λ
    [exp(a * t),
        a * exp(a * t)]
end
function prob_bvp_linear_18_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = -1 / p * u[2]
end
function prob_bvp_linear_18_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1
end
function prob_bvp_linear_18_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - exp(-1 / p)
end
prob_bvp_linear_18_function = BVPFunction(
    prob_bvp_linear_18_f!, (prob_bvp_linear_18_bca!, prob_bvp_linear_18_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_linear_18_analytic, twopoint = Val(true))
prob_bvp_linear_18_tspan = (0, 1)
@doc raw"""
    prob_bvp_linear_18

Linear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = -\frac{1}{\lambda}y_2
```

where

```math
f(y_2)=-y_1
```

with boundary condition

```math
y_1(0)=1, y_1(1)=0.1/\sqrt{\lambda+0.01}
```

# Solution

The analytical solution for ``t \in [0, 1]`` is

```math
y_1(t) = \exp(-t/\lambda)
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=227)
"""
prob_bvp_linear_18 = BVProblem(prob_bvp_linear_18_function,
    [1.0, 0.0],
    prob_bvp_linear_18_tspan,
    λ)

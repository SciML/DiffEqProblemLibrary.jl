λ = 1

################### nonlinear_bvp1 ############################
# No analytical solution
prob_bvp_nonlinear_1_f(t, u1, u2) = -exp(u1) * u2 + π / 2 * sin(π * t / 2) * exp(2 * u1)
function prob_bvp_nonlinear_1_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_nonlinear_1_f(t, u[1], u[2])
end
function prob_bvp_nonlinear_1_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1
end
function prob_bvp_nonlinear_1_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - exp(-1 / p)
end
prob_bvp_nonlinear_1_function = BVPFunction(
    prob_bvp_nonlinear_1_f!, (prob_bvp_nonlinear_1_bca!, prob_bvp_nonlinear_1_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_nonlinear_1_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_1

Nonlinear boundary value problem with no analytical solution, given by

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

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=230)
"""
prob_bvp_nonlinear_1 = BVProblem(prob_bvp_nonlinear_1_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_1_tspan,
    λ)

################### nonlinear_bvp2 ############################
function prob_bvp_nonlinear_2_analytic(u, λ, t)
    [1 + λ * log(cosh((t - 0.745) / λ)),
        tanh((t - 0.745) / λ) / λ]
end
prob_bvp_nonlinear_2_f(u2) = -u2^2 + 1
function prob_bvp_nonlinear_2_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_nonlinear_2_f(u[2])
end
function prob_bvp_nonlinear_2_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1 - p * log(cosh(-0.745 / p))
end
function prob_bvp_nonlinear_2_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1 - p * log(cosh(0.255 / p))
end
prob_bvp_nonlinear_2_function = BVPFunction(
    prob_bvp_nonlinear_2_f!, (prob_bvp_nonlinear_2_bca!, prob_bvp_nonlinear_2_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_nonlinear_2_analytic, twopoint = Val(true))
prob_bvp_nonlinear_2_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_2

Nonlinear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = -\frac{1}{\lambda}f(y_2)
```

where

```math
f(y_2)=--y_2^2+1
```

with boundary condition

```math
y_1(0)=1+\lambda\ln\cosh(-0.745/\lambda), y_1(1)=1+\lambda\ln\cosh(0.255/\lambda)
```

# Solution

The analytical solution for ``t \in [0, 1]`` is

```math
y_1(t) = 1+\lambda\ln\cosh((t-0.745)/\lambda)
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=234)
"""
prob_bvp_nonlinear_2 = BVProblem(prob_bvp_nonlinear_2_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_2_tspan,
    λ)

################### nonlinear_bvp3 ############################
function prob_bvp_nonlinear_3_analytic(u, λ, t)
    a = -1 / sqrt(λ)
    [exp(a * t),
        a * exp(a * t)]
end
prob_bvp_nonlinear_3_f(t, u1, p) = u1 + u1^2 - exp(-2 * t / sqrt(p))
function prob_bvp_nonlinear_3_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_nonlinear_3_f(t, u[1], p)
end
function prob_bvp_nonlinear_3_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1
end
function prob_bvp_nonlinear_3_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - exp(-1 / sqrt(λ))
end
prob_bvp_nonlinear_3_function = BVPFunction(
    prob_bvp_nonlinear_3_f!, (prob_bvp_nonlinear_3_bca!, prob_bvp_nonlinear_3_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), analytic = prob_bvp_nonlinear_3_analytic, twopoint = Val(true))
prob_bvp_nonlinear_3_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_3

Nonlinear boundary value problem with analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = -\frac{1}{\lambda}f(y, y_1)
```

where

```math
f(y_1)=y_1+y_1^2-\exp(-2t/\sqrt{\lambda})
```

with boundary condition

```math
y_1(0)=1, y_1(1)=\exp(-1/\sqrt{\lambda})
```

# Solution

The analytical solution for ``t \in [0, 1]`` is

```math
y_1(t) = \exp(-t/\sqrt{\lambda})
```
```math
y_2(t) = y_1'(t)
```

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=237)
"""
prob_bvp_nonlinear_3 = BVProblem(prob_bvp_nonlinear_3_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_3_tspan,
    λ)

################### nonlinear_bvp4 ############################
# No analytical solution
prob_bvp_nonlinear_4_f(u1, u2) = -u2 - u1^2
function prob_bvp_nonlinear_4_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_nonlinear_4_f(u[1], u[2])
end
function prob_bvp_nonlinear_4_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
end
function prob_bvp_nonlinear_4_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1 / 2
end
prob_bvp_nonlinear_4_function = BVPFunction(
    prob_bvp_nonlinear_4_f!, (prob_bvp_nonlinear_4_bca!, prob_bvp_nonlinear_4_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_nonlinear_4_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_4

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = -\frac{1}{\lambda}f(y_1, y_2)
```

where

```math
f(y_1, y_2)=-y_2-y_1^2
```

with boundary condition

```math
y_1(0)=0, y_1(1)=1/2
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=240)
"""
prob_bvp_nonlinear_4 = BVProblem(prob_bvp_nonlinear_4_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_4_tspan,
    λ)

################### nonlinear_bvp5 ############################
# No analytical solution
prob_bvp_nonlinear_5_f(u1, p) = p * sinh(p * u1)
function prob_bvp_nonlinear_5_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = prob_bvp_nonlinear_5_f(u[1], p)
end
function prob_bvp_nonlinear_5_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
end
function prob_bvp_nonlinear_5_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1
end
prob_bvp_nonlinear_5_function = BVPFunction(
    prob_bvp_nonlinear_5_f!, (prob_bvp_nonlinear_5_bca!, prob_bvp_nonlinear_5_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_nonlinear_5_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_5

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = -\frac{1}{\lambda}f(y_1)
```

where

```math
f(y_1)=\lambda\sinh(\lambda z)
```

with boundary condition

```math
y_1(0)=0, y_1(1)=1
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=243)
"""
prob_bvp_nonlinear_5 = BVProblem(prob_bvp_nonlinear_5_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_5_tspan,
    λ)

################### nonlinear_bvp6 ############################
# No analytical solution
γ = 1.4
A(t) = 1 + t^2
function prob_bvp_nonlinear_6_f(t, u1, u2, p)
    ((1 + p) / 2 - p * 2 * t) * u1 * u2 - u2 / u1 - 2 * t / A(t) * (1 - (p - 1) / 2 * u1^2)
end
function prob_bvp_nonlinear_6_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = prob_bvp_nonlinear_6_f(t, u[1], u[2], p) / (p * A(t) * u[1])
end
function prob_bvp_nonlinear_6_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 0.9129
end
function prob_bvp_nonlinear_6_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 0.375
end
prob_bvp_nonlinear_6_function = BVPFunction(
    prob_bvp_nonlinear_6_f!, (prob_bvp_nonlinear_6_bca!, prob_bvp_nonlinear_6_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_nonlinear_6_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_6

This problem describes a shock wave in a one dimension nozzle flow.

The steady state Navier-Stokes equations generate a second order differential equations which can be reduced to a first order system described by nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = -\frac{1}{\lambda}f(y_1)
```

where

```math
f(t, y_1, y_2)=(\frac{1+\gamma}{2}-\lambda A'(t))y_1y_2-\frac{y_2}{y_1}-\frac{A'(t)}{A(t)}(1-(\frac{\gamma-1}{2})y_1^2)
```

with boundary condition

```math
y_1(0)=0.9129, y_1(1)=0.375
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=247)
"""
prob_bvp_nonlinear_6 = BVProblem(prob_bvp_nonlinear_6_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_6_tspan,
    λ)

################### nonlinear_bvp7 ############################
# No analytical solution
prob_bvp_nonlinear_7_f(u1, u2) = -u1 * u2 + u1
function prob_bvp_nonlinear_7_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_nonlinear_7_f(u[1], u[2])
end
function prob_bvp_nonlinear_7_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] + 1 / 3
end
function prob_bvp_nonlinear_7_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1 / 3
end
prob_bvp_nonlinear_7_function = BVPFunction(
    prob_bvp_nonlinear_7_f!, (prob_bvp_nonlinear_7_bca!, prob_bvp_nonlinear_7_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_nonlinear_7_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_7

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(y_1, y_2)
```

where

```math
f(y_1, y_2)=-y_1y_2+y_1
```

with boundary condition

```math
y_1(0)=-1/3, y_1(1)=1/3
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=250)
"""
prob_bvp_nonlinear_7 = BVProblem(prob_bvp_nonlinear_7_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_7_tspan,
    λ)

################### nonlinear_bvp8 ############################
# No analytical solution
prob_bvp_nonlinear_8_f(u1, u2) = -u1 * u2 + u1
function prob_bvp_nonlinear_8_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_nonlinear_8_f(u[1], u[2])
end
function prob_bvp_nonlinear_8_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1
end
function prob_bvp_nonlinear_8_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] + 1 / 3
end
prob_bvp_nonlinear_8_function = BVPFunction(
    prob_bvp_nonlinear_8_f!, (prob_bvp_nonlinear_8_bca!, prob_bvp_nonlinear_8_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_nonlinear_8_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_8

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(y_1, y_2)
```

where

```math
f(y_1, y_2)=-y_1y_2+y_1
```

with boundary condition

```math
y_1(0)=1, y_1(1)=-1/3
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=252)
"""
prob_bvp_nonlinear_8 = BVProblem(prob_bvp_nonlinear_8_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_8_tspan,
    λ)

################### nonlinear_bvp9 ############################
# No analytical solution
prob_bvp_nonlinear_9_f(u1, u2) = -u1 * u2 + u1
function prob_bvp_nonlinear_9_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_nonlinear_9_f(u[1], u[2])
end
function prob_bvp_nonlinear_9_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1
end
function prob_bvp_nonlinear_9_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1 / 3
end
prob_bvp_nonlinear_9_function = BVPFunction(
    prob_bvp_nonlinear_9_f!, (prob_bvp_nonlinear_9_bca!, prob_bvp_nonlinear_9_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_nonlinear_9_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_9

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(y_1, y_2)
```

where

```math
f(y_1, y_2)=-y_1y_2+y_1
```

with boundary condition

```math
y_1(0)=1, y_1(1)=1/3
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=252)
"""
prob_bvp_nonlinear_9 = BVProblem(prob_bvp_nonlinear_9_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_9_tspan,
    λ)

################### nonlinear_bvp10 ############################
# No analytical solution
prob_bvp_nonlinear_10_f(u1, u2) = -u1 * u2 + u1
function prob_bvp_nonlinear_10_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_nonlinear_10_f(u[1], u[2])
end
function prob_bvp_nonlinear_10_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] - 1
end
function prob_bvp_nonlinear_10_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 3 / 2
end
prob_bvp_nonlinear_10_function = BVPFunction(
    prob_bvp_nonlinear_10_f!, (prob_bvp_nonlinear_10_bca!, prob_bvp_nonlinear_10_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_nonlinear_10_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_10

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(y_1, y_2)
```

where

```math
f(y_1, y_2)=-y_1y_2+y_1
```

with boundary condition

```math
y_1(0)=1, y_1(1)=3/2
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=256)
"""
prob_bvp_nonlinear_10 = BVProblem(prob_bvp_nonlinear_10_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_10_tspan,
    λ)

################### nonlinear_bvp11 ############################
# No analytical solution
prob_bvp_nonlinear_11_f(u1, u2) = -u1 * u2 + u1
function prob_bvp_nonlinear_11_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_nonlinear_11_f(u[1], u[2])
end
function prob_bvp_nonlinear_11_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
end
function prob_bvp_nonlinear_11_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 3 / 2
end
prob_bvp_nonlinear_11_function = BVPFunction(
    prob_bvp_nonlinear_11_f!, (prob_bvp_nonlinear_11_bca!, prob_bvp_nonlinear_11_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_nonlinear_11_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_11

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(y_1, y_2)
```

where

```math
f(y_1, y_2)=-y_1y_2+y_1
```

with boundary condition

```math
y_1(0)=0, y_1(1)=3/2
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=258)
"""
prob_bvp_nonlinear_11 = BVProblem(prob_bvp_nonlinear_11_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_11_tspan,
    λ)

################### nonlinear_bvp12 ############################
# No analytical solution
prob_bvp_nonlinear_12_f(u1, u2) = -u1 * u2 + u1
function prob_bvp_nonlinear_12_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * prob_bvp_nonlinear_12_f(u[1], u[2])
end
function prob_bvp_nonlinear_12_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] + 7 / 6
end
function prob_bvp_nonlinear_12_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 3 / 2
end
prob_bvp_nonlinear_12_function = BVPFunction(
    prob_bvp_nonlinear_12_f!, (prob_bvp_nonlinear_12_bca!, prob_bvp_nonlinear_12_bcb!),
    bcresid_prototype = (zeros(1), zeros(1)), twopoint = Val(true))
prob_bvp_nonlinear_12_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_12

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}f(y_1, y_2)
```

where

```math
f(y_1, y_2)=-y_1y_2+y_1
```

with boundary condition

```math
y_1(0)=-7/6, y_1(1)=3/2
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=260)
"""
prob_bvp_nonlinear_12 = BVProblem(prob_bvp_nonlinear_12_function,
    [1.0, 0.0],
    prob_bvp_nonlinear_12_tspan,
    λ)

################### nonlinear_bvp13 ############################
# No analytical solution
function prob_bvp_nonlinear_13_f(z, θ, M, Q, p)
    1 / p * ((z - 1) * cos(θ) - M * sec(θ)) + p * Q * tan(θ)
end
function prob_bvp_nonlinear_13_f!(du, u, p, t)
    du[1] = sin(u[2])
    du[2] = u[3]
    du[3] = -u[4] / p
    du[4] = prob_bvp_nonlinear_13_f(u[1], u[2], u[3], u[4], p)
end
function prob_bvp_nonlinear_13_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
    res_a[2] = u_a[3]
end
function prob_bvp_nonlinear_13_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1]
    res_b[2] = u_b[3]
end
prob_bvp_nonlinear_13_function = BVPFunction(
    prob_bvp_nonlinear_13_f!, (prob_bvp_nonlinear_13_bca!, prob_bvp_nonlinear_13_bcb!),
    bcresid_prototype = (zeros(2), zeros(2)), twopoint = Val(true))
prob_bvp_nonlinear_13_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_13

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = \sin(y_2)
```
```math
\frac{dy_2}{dt} = y_3
```
```math
\frac{dy_3}{dt} = -y_4/\lambda
```
```math
\frac{dy_4}{dt} = f(y_1, y_2, y_3, y_4)
```

where

```math
f(z, \theta, M, Q)=\frac{1}{\lambda}((z-1)\cos\theta-M\sec\theta)+\lambda Q\tan\theta
```

with boundary condition

```math
y_1(0)=0, y_3(0)=0, y_1(1)=0, y_3(1)=0
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=262)
"""
prob_bvp_nonlinear_13 = BVProblem(prob_bvp_nonlinear_13_function,
    [1.0, 0.0, 0.0, 0.0],
    prob_bvp_nonlinear_13_tspan,
    λ)

################### nonlinear_bvp14 ############################
# No analytical solution
prob_bvp_nonlinear_14_f(y1, y2, y3, y4, p) = p * (y2 * y3 - y1 * y4)
function prob_bvp_nonlinear_14_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = u[3]
    du[3] = u[4]
    du[4] = prob_bvp_nonlinear_14_f(u[1], u[2], u[3], u[4], p)
end
function prob_bvp_nonlinear_14_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
    res_a[2] = u_a[2]
end
function prob_bvp_nonlinear_14_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1
    res_b[2] = u_b[2]
end
prob_bvp_nonlinear_14_function = BVPFunction(
    prob_bvp_nonlinear_14_f!, (prob_bvp_nonlinear_14_bca!, prob_bvp_nonlinear_14_bcb!),
    bcresid_prototype = (zeros(2), zeros(2)), twopoint = Val(true))
prob_bvp_nonlinear_14_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_14

This problem arises from fluid injection through one side of a long vertical channel

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = y_3
```
```math
\frac{dy_3}{dt} = y_4
```
```math
\frac{dy_4}{dt} = f(y_1, y_2, y_3, y_4)
```

where

```math
f(y_1, y_2, y_3, y_4)=\lambda(y_2y_3-y_1y_4)
```

with boundary condition

```math
y_1(0)=0, y_2(0)=0, y_1(1)=1, y_2(1)=0
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=264)
"""
prob_bvp_nonlinear_14 = BVProblem(prob_bvp_nonlinear_14_function,
    [1.0, 0.0, 0.0, 0.0],
    prob_bvp_nonlinear_14_tspan,
    λ)

################### nonlinear_bvp15 ############################
# No analytical solution
function prob_bvp_nonlinear_15_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = 1 / p * u[1] * u[4] - u[3] * u[2]
    du[3] = u[4]
    du[4] = u[5]
    du[5] = u[6]
    du[6] = 1 / p * (-u[3] * u[6] - u[1] * u[2])
end
function prob_bvp_nonlinear_15_bca!(res_a, u_a, p)
    res_a[1] = u_a[1] + 1
    res_a[2] = u_a[3]
    res_a[3] = u_a[4]
end
function prob_bvp_nonlinear_15_bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 1
    res_b[2] = u_b[3]
    res_b[3] = u_b[4]
end
prob_bvp_nonlinear_15_function = BVPFunction(
    prob_bvp_nonlinear_15_f!, (prob_bvp_nonlinear_15_bca!, prob_bvp_nonlinear_15_bcb!),
    bcresid_prototype = (zeros(3), zeros(3)), twopoint = Val(true))
prob_bvp_nonlinear_15_tspan = (0, 1)
@doc raw"""
    prob_bvp_nonlinear_15

This problem arises from fluid injection through one side of a long vertical channel

Nonlinear boundary value problem with no analytical solution, given by

```math
\frac{dy_1}{dt} = y_2
```
```math
\frac{dy_2}{dt} = \frac{1}{\lambda}y_1y_4-y_3y_2
```
```math
\frac{dy_3}{dt} = y_4
```
```math
\frac{dy_4}{dt} = y_5
```
```math
\frac{dy_5}{dt} = y_6
```
```math
\frac{dy_6}{dt} = \frac{1}{\lambda}(-y_3y_6-y_1y_2)
```

with boundary condition

```math
y_1(0)=-1, y_3(0)=0, y_4(0)=0, y_1(1)=1, y_3(1)=0, y_4(1)=0
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=266)
"""
prob_bvp_nonlinear_15 = BVProblem(prob_bvp_nonlinear_15_function,
    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    prob_bvp_nonlinear_15_tspan,
    λ)

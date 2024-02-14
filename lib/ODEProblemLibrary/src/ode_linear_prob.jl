# Linear ODE
linear = (u, p, t) -> (p * u)
linear_analytic = (u0, p, t) -> u0 * exp(p * t)
"""
Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u_0=\\frac{1}{2}``, ``α=1.01``, and solution

```math
u(t) = u_0e^{αt}
```

with Float64s. The parameter is ``α``
"""
prob_ode_linear = ODEProblem(ODEFunction(linear, analytic = linear_analytic),
    1 / 2, (0.0, 1.0), 1.01)

"""
Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u_0=\\frac{1}{2}``, ``α=1.01``, and solution

```math
u(t) = u_0e^{αt}
```

with BigFloats
"""
prob_ode_bigfloatlinear = ODEProblem(ODEFunction(linear, analytic = linear_analytic),
    big(0.5), (0.0, 1.0), big(1.01))

f_2dlinear = (du, u, p, t) -> (@. du = p * u)
f_2dlinear_analytic = (u0, p, t) -> @. u0 * exp(p * t)
"""
4x2 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u_0`` as all uniformly distributed random numbers, 
``α=1.01``, and solution

```math
u(t) = u_0e^{αt}
```

with Float64s
"""
prob_ode_2Dlinear = ODEProblem(ODEFunction(f_2dlinear, analytic = f_2dlinear_analytic),
    rand(4, 2), (0.0, 1.0), 1.01)

"""
100x100 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u_0`` as all uniformly distributed random numbers, 
``α=1.01``, and solution

```math
u(t) = u_0e^{αt}
```

with Float64s
"""
prob_ode_large2Dlinear = ODEProblem(
    ODEFunction(f_2dlinear, analytic = f_2dlinear_analytic),
    rand(100, 100), (0.0, 1.0), 1.01)

"""
4x2 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u_0`` as all uniformly distributed random numbers, 
``α=1.01``, and solution

```math
u(t) = u_0e^{αt}
```

with BigFloats
"""
prob_ode_bigfloat2Dlinear = ODEProblem(
    ODEFunction(f_2dlinear,
        analytic = f_2dlinear_analytic),
    BigFloat.(rand(4, 2)) .* ones(4, 2) / 2, (0.0, 1.0),
    big(1.01))

f_2dlinear_notinplace = (u, p, t) -> p * u
"""
4x2 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u_0`` as all uniformly distributed random numbers, 
``α=1.01``, and solution

```math
u(t) = u_0e^{αt}
```

on Float64. Purposefully not in-place as a test.
"""
prob_ode_2Dlinear_notinplace = ODEProblem(
    ODEFunction(f_2dlinear_notinplace,
        analytic = f_2dlinear_analytic),
    rand(4, 2), (0.0, 1.0), 1.01)

# Linear ODE
linear = (u,p,t) -> (1.01*u)
(f::typeof(linear))(::Type{Val{:analytic}},u0,p,t) = u0*exp(1.01*t)
"""
Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u0=1/2``, ``α=1.01``, and solution

```math
u(t) = u0e^{αt}
```

with Float64s
"""
prob_ode_linear = ODEProblem(linear,1/2,(0.0,1.0))

const linear_bigα = parse(BigFloat,"1.01")
f_linearbig = (u,p,t) -> (linear_bigα*u)
(f::typeof(f_linearbig))(::Type{Val{:analytic}},u0,p,t) = u0*exp(linear_bigα*t)
"""
Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u0=1/2``, ``α=1.01``, and solution

```math
u(t) = u0e^{αt}
```

with BigFloats
"""
prob_ode_bigfloatlinear = ODEProblem(f_linearbig,parse(BigFloat,"0.5"),(0.0,1.0))

f_2dlinear = (du,u,p,t) -> begin
  for i in 1:length(u)
    du[i] = 1.01*u[i]
  end
end
(f::typeof(f_2dlinear))(::Type{Val{:analytic}},u0,p,t) = u0*exp.(1.01*t)
"""
4x2 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u0=1/2``, ``α=1.01``, and solution

```math
u(t) = u0e^{αt}
```

with Float64s
"""
prob_ode_2Dlinear = ODEProblem(f_2dlinear,rand(4,2),(0.0,1.0))

"""
100x100 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u0=1/2``, ``α=1.01``, and solution

```math
u(t) = u0e^{αt}
```

with Float64s
"""
prob_ode_large2Dlinear = ODEProblem(f_2dlinear,rand(100,100),(0.0,1.0))

f_2dlinearbig = (du,u,p,t) -> begin
  for i in 1:length(u)
    du[i] = linear_bigα*u[i]
  end
end
(f::typeof(f_2dlinearbig))(::Type{Val{:analytic}},u0,p,t) = u0*exp.(1.01*t)
"""
4x2 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u0=1/2``, ``α=1.01``, and solution

```math
u(t) = u0e^{αt}
```

with BigFloats
"""
prob_ode_bigfloat2Dlinear = ODEProblem(f_2dlinearbig,map(BigFloat,rand(4,2)).*ones(4,2)/2,(0.0,1.0))
f_2dlinear_notinplace = (u,p,t) -> 1.01*u
(f::typeof(f_2dlinear_notinplace))(::Type{Val{:analytic}},u0,p,t) = u0*exp.(1.01*t)
"""
4x2 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u0=1/2``, ``α=1.01``, and solution

```math
u(t) = u0e^{αt}
```

on Float64. Purposefully not in-place as a test.
"""
prob_ode_2Dlinear_notinplace = ODEProblem(f_2dlinear_notinplace,rand(4,2),(0.0,1.0))

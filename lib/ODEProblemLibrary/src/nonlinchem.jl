function nonLinChem(dy, y, p, t)
    dy[1] = -y[1]
    dy[2] = y[1] - (y[2])^2
    dy[3] = (y[2])^2
end
y0 = [1.0; 0.0; 0.0]
tspan = (0.0, 20.0)
function nlc_analytic(u0, p, t)
    [exp(-t);
     (2sqrt(exp(-t))besselk(1, 2sqrt(exp(-t))) -
      2besselk(1, 2) / besseli(1, 2) * sqrt(exp(-t))besseli(1, 2sqrt(exp(-t)))) /
     (2besselk(0, 2sqrt(exp(-t))) +
      (2besselk(1, 2) / besseli(1, 2))besseli(0, 2sqrt(exp(-t))))
     -exp(-t) + 1 +
     (-2sqrt(exp(-t)) * besselk(1, 2sqrt(exp(-t))) +
      sqrt(exp(-t)) * besseli(1, 2sqrt(exp(-t))) * 2besselk(1, 2) / besseli(1, 2)) /
     (2besselk(0, 2sqrt(exp(-t))) +
      2besselk(1, 2) / besseli(1, 2) * besseli(0, 2sqrt(exp(-t))))]
end
nonLinChem_f = ODEFunction(nonLinChem, analytic = nlc_analytic)

@doc doc"""
Nonlinear system of reactions with an analytical solution

```math
\frac{dy_1}{dt} = -y_1
```

```math
\frac{dy_2}{dt} = y_1 - y_2^2
```

```math
\frac{dy_3}{dt} = y_2^2
```

with initial condition ``y=[1;0;0]`` on a time span of ``t \in (0,20)``

From

Liu, L. C., Tian, B., Xue, Y. S., Wang, M., & Liu, W. J. (2012). Analytic solution 
for a nonlinear chemistry system of ordinary differential equations. Nonlinear 
Dynamics, 68(1-2), 17-21.

The analytical solution is implemented, allowing easy testing of ODE solvers.
"""
prob_ode_nonlinchem = ODEProblem(nonLinChem, y0, tspan)

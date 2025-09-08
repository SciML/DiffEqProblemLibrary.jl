# Strange attractors from https://www.dynamicmath.xyz/strange-attractors/ and
# https://en.wikipedia.org/wiki/List_of_chaotic_maps
# Opted for the equations as reported in papers

## Thomas

function thomas_eqs(du, u, p, t)
    x, y, z = u
    b = p[1]
    du[1] = sin(y) - b * x
    du[2] = sin(z) - b * y
    du[3] = sin(x) - b * z
end

@doc doc"""
Thomas' cyclically symmetric attractor equations

```math
\begin{align*}
\frac{dx}{dt} &= \sin(y) - bx \\
\frac{dy}{dt} &= \sin(z) - by \\
\frac{dz}{dt} &= \sin(x) - bz
\end{align*}
```

with parameter ``b = 0.208186`` and initial conditions ``x(0)=1, y(0)=0, z(0)=0``

[Reference](https://www.worldscientific.com/doi/abs/10.1142/S0218127499001383)

[Wikipedia](https://en.wikipedia.org/wiki/Thomas%27_cyclically_symmetric_attractor)
"""
thomas = ODEFunction(thomas_eqs)
prob_ode_thomas = ODEProblem(thomas, [1.0, 0.0, 0.0], (0.0, 1.0), [0.208186])

## Lorenz

function lorenz_eqs(du, u, p, t)
    x, y, z = u
    σ, ρ, β = p
    du[1] = σ * (y - x)
    du[2] = x * (ρ - z) - y
    du[3] = x * y - β * z
end

@doc doc"""
Lorenz equations

```math
\begin{align*}
\frac{dx}{dt} &= σ(y - x) \\
\frac{dy}{dt} &= x(ρ - z) - y \\
\frac{dz}{dt} &= xy - βz
\end{align*}
```

with parameters ``σ=10, ρ=28, β=8/3`` and initial conditions ``x(0)=1, y(0)=0, z(0)=0``

[Reference](https://journals.ametsoc.org/view/journals/atsc/20/2/1520-0469_1963_020_0130_dnf_2_0_co_2.xml)

[Wikipedia](https://en.wikipedia.org/wiki/Lorenz_system)
"""
lorenz = ODEFunction(lorenz_eqs)
prob_ode_lorenz = ODEProblem(lorenz, [1.0, 0.0, 0.0], (0.0, 1.0), [10.0, 28.0, 8/3])

## Aizawa

function aizawa_eqs(du, u, p, t)
    x, y, z = u
    a, b, c, d, e, f = p
    du[1] = (z - b) * x - d * y
    du[2] = d * x + (z - b) * y
    du[3] = c + a * z - z^3 / 3 - (x^2 + y^2) * (1 + e * z) + f * z * x^3
end

@doc doc"""
Aizawa equations

```math
\begin{align*}
\frac{dx}{dt} &= (z - b)x - dy \\
\frac{dy}{dt} &= dx + (z - b)y \\
\frac{dz}{dt} &= c + az - \frac{z^3}{3} - (x^2 + y^2)(1 + ez) + fzx^3
\end{align*}
```

with parameters ``a=0.95, b=0.7, c=0.6, d=3.5, e=0.25, f=0.1`` and initial conditions ``x(0)=1, y(0)=0, z(0)=0``

[Reference](https://journals.ametsoc.org/view/journals/atsc/20/2/1520-0469_1963_020_0130_dnf_2_0_co_2.xml)
"""
aizawa = ODEFunction(aizawa_eqs)
prob_ode_aizawa = ODEProblem(aizawa, [1.0, 0.0, 0.0], (0.0, 1.0), [0.95, 0.7, 0.6, 3.5, 0.25, 0.1])

## Dadras

function dadras_eqs(du, u, p, t)
    x, y, z = u
    a, b, c, d, e = p
    du[1] = y - a * x + b * y * z
    du[2] = c * y - x * z + z
    du[3] = d * x * y - e * z
end

@doc doc"""
Dadras equations

```math
\begin{align*}
\frac{dx}{dt} &= y - ax + byz \\
\frac{dy}{dt} &= cy - xz + z \\
\frac{dz}{dt} &= dxy - ez
\end{align*}
```

with parameters ``a=3, b=2.7, c=1.7, d=2, e=9`` and initial conditions ``x(0)=1, y(0)=0, z(0)=0``

[Reference](https://www.sciencedirect.com/science/article/abs/pii/S0375960109009591)
"""
dadras = ODEFunction(dadras_eqs)
prob_ode_dadras = ODEProblem(dadras, [1.0, 0.0, 0.0], (0.0, 1.0), [3.0, 2.7, 1.7, 2.0, 9.0])

## Chen

function chen_eqs(du, u, p, t)
    x, y, z = u
    a, b, c = p
    du[1] = a * (y - x)
    du[2] = (c - a) * x - x * z + c * y
    du[3] = x * y - b * z
end

@doc doc"""
Chen equations

```math
\begin{align*}
\frac{dx}{dt} &= a(y - x) \\
\frac{dy}{dt} &= (c - a)x - xz + cy \\
\frac{dz}{dt} &= xy - bz
\end{align*}
```

with parameters ``a=35, b=3, c=28`` and initial conditions ``x(0)=1, y(0)=0, z(0)=0``

[Reference](https://www.worldscientific.com/doi/abs/10.1142/S0218127499001024)
"""
chen = ODEFunction(chen_eqs)
prob_ode_chen = ODEProblem(chen, [1.0, 0.0, 0.0], (0.0, 1.0), [35.0, 3.0, 28.0])

## Rössler

function rossler_eqs(du, u, p, t)
    x, y, z = u
    a, b, c = p
    du[1] = -(y + z)
    du[2] = x + a * y
    du[3] = b + z * (x - c)
end

@doc doc"""
Rössler equations

```math
\begin{align*}
\frac{dx}{dt} &= -(y + z) \\
\frac{dy}{dt} &= x + ay \\
\frac{dz}{dt} &= b + z(x - c)
\end{align*}
```

with parameters ``a=0.2, b=0.2, c=5.7`` and initial conditions ``x(0)=1, y(0)=0, z(0)=0``

[Reference](https://www.sciencedirect.com/science/article/abs/pii/0375960176901018)
[Wikipedia](https://en.wikipedia.org/wiki/R%C3%B6ssler_attractor)
"""
rossler = ODEFunction(rossler_eqs)
prob_ode_rossler = ODEProblem(rossler, [1.0, 0.0, 0.0], (0.0, 1.0), [0.2, 0.2, 5.7])

## Rabinovich-Fabrikant

function rabinovich_fabrikant_eqs(du, u, p, t)
    x, y, z = u
    a, b = p
    du[1] = y * (z - 1 + x^2) + b * x
    du[2] = x * (3 * z + 1 - x^2) + b * y
    du[3] = -2 * z * (a + x * y)
end

@doc doc"""
Rabinovich-Fabrikant equations

```math
\begin{align*}
\frac{dx}{dt} &= y(z - 1 + x^2) + bx \\
\frac{dy}{dt} &= x(3z + 1 - x^2) + by \\
\frac{dz}{dt} &= -2z(a + xy)
\end{align*}
```

with parameters ``a=0.14, b=0.10`` and initial conditions ``x(0)=1, y(0)=0, z(0)=0``

[Reference](https://en.wikipedia.org/wiki/Rabinovich%E2%80%93Fabrikant_equations)
"""
rabinovich_fabrikant = ODEFunction(rabinovich_fabrikant_eqs)
prob_ode_rabinovich_fabrikant = ODEProblem(rabinovich_fabrikant, [1.0, 0.0, 0.0], (0.0, 1.0), [0.14, 0.10])

## Sprott

function sprott_eqs(du, u, p, t)
    x, y, z = u
    a, b = p
    du[1] = y + a * x * y + x * z
    du[2] = 1 - b * x^2 + y * z
    du[3] = x - x^2 - y^2
end

@doc doc"""
Sprott equations

```math
\begin{align*}
\frac{dx}{dt} &= y + axy + xz \\
\frac{dy}{dt} &= 1 - bx^2 + yz \\
\frac{dz}{dt} &= x - x^2 - y^2
\end{align*}
```

with parameters ``a=2.07, b=1.79`` and initial conditions ``x(0)=1, y(0)=0, z(0)=0``

[Reference](https://sprott.physics.wisc.edu/pubs/paper423.pdf)
"""
sprott = ODEFunction(sprott_eqs)
prob_ode_sprott = ODEProblem(sprott, [1.0, 0.0, 0.0], (0.0, 1.0), [2.07, 1.79])

## Hindmarsh-Rose

function hindmarsh_rose_eqs(du, u, p, t)
    x, y, z = u
    a, b, c, d, r, s, xr, i = p
    du[1] = y - a * x^3 + b * x^2 - z + i
    du[2] = c - d * x^2 - y
    du[3] = r * (s * (x - xr) - z)
end

@doc doc"""
Hindmarsh-Rose equations

```math
\begin{align*}
\frac{dx}{dt} &= y - ax^3 + bx^2 - z + i \\
\frac{dy}{dt} &= c - dx^2 - y \\
\frac{dz}{dt} &= r(s(x - x_r) - z)
\end{align*}
```

with parameters ``a=1, b=3, c=1, d=5, r=1e-2, s=4, x_r=-8/5, i=5`` and initial conditions ``x(0)=1, y(0)=0, z(0)=0``

[Reference](https://en.wikipedia.org/wiki/Hindmarsh%E2%80%93Rose_model)
"""
hindmarsh_rose = ODEFunction(hindmarsh_rose_eqs)
prob_ode_hindmarsh_rose = ODEProblem(hindmarsh_rose, [1.0, 0.0, 0.0], (0.0, 1.0), [1.0, 3.0, 1.0, 5.0, 1e-2, 4.0, -8/5, 5.0])

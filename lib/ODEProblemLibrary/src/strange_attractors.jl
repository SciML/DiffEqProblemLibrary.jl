# Strange attractors from https://www.dynamicmath.xyz/strange-attractors/ and
# https://en.wikipedia.org/wiki/List_of_chaotic_maps
# Opted for the equations as reported in papers

# Thomas
@parameters b = 0.208186
@variables x(t)=1 y(t)=0 z(t)=0

eqs = [D(x) ~ sin(y) - b * x,
    D(y) ~ sin(z) - b * y,
    D(z) ~ sin(x) - b * z]

@mtkbuild thomas = ODESystem(eqs, t)

"""
Thomas' cyclically symmetric attractor equations

```math
$(latexify(thomas))
```

[Reference](https://www.worldscientific.com/doi/abs/10.1142/S0218127499001383)

[Wikipedia](https://en.wikipedia.org/wiki/Thomas%27_cyclically_symmetric_attractor)
"""
prob_ode_thomas = ODEProblem(thomas, [], (0.0, 1.0))

# Lorenz
@parameters σ=10 ρ=28 β=8 / 3
@variables x(t)=1 y(t)=0 z(t)=0

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@mtkbuild lorenz = ODESystem(eqs, t)

"""
Lorenz equations

```math
$(latexify(lorenz))
```

[Reference](https://journals.ametsoc.org/view/journals/atsc/20/2/1520-0469_1963_020_0130_dnf_2_0_co_2.xml)

[Wikipedia](https://en.wikipedia.org/wiki/Lorenz_system)
"""
prob_ode_lorenz = ODEProblem(lorenz, [], (0.0, 1.0))

# Aizawa
@parameters a=0.95 b=0.7 c=0.6 d=3.5 e=0.25 f=0.1
@variables x(t)=1 y(t)=0 z(t)=0

eqs = [D(x) ~ (z - b) * x - d * y,
    D(y) ~ d * x + (z - b) * y,
    D(z) ~ c + a * z - z^3 / 3 - (x^2 + y^2) * (1 + e * z) + f * z * x^3]

@mtkbuild aizawa = ODESystem(eqs, t)

"""
Aizawa equations

```math
$(latexify(aizawa))
```

[Reference](https://journals.ametsoc.org/view/journals/atsc/20/2/1520-0469_1963_020_0130_dnf_2_0_co_2.xml)

"""
prob_ode_aizawa = ODEProblem(aizawa, [], (0.0, 1.0))

# Dadras
@parameters a=3 b=2.7 c=1.7 d=2 e=9
@variables x(t)=1 y(t)=0 z(t)=0

eqs = [D(x) ~ y - a * x + b * y * z,
    D(y) ~ c * y - x * z + z,
    D(z) ~ d * x * y - e * z]

@mtkbuild dadras = ODESystem(eqs, t)

"""
Dadras equations

```math
$(latexify(dadras))
```

[Reference](https://www.sciencedirect.com/science/article/abs/pii/S0375960109009591)

"""
prob_ode_dadras = ODEProblem(dadras, [], (0.0, 1.0))

# chen
@parameters a=35 b=3 c=28
@variables x(t)=1 y(t)=0 z(t)=0

eqs = [D(x) ~ a * (y - x),
    D(y) ~ (c - a) * x - x * z + c * y,
    D(z) ~ x * y - b * z]

@mtkbuild chen = ODESystem(eqs, t)

"""
chen equations

```math
$(latexify(chen))
```

[Reference](https://www.worldscientific.com/doi/abs/10.1142/S0218127499001024)

"""
prob_ode_chen = ODEProblem(chen, [], (0.0, 1.0))

# rossler
@parameters a=0.2 b=0.2 c=5.7
@variables x(t)=1 y(t)=0 z(t)=0

eqs = [D(x) ~ -(y + z),
    D(y) ~ x + a * y,
    D(z) ~ b + z * (x - c)]

@mtkbuild rossler = ODESystem(eqs, t)

"""
rossler equations

```math
$(latexify(rossler))
```

[Reference](https://www.sciencedirect.com/science/article/abs/pii/0375960176901018)
[Wikipedia](https://en.wikipedia.org/wiki/R%C3%B6ssler_attractor)

"""
prob_ode_rossler = ODEProblem(rossler, [], (0.0, 1.0))

# rabinovich_fabrikant
@parameters a=0.14 b=0.10
@variables x(t)=1 y(t)=0 z(t)=0

eqs = [D(x) ~ y * (z - 1 + x^2) + b * x,
    D(y) ~ x * (3 * z + 1 - x^2) + b * y,
    D(z) ~ -2 * z * (a + x * y)]

@mtkbuild rabinovich_fabrikant = ODESystem(eqs, t)

"""
rabinovich_fabrikant equations

```math
$(latexify(rabinovich_fabrikant))
```

[Reference](https://en.wikipedia.org/wiki/Rabinovich%E2%80%93Fabrikant_equations)

"""
prob_ode_rabinovich_fabrikant = ODEProblem(rabinovich_fabrikant, [], (0.0, 1.0))

# sprott
@parameters a=2.07 b=1.79
@variables x(t)=1 y(t)=0 z(t)=0

eqs = [D(x) ~ y + a * x * y + x * z,
    D(y) ~ 1 - b * x^2 + y * z,
    D(z) ~ x - x^2 - y^2]

@mtkbuild sprott = ODESystem(eqs, t)

"""
sprott equations

```math
$(latexify(sprott))
```

[Reference](https://sprott.physics.wisc.edu/pubs/paper423.pdf)

"""
prob_ode_sprott = ODEProblem(sprott, [], (0.0, 1.0))

# hindmarsh_rose
@parameters a=1 b=3 c=1 d=5 r=1e-2 s=4 xr=-8 / 5 i=5
@variables x(t)=1 y(t)=0 z(t)=0

eqs = [D(x) ~ y - a * x^3 + b * x^2 - z + i,
    D(y) ~ c - d * x^2 - y,
    D(z) ~ r * (s * (x - xr) - z)]

@mtkbuild hindmarsh_rose = ODESystem(eqs, t)

"""
hindmarsh_rose equations

```math
$(latexify(hindmarsh_rose))
```

[Reference](https://en.wikipedia.org/wiki/Hindmarsh%E2%80%93Rose_model)

"""
prob_ode_hindmarsh_rose = ODEProblem(hindmarsh_rose, [], (0.0, 1.0))

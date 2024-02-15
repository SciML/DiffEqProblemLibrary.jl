#=
The following examples are taken from:
W.H. Enright and H. Hayashi, The evaluation of numerical software for delay
differential equations.
=#

# Problem A1
@doc raw"""
    prob_dde_DDETST_A1

Delay differential equation model of blood production, given by

```math
u'(t) = \frac{0.2 u(t - 14)}{1 + u(t - 14)^{10}} - 0.1 u(t)
```

for ``t \in [0, 500]`` and history function ``\phi(t) = 0.5`` for ``t \leq 0``.

# References

Mackey, M. C. and Glass, L. (1977). Oscillation and chaos in physiological control systems,
Science (197), pp. 287-289.
"""
prob_dde_DDETST_A1

function f_dde_DDETST_A1(u, h, p, t)
    z = h(p, t - 14)

    0.2 * z / (1 + z^10) - 0.1 * u
end

function h_dde_DDETST_A1(p, t)
    t ≤ 0 || error("history function is only implemented for t ≤ 0")

    0.5
end

const prob_dde_DDETST_A1 = DDEProblem(f_dde_DDETST_A1, h_dde_DDETST_A1, (0.0, 500.0);
    constant_lags = [14])

# Problem A2
@doc raw"""
    prob_dde_DDETST_A2

Delay differential equation model of chronic granulocytic leukemia, given by

```math
u_1'(t) = \frac{1.1}{1 + \sqrt{10} u_1(t - 20)^{5/4}} - \frac{10 u_1(t)}{1 + 40 u_2(t)},
```
```math
u_2'(t) = \frac{100 u_1(t)}{1 + 40 u_2(t)} - 2.43 u_2(t),
```

for ``t \in [0, 100]`` and history function

```math
\phi_1(t) = 1.05767027/3,
```
```math
\phi_2(t) = 1.030713491/3,
```

for ``t \leq 0``.

# References

Wheldon, T., Kirk, J. and Finlay, H. (1974). Cyclical granulopoiesis in chronic granulocytic
leukemia: A simulation study., Blood (43), pp. 379-387.
"""
prob_dde_DDETST_A2

function f_dde_DDETST_A2!(du, u, h, p, t)
    z = 10 * u[1] / (1 + 40 * u[2])

    du[1] = 1.1 / (1 + sqrt(10) * (h(p, t - 20; idxs = 1)^(5 / 4))) - z
    du[2] = 10 * z - 2.43 * u[2]

    nothing
end

function h_dde_DDETST_A2(p, t; idxs::Union{Nothing, Int} = nothing)
    t ≤ 0 || error("history function is only implemented for t ≤ 0")

    if idxs === nothing
        [1.05767027 / 3, 1.030713491 / 3]
    elseif idxs == 1
        1.05767027 / 3
    elseif idxs == 2
        1.030713491 / 3
    else
        error("delay differential equation consists of two components")
    end
end

const prob_dde_DDETST_A2 = DDEProblem(f_dde_DDETST_A2!, h_dde_DDETST_A2, (0.0, 100.0);
    constant_lags = [20])

# Problem B2
@doc raw"""
    prob_dde_DDETST_B1

Delay differential equation

```math
u'(t) = 1 - u(\exp(1 - 1/t))
```

for ``t \in [0.1, 10]`` with history function ``\phi(t) = \log t`` for ``t \in (0, 0.1]``.

# Solution

The analytical solution for ``t \in [0.1, 10]`` is

```math
u(t) = \log t.
```

# References

Neves, K. W. (1975). Automatic integration of functional differential equations: An
approach, ACM Trans. Math. Soft. (1), pp. 357-368.
"""
prob_dde_DDETST_B1

function f_dde_DDETST_B1(u, h, p, t)
    1 - h(p, exp(1 - 1 / t))
end

function h_dde_DDETST_B1(p, t)
    0 < t ≤ 0.1 || error("history function is only implemented for t ∈ (0, 0.1]")

    log(t)
end

function fanalytic_dde_DDETST_B1(u₀, ::typeof(h_dde_DDETST_B1), p, t)
    0.1 ≤ t ≤ 10 && u₀ == log(0.1) ||
        error("analytical solution is only implemented for t ∈ [0.1, 10] and u(0) = log 0.1")

    log(t)
end

const prob_dde_DDETST_B1 = DDEProblem(
    DDEFunction(f_dde_DDETST_B1;
        analytic = fanalytic_dde_DDETST_B1),
    h_dde_DDETST_B1, (0.1, 10.0);
    dependent_lags = ((u, p, t) -> t - exp(1 - 1 / t),))

# Problem B2
@doc raw"""
    prob_dde_DDETST_B2

Delay differential equation

```math
u'(t) = - 1 - u(t) + 2 [u(t / 2) < 0]
```

for ``t \in [0, 2 \log 66]`` with history function ``\phi(0) = 1``.

# Solution

The analytical solution for ``t \in [0, 2 \log 66]`` is

```math
u(t) = \begin{cases}
  2 \exp(-t) - 1 & \text{if } t \in [0, 2 \log 2], \\
  1 - 6 \exp(-t) & \text{if  } t \in (2 \log 2, 2 \log 6], \\
  66 \exp(-t) - 1 & \text{if } t \in (2 \log 6, 2 \log 66].
\end{cases}
```

# References

Neves, K. W. and Thompson, S. (1992). Solution of systems of functional differential
equations with state dependent delays, Technical Report TR-92-009, Computer Science,
Radford University.
"""
prob_dde_DDETST_B2

function f_dde_DDETST_B2(u, h, p, t)
    if h(p, t / 2) < 0
        1 - u
    else
        -1 - u
    end
end

function h_dde_DDETST_B2(p, t)
    t == 0 || error("history function is only implemented for t = 0")

    1.0
end

function fanalytic_dde_DDETST_B2(u₀, ::typeof(h_dde_DDETST_B2), p, t)
    0 ≤ t ≤ 2 * log(66) && u₀ == 1 ||
        error("analytical solution is only implemented for t ∈ [0, 2 log(66)] and u(0) = 1")

    if t ≤ 2 * log(2)
        2 * exp(-t) - 1
    elseif t ≤ 2 * log(6)
        1 - 6 * exp(-t)
    else
        66 * exp(-t) - 1
    end
end

const prob_dde_DDETST_B2 = DDEProblem(
    DDEFunction(f_dde_DDETST_B2;
        analytic = fanalytic_dde_DDETST_B2),
    h_dde_DDETST_B2, (0.0, 2 * log(66));
    dependent_lags = ((u, p, t) -> t / 2,))

# Problem C1
@doc raw"""
    prob_dde_DDETST_C1

Delay differential equation

```math
u'(t) = - 2 u(t - 1 - |u(t)|) (1 - u(t)^2)
```

for ``t \in [0, 30]`` with history function ``\phi(t) = 0.5`` for ``t \leq 0``.

# References

Paul, C. A. H. (1994). A test set of functional differential equations, Technical Report
249, The Department of Mathematics, The University of Manchester, Manchester, England.
"""
prob_dde_DDETST_C1

f_dde_DDETST_C1(u, h, p, t) = -2 * h(p, t - 1 - abs(u)) * (1 - u^2)

function h_dde_DDETST_C1(p, t)
    t ≤ 0 || error("history function is only implemented for t ≤ 0")

    0.5
end

const prob_dde_DDETST_C1 = DDEProblem(f_dde_DDETST_C1, h_dde_DDETST_C1, (0.0, 30.0);
    dependent_lags = ((u, p, t) -> 1 + abs(u),))

# Problem C2
@doc raw"""
    prob_dde_DDETST_C2

Delay differential equation

```math
u_1'(t) = - 2 u_1(t - u_2(t)),
```
```math
u_₂'(t) = \frac{|u_1(t - u_2(t))| - |u_1(t)|}{1 + |u_1(t - u_2(t))|},
```

for ``t \in [0, 40]`` with history function

```math
\phi_1(t) = 1,
```
```math
\phi_2(t) = 0.5,
```

for ``t \leq 0``.

# References

Paul, C. A. H. (1994). A test set of functional differential equations, Technical Report
249, The Department of Mathematics, The University of Manchester, Manchester, England.
"""
prob_dde_DDETST_C2

function f_dde_DDETST_C2!(du, u, h, p, t)
    z = h(p, t - u[2]; idxs = 1)
    absz = abs(z)

    du[1] = -2 * z
    du[2] = (absz - abs(u[1])) / (1 + absz)

    nothing
end

function h_dde_DDETST_C2(p, t; idxs::Union{Nothing, Int} = nothing)
    t ≤ 0 || error("history function is only implemented for t ≤ 0")

    if idxs === nothing
        [1.0, 0.5]
    elseif idxs == 1
        1.0
    elseif idxs == 2
        0.5
    else
        error("delay differential equation consists of two components")
    end
end

const prob_dde_DDETST_C2 = DDEProblem(f_dde_DDETST_C2!, h_dde_DDETST_C2, (0.0, 30.0);
    dependent_lags = ((u, p, t) -> u[2],))

# Problem C3
@doc raw"""
    prob_dde_DDETST_C3

Delay differential equation model of hematopoiesis, given by

```math
u_1'(t) = \hat{s}_0 u_2(t - T_1) - \gamma u_1(t) - Q,
```
```math
u_2'(t) = f(u_1(t)) - k u_2(t),
```
```math
u_3'(t) = 1 - \frac{Q \exp(\gamma u_3(t))}{\hat{s}_0 u_2(t - T_1 - u_3(t))},
```

for ``t \in [0, 300]`` with history function ``\phi_1(0) = 3.325``, ``\phi_3(0) = 120``, and

```math
\phi_2(t) = \begin{cases}
  10 & \text{if } t \in [- T_1, 0],\\
  9.5 & \text{if } t < - T_1,
\end{cases}
```

where ``f(y) = a / (1 + K y^r)``, ``\hat{s}_0 = 0.0031``, ``T_1 = 6``, ``\gamma = 0.001``,
``Q = 0.0275``, ``k = 2.8``, ``a = 6570``, ``K = 0.0382``, and ``r = 6.96``.

# References

Mahaffy, J. M., Belair, J. and Mackey, M. C. (1996). Hematopoietic model with moving
boundary condition and state dependent delay, Private communication.
"""
const prob_dde_DDETST_C3 = let s₀ = 0.0031, T₁ = 6, γ = 0.001, Q = 0.0275, k = 2.8,
    a = 6570, K = 0.0382, r = 6.96

    global function f_dde_DDETST_C3!(du, u, h, p, t)
        du[1] = s₀ * h(p, t - T₁; idxs = 2) - γ * u[1] - Q
        du[2] = a / (1 + K * u[1]^r) - k * u[2]
        du[3] = 1 - Q * exp(γ * u[3]) / (s₀ * h(p, t - T₁ - u[3]; idxs = 2))

        nothing
    end

    global function h_dde_DDETST_C3(p, t; idxs::Union{Nothing, Int} = nothing)
        t ≤ 0 || error("history function is only implemented for t ≤ 0")

        if idxs === nothing
            [3.325, t < -T₁ ? 9.5 : 10.0, 120.0]
        elseif idxs == 1
            3.325
        elseif idxs == 2
            t < -T₁ ? 9.5 : 10.0
        elseif idxs == 3
            120.0
        else
            error("delay differential equation consists of three components")
        end
    end

    DDEProblem(f_dde_DDETST_C3!, h_dde_DDETST_C3, (0.0, 300.0);
        constant_lags = [T₁], dependent_lags = ((u, p, t) -> T₁ + u[3],))
end

# Problem C4
@doc raw"""
    prob_dde_DDETST_C4

Delay differential equation model of hematopoiesis, given by the same delay differential
equation as [`prob_dde_DDETST_C3`](@ref)

```math
u_1'(t) = \hat{s}_0 u_2(t - T_1) - \gamma u_1(t) - Q,
```
```math
u_2'(t) = f(u_1(t)) - k u_2(t),
```
```math
u_3'(t) = 1 - \frac{Q \exp(\gamma u_3(t))}{\hat{s}_0 u_2(t - T_1 - u_3(t))},
```

for ``t \in [0, 100]`` with history function
``\phi_1(0) = 3.5``, ``\phi_3(0) = 50``, and ``\phi_2(t) = 10`` for ``t \leq 0``, where
``f(y) = a / (1 + K y^r)``, ``\hat{s}_0 = 0.00372``, ``T_1 = 3``, ``\gamma = 0.1``,
``Q = 0.00178``, ``k = 6.65``, ``a = 15600``, ``K = 0.0382``, and ``r = 6.96``.

# References

Mahaffy, J. M., Belair, J. and Mackey, M. C. (1996). Hematopoietic model with moving
boundary condition and state dependent delay, Private communication.
"""
const prob_dde_DDETST_C4 = let s₀ = 0.00372, T₁ = 3, γ = 0.01, Q = 0.00178, k = 6.65,
    a = 15600, K = 0.0382, r = 6.96

    global function f_dde_DDETST_C4!(du, u, h, p, t)
        du[1] = s₀ * h(p, t - T₁; idxs = 2) - γ * u[1] - Q
        du[2] = a / (1 + K * u[1]^r) - k * u[2]
        du[3] = 1 - Q * exp(γ * u[3]) / (s₀ * h(p, t - T₁ - u[3]; idxs = 2))

        nothing
    end

    global function h_dde_DDETST_C4(p, t; idxs::Union{Nothing, Int} = nothing)
        t ≤ 0 || error("history function is only implemented for t ≤ 0")

        if idxs === nothing
            [3.5, 10.0, 50.0]
        elseif idxs == 1
            3.5
        elseif idxs == 2
            10.0
        elseif idxs == 3
            50.0
        else
            error("delay differential equation consists of three components")
        end
    end

    DDEProblem(f_dde_DDETST_C4!, [3.5, 10.0, 50.0], h_dde_DDETST_C4, (0.0, 100.0);
        constant_lags = [T₁], dependent_lags = ((u, p, t) -> T₁ + u[3],))
end

# Problem D1
@doc raw"""
    prob_dde_DDETST_D1

Delay differential equation

```math
u_1'(t) = u_2(t), 
```
```math
u_2'(t) = - u_2(\exp(1 - u_2(t))) u_2(t)^2 \exp(1 - u_2(t)),
```

for ``t \in [0.1, 5]`` with history function

```math
\phi_1(t) = \log t, 
```
```math
\phi_2(t) = 1 / t,
```

for ``t \in (0, 0.1]``.

# Solution

The analytical solution for ``t \in [0.1, 5]`` is

```math
u_1(t) = \log t, 
```
```math
u_2(t) = 1 / t.
```

# References

Neves, K. W. (1975). Automatic integration of functional differential equations: An
approach, ACM Trans. Math. Soft. (1), pp. 357-368.
"""
prob_dde_DDETST_D1

function f_dde_DDETST_D1!(du, u, h, p, t)
    s = exp(1 - u[2])

    du[1] = u[2]
    du[2] = -h(p, s; idxs = 2) * u[2]^2 * s

    nothing
end

function h_dde_DDETST_D1(p, t; idxs::Union{Nothing, Int} = nothing)
    0 < t ≤ 0.1 || error("history function is only implemented for 0 < t ≤ 0.1")

    if idxs === nothing
        [log(t), 1 / t]
    elseif idxs == 1
        log(t)
    elseif idxs == 2
        1 / t
    else
        error("delay differential equation consists of two components")
    end
end

function fanalytic_dde_DDETST_D1(u₀, ::typeof(h_dde_DDETST_D1), p, t)
    0.1 ≤ t ≤ 5 && u₀ == [log(0.1), 10] ||
        error("analytical solution is only implemented for t ∈ [0.1, 5] and u(0.1) = [log 0.1, 10]")

    [log(t), 1 / t]
end

const prob_dde_DDETST_D1 = DDEProblem(
    DDEFunction(f_dde_DDETST_D1!;
        analytic = fanalytic_dde_DDETST_D1),
    h_dde_DDETST_D1, (0.1, 5.0);
    dependent_lags = ((u, p, t) -> t - exp(1 - u[2]),))

# Problem D2
@doc raw"""
    prob_dde_DDETST_D2

Delay differential equation model of antigen antibody dynamics with fading memory, given by

```math
u_1'(t) = - r_1 u_1(t) u_2(t) + r_2 u_3(t), 
```
```math
u_2'(t) = - r_1 u_1(t) u_2(t) + \alpha r_1 u_1(t - u_4(t)) u_2(t - u_4(t)),
```
```math
u_3'(t) = r_1 u_1(t) u_2(t) - r_2 u_3(t), 
```
```math
u_4'(t) = 1 + \frac{3 \delta - u_1(t) u_2(t) - u_3(t)}{u_1(t - u_4(t)) u_2(t - u_4(t)) + u_3(t - u_4(t))} \exp(\delta u_4(t)),
```

for ``t \in [0, 40]`` with history function

```math
\phi_1(t) = 5, 
```
```math
\phi_2(t) = 0.1, 
```
```math
\phi_3(t) = 0, 
```
```math
\phi_4(t) = 0,
```

for ``t \leq 0``, where ``r_1 = 0.02``, ``r_2 = 0.005``, ``\alpha = 3``, and ``\delta = 0.01``.

# References

Gatica, J. and Waltman, P. (1982). A threshold model of antigen antibody dynamics with
fading memory, in Lakshmikantham (ed.), Nonlinear phenomena in mathematical science,
Academic Press, New York, pp. 425-439.
"""
const prob_dde_DDETST_D2 = let r₁ = 0.02, r₂ = 0.005, α = 3, δ = 0.01
    global function f_dde_DDETST_D2!(du, u, h, p, t)
        u₁u₂ = u[1] * u[2]
        r₁u₁u₂ = r₁ * u₁u₂
        r₂u₃ = r₂ * u[3]

        v = h(p, t - u[4])
        v₁v₂ = v[1] * v[2]

        du[1] = -r₁u₁u₂ + r₂u₃
        du[2] = -r₁u₁u₂ + α * r₁ * v₁v₂
        du[3] = -du[1]
        du[4] = 1 + (3 * δ - u₁u₂ - u[3]) * exp(δ * u[4]) / (v₁v₂ + v[3])

        nothing
    end

    global function h_dde_DDETST_D2(p, t)
        t ≤ 0 || error("history function is only implemented for t ≤ 0")

        [5.0, 0.1, 0.0, 0.0]
    end

    DDEProblem(f_dde_DDETST_D2!, h_dde_DDETST_D2, (0.0, 40.0);
        dependent_lags = ((u, p, t) -> u[4],))
end

# Problem E1
@doc raw"""
    prob_dde_DDETST_E1

Delay differential equation model of a food-limited population, given by

```math
u(t) = r u(t) (1 - u(t - 1) - c u'(t - 1))
```

for ``t \in [0, 40]`` with history function ``\phi(t) = 2 + t`` for ``t \leq 0``, where
``r = \pi / \sqrt{3} + 1/20`` and ``c = \sqrt{3} / (2 \pi) - 1 / 25``.

# References

Kuang, Y. and Feldstein, A. (1991). Boundedness of solutions of a nonlinear nonautonomous
neutral delay equation, J. Math. Anal. Appl. (156), pp. 293-304.
"""
const prob_dde_DDETST_E1 = let r = π / sqrt(3) + 1 / 20, c = sqrt(3) / (2 * π) - 1 / 25
    global f_dde_DDETST_E1(u, h, p, t) = r * u * (1 - h(p, t - 1) - c * h(p, t - 1, Val{1}))

    global function h_dde_DDETST_E1(p, t)
        t ≤ 0 || error("history function is only implemented for t ≤ 0")

        2 + t
    end

    global function h_dde_DDETST_E1(p, t, ::Type{Val{1}})
        t ≤ 0 || error("history function is only implemented for t ≤ 0")

        1.0
    end

    DDEProblem(f_dde_DDETST_E1, h_dde_DDETST_E1, (0.0, 40.0);
        constant_lags = [1], neutral = true)
end

# Problem E2
@doc raw"""
    prob_dde_DDETST_E2

Delay differential equation model of a logistic Gauss-type predator-prey system, given by

```math
u_1'(t) = u_1(t) (1 - u_1(t - \tau) - \rho u_1'(t - \tau)) - \frac{u_2(t) u_1(t)^2}{u_1(t)^2 + 1}, 
```
```math
u_2'(t) = u_2(t) \left(\frac{u_1(t)^2}{u_1(t)^2 + 1} - \alpha\right),
```

for ``t \in [0, 2]`` with history function

```math
\phi_1(t) = 0.33 - t / 10, 
```
```math
\phi_2(t) = 2.22 + t / 10,
```

for ``t \leq 0``, where ``\alpha = 0.1``, ``\rho = 2.9``, and ``\tau = 0.42``.

# References

Kuang, Y. (1991). On neutral delay logistics Gauss-type predator-prey systems, Dyn. Stab.
Systems (6), pp. 173-189.
"""
const prob_dde_DDETST_E2 = let α = 0.1, ρ = 2.9, τ = 0.42
    global function f_dde_DDETST_E2!(du, u, h, p, t)
        v = u[1]^2
        z = u[2] * v / (v + 1)

        du[1] = u[1] * (1 - h(p, t - τ; idxs = 1) - ρ * h(p, t - 1, Val{1}; idxs = 1)) - z
        du[2] = z - α * u[2]

        nothing
    end

    global function h_dde_DDETST_E2(p, t; idxs::Union{Nothing, Int} = nothing)
        t ≤ 0 || error("history function is only implemented for t ≤ 0")

        if idxs === nothing
            [0.33 - 0.1 * t, 2.22 + 0.1 * t]
        elseif idxs == 1
            0.33 - 0.1 * t
        elseif idxs == 2
            2.22 + 0.1 * t
        else
            error("delay differential equation consists of two components")
        end
    end

    global function h_dde_DDETST_E2(p, t, ::Type{Val{1}};
            idxs::Union{Nothing, Int} = nothing)
        t ≤ 0 || error("history function is only implemented for t ≤ 0")

        if idxs === nothing
            [-0.1, -0.1]
        elseif idxs == 1 || idxs == 2
            -0.1
        else
            error("delay differential equation consists of two components")
        end
    end

    DDEProblem(f_dde_DDETST_E2!, h_dde_DDETST_E2, (0.0, 2.0);
        constant_lags = [τ], neutral = true)
end

# Problem F1
@doc raw"""
    prob_dde_DDETST_F1

Delay differential equation

```math
u'(t) = 2 \cos(2t) u(t / 2)^{2 \cos t} + \log(u'(t / 2)) - \log(2 \cos t) - \sin t
```

for ``t \in [0, 1]`` with history function ``\phi(0) = 1`` and ``\phi'(0) = 2``.

# Solution

The analytical solution for ``t \in [0, 1]`` is

```math
u(t) = \exp(\sin(2t)).
```

# References

Jackiewicz, Z. (1981). One step methods for the numerical solution of Volterra functional
differential equations of neutral type, Applicable Anal. (12), pp. 1-11.
"""
prob_dde_DDETST_F1

function f_dde_DDETST_F1(u, h, p, t)
    s, c = sincos(t)
    c2 = 2 * c

    to2 = t / 2
    h0 = h(p, to2)
    h1 = h(p, to2, Val{1})

    2 * cos(2 * t) * h0^c2 + log(h1) - log(c2) - s
end

function h_dde_DDETST_F1(p, t)
    iszero(t) || error("history function is only implemented for t = 0")

    1.0
end

function h_dde_DDETST_F1(p, t, ::Type{Val{1}})
    iszero(t) || error("history function is only implemented for t = 0")

    2.0
end

function fanalytic_dde_DDETST_F1(u₀, ::typeof(h_dde_DDETST_F1), p, t)
    0 ≤ t ≤ 1 && u₀ == 1 ||
        error("analytical solution is only implemented for t ∈ [0, 1] and u(0) = 1")

    exp(sin(2 * t))
end

const prob_dde_DDETST_F1 = DDEProblem(
    DDEFunction(f_dde_DDETST_F1;
        analytic = fanalytic_dde_DDETST_F1),
    h_dde_DDETST_F1, (0.0, 0.1);
    dependent_lags = ((u, p, t) -> t / 2,),
    neutral = true)

# Problem F2
@doc raw"""
    prob_dde_DDETST_F2

Delay differential equation

```math
u'(t) = u'(2t - 0.5)
```

for ``t \in [0.25, 0.499]`` with history function ``\phi(t) = \exp(-t^2)`` and
``\phi'(t) = -2t \exp(-t^2)`` for ``t \leq 0.25``.

# Solution

The analytical solution for ``t \in [0.25, 0.499]`` is

```math
u(t) = u_i(t) = \exp(-4^i t^2 + B_i t + C_i) / 2^i + K_i
```

if ``t \in [x_i, x_{i + 1}]``, where

```math
x_i = (1 - 2^{-i}) / 2, 
```
```math
B_i = 2 (4^{i-1} + B_{i-1}), 
```
```math
C_i = - 4^{i-2} - B_{i-1} / 2 + C_{i-1}, 
```
```math
K_i = - \exp(-4^i x_i^2 + B_i x_i + C_i) / 2^i + u_{i-1}(x_i),
```

and ``B_0 = C_0 = K_0 = 0``.

# References

Neves, K. W. and Thompson, S. (1992). Solution of systems of functional differential
equations with state dependent delays, Technical Report TR-92-009, Computer Science,
Radford University.
"""
prob_dde_DDETST_F2

f_dde_DDETST_F2(u, h, p, t) = h(p, 2 * t - 0.5, Val{1})

function h_dde_DDETST_F2(p, t)
    t ≤ 0.25 || error("history function is only implemented for t ≤ 0.25")

    exp(-t^2)
end

function h_dde_DDETST_F2(p, t, ::Type{Val{1}})
    t ≤ 0.25 || error("history function is only implemented for t ≤ 0.25")

    -2 * t * exp(-t^2)
end

function fanalytic_dde_DDETST_F2(u₀, ::typeof(h_dde_DDETST_F2), p, t)
    0.25 ≤ t ≤ 0.499 && u₀ == exp(-0.25^2) ||
        error("analytical solution is only implemented for t ∈ [0.25, 0.499] and u(0.25) = exp(-0.0625)")

    i = 0
    x = 0.0
    B = C = K = 0.0

    while true
        i += 1
        y = (1 - 2^(-i)) / 2

        if y ≥ t
            return exp(-4^(i - 1) * t^2 + B * t + C) / 2^(i - 1) + K
        else
            K += exp(-4^(i - 1) * x^2 + B * x + C) / 2^(i - 1)
            C += -4^(i - 2) - B / 2
            B = 2 * (4^(i - 1) + B)
            K -= exp(-4^i * x^2 + B * x + C) / 2^i
            x = y
        end
    end
end

const prob_dde_DDETST_F2 = DDEProblem(
    DDEFunction(f_dde_DDETST_F2;
        analytic = fanalytic_dde_DDETST_F2),
    h_dde_DDETST_F2, (0.25, 0.499);
    dependent_lags = ((u, p, t) -> 1 / 2 - t,),
    neutral = true)

# Problem F3
@doc raw"""
    prob_dde_DDETST_F3

Delay differential equation

```math
u'(t) = \exp(-u(t)) + L_3 \left[\sin(u'(\alpha(t))) - \sin\left(\frac{1}{3 + \alpha(t)}\right)\right]
```

for ``t \in [0, 10]`` with history function ``\phi(0) = \log 3`` and ``\phi'(0) = 1 / 3``,
where ``\alpha(t) = 0.5 t (1 - \cos(2 \pi t))`` and ``L_3 = 0.2``.

# Solution

The analytical solution for ``t \in [0, 10]`` is

```math
u(t) = \log(t + 3).
```
"""
prob_dde_DDETST_F3

let L₃ = 0.2
    global function f_dde_DDETST_F3(u, h, p, t)
        α = 0.5 * t * (1 - cos(2 * π * t))

        exp(-u) + L₃ * (sin(h(p, α, Val{1})) - sin(inv(3 + α)))
    end
end

function h_dde_DDETST_F345(p, t)
    iszero(t) || error("history function is only implemented for t = 0")

    log(3)
end

function h_dde_DDETST_F345(p, t, ::Type{Val{1}})
    iszero(t) || error("history function is only implemented for t = 0")

    inv(3)
end

function fanalytic_dde_DDETST_F345(u₀, ::typeof(h_dde_DDETST_F345), p, t)
    0 ≤ t ≤ 10 && u₀ == log(3) ||
        error("analytical solution is only implemented for t ∈ [0, 10] and u(0) = log 3")

    log(t + 3)
end

const prob_dde_DDETST_F3 = DDEProblem(
    DDEFunction(f_dde_DDETST_F3;
        analytic = fanalytic_dde_DDETST_F345),
    h_dde_DDETST_F345, (0.0, 10.0);
    dependent_lags = ((u, p, t) -> 0.5 * t *
                                   (1 + cos(2 * π * t)),),
    neutral = true)

# Problem F4
"""
    prob_dde_DDETST_F4

Same delay differential equation as [`prob_dde_DDETST_F3`](@ref) with ``L_3 = 0.4``.
"""
prob_dde_DDETST_F4

let L₃ = 0.4
    global function f_dde_DDETST_F4(u, h, p, t)
        α = 0.5 * t * (1 - cos(2 * π * t))

        exp(-u) + L₃ * (sin(h(p, α, Val{1})) - sin(inv(3 + α)))
    end
end

const prob_dde_DDETST_F4 = remake(prob_dde_DDETST_F3;
    f = DDEFunction(f_dde_DDETST_F4;
        analytic = fanalytic_dde_DDETST_F345))

# Problem F5
"""
    prob_dde_DDETST_F5

Same delay differential equation as [`prob_dde_DDETST_F3`](@ref) with ``L_3 = 0.6``.
"""
prob_dde_DDETST_F5

let L₃ = 0.6
    global function f_dde_DDETST_F5(u, h, p, t)
        α = 0.5 * t * (1 - cos(2 * π * t))

        exp(-u) + L₃ * (sin(h(p, α, Val{1})) - sin(inv(3 + α)))
    end
end

const prob_dde_DDETST_F5 = remake(prob_dde_DDETST_F3;
    f = DDEFunction(f_dde_DDETST_F5;
        analytic = fanalytic_dde_DDETST_F345))

# Problem G1
@doc raw"""
    prob_dde_DDETST_G1

Delay differential equation

```math
u'(t) = - u'(t - u(t)^2 / 4)
```

for ``t \in [0, 1]`` with history function ``\phi(t) = 1 - t`` for ``t \leq 0`` and
``\phi'(t) = -1`` for ``t < 0``.

# Solution

The analytical solution for ``t \in [0, 1]`` is

```math
u(t) = t + 1.
```

# References

El'sgol'ts, L. E. and Norkin, S. B. (1973). Introduction to the Theory and Application of
Differential Equations with Deviating Arguments, Academic Press, New York, p. 44.
"""
prob_dde_DDETST_G1

f_dde_DDETST_G1(u, h, p, t) = -h(p, t - u^2 / 4, Val{1})

function h_dde_DDETST_G1(p, t)
    t ≤ zero(t) || error("history function is only implemented for t ≤ 0")

    1 - t
end

function h_dde_DDETST_G1(p, t, ::Type{Val{1}})
    t < zero(t) || error("history function is only implemented for t < 0")

    -1.0
end

function fanalytic_dde_DDETST_G1(u₀, ::typeof(h_dde_DDETST_G1), p, t)
    0 ≤ t ≤ 1 && u₀ == 1 ||
        error("analytical solution is only implemented for t ∈ [0, 1] and u(0) = 1")

    t + 1
end

const prob_dde_DDETST_G1 = DDEProblem(
    DDEFunction(f_dde_DDETST_G1;
        analytic = fanalytic_dde_DDETST_G1),
    h_dde_DDETST_G1, (0.0, 1.0);
    dependent_lags = ((u, p, t) -> u^2 / 4,),
    neutral = true)

# Problem G2
@doc raw"""
    prob_dde_DDETST_G2

Delay differential equation

```math
u'(t) = - u'(u(t) - 2)
```

for ``t \in [0, 1]`` with history function ``\phi(t) = 1 - t`` for ``t \leq 0`` and
``\phi'(t) = -1`` for ``t < 0``.

# Solution

The analytical solution for ``t \in [0, 1]`` is

```math
u(t) = t + 1.
```

El'sgol'ts, L. E. and Norkin, S. B. (1973). Introduction to the Theory and Application of
Differential Equations with Deviating Arguments, Academic Press, New York, pp. 44-45.
"""
prob_dde_DDETST_G2

f_dde_DDETST_G2(u, h, p, t) = -h(p, u - 2, Val{1})

function h_dde_DDETST_G2(p, t)
    t ≤ zero(t) || error("history function is only implemented for t ≤ 0")

    1 - t
end

function h_dde_DDETST_G2(p, t, ::Type{Val{1}})
    t < zero(t) || error("history function is only implemented for t < 0")

    -1.0
end

function fanalytic_dde_DDETST_G2(u₀, ::typeof(h_dde_DDETST_G2), p, t)
    0 ≤ t ≤ 1 && u₀ == 1.0 ||
        error("analytical solution is only implemented for t ∈ [0, 1] and u(0) = 1")

    t + 1
end

const prob_dde_DDETST_G2 = DDEProblem(
    DDEFunction(f_dde_DDETST_G2;
        analytic = fanalytic_dde_DDETST_G2),
    h_dde_DDETST_G2, (0.0, 1.0);
    dependent_lags = ((u, p, t) -> t + 2 - u,),
    neutral = true)

# Problem H1
@doc raw"""
    prob_dde_DDETST_H1

Delay differential equation

```math
u'(t) = - \frac{4 t u(t)^2}{4 + \log(\cos(2t))^2} + \tan(2t) + 0.5 \arctan\left(u'\left(\frac{t u(t)^2}{1 + u(t)^2}\right)\right)
```

for ``t \in [0, 0.225 \pi]`` with history function ``\phi(0) = 0`` and ``\phi'(0) = 0``.

# Solution

The analytical solution for ``t \in [0, 0.225 \pi]`` is

```math
u(t) = - \log(\cos(2t)) / 2.
```

# References

Castleton, R. N. and Grimm, L. J. (1973). A first order method for differential equations of
neutral type, Math. Comput. (27), pp. 571-577.
"""
prob_dde_DDETST_H1

function f_dde_DDETST_H1(u, h, p, t)
    v = u^2
    z = t * v

    -4 * z / (4 + log(cos(2 * t))^2) + tan(2 * t) + 0.5 * atan(h(p, z / (1 + v), Val{1}))
end

function h_dde_DDETST_H1(p, t)
    iszero(t) || error("history function is only implemented for t = 0")

    0.0
end

function h_dde_DDETST_H1(p, t, ::Type{Val{1}})
    iszero(t) || error("history function is only implemented for t = 0")

    0.0
end

function fanalytic_dde_DDETST_H1(u₀, ::typeof(h_dde_DDETST_H1), p, t)
    0 ≤ t ≤ 0.225 * π && iszero(u₀) ||
        error("analytical solution is only implemented for t ∈ [0, 0.225π] and u(0) = 0")

    -log(cos(2 * t)) / 2
end

const prob_dde_DDETST_H1 = DDEProblem(
    DDEFunction(f_dde_DDETST_H1;
        analytic = fanalytic_dde_DDETST_H1),
    h_dde_DDETST_H1, (0.0, 0.225 * π);
    dependent_lags = ((u, p, t) -> t / (1 + u^2),),
    neutral = true)

# Problem H2
@doc raw"""
    prob_dde_DDETST_H2

Delay differential equation

```math
u'(t) = \cos(t) (1 + u(t u(t)^2)) + L_3 u(t) u'(t u(t)^2) + (1 - L_3) \sin(t) \cos(t \sin(t)^2) - \sin(t + t \sin(t)^2)
```

for ``t \in [0, \pi]`` with history function ``\phi(0) = 0`` and ``\phi'(0) = 1``, where
``L_3 = 0.1``.

# Solution

The analytical solution for ``t \in [0, \pi]`` is

```math
u(t) = \sin(t).
```

# References

Hayashi, H. (1996). Numerical solution of retarded and neutral delay differential equations
using continuous Runge-Kutta methods, PhD thesis, Department of Computer Science, University
of Toronto, Toronto, Canada.
"""
prob_dde_DDETST_H2

let L₃ = 0.1
    global function f_dde_DDETST_H2(u, h, p, t)
        s, c = sincos(t)
        a = t * s^2

        z = t * u^2

        c * (1 + h(p, z)) + L₃ * u * h(p, z, Val{1}) + (1 - L₃) * s * cos(a) - sin(t + a)
    end
end

function h_dde_DDETST_H234(p, t)
    iszero(t) || error("history function is only implemented for t = 0")

    0.0
end

function h_dde_DDETST_H234(p, t, ::Type{Val{1}})
    iszero(t) || error("history function is only implemented for t = 0")

    1.0
end

function fanalytic_dde_DDETST_H234(u₀, ::typeof(h_dde_DDETST_H234), p, t)
    0 ≤ t ≤ π && iszero(u₀) ||
        error("analytical solution is only implemented for t ∈ [0, π] and u(0) = 0")

    sin(t)
end

const prob_dde_DDETST_H2 = DDEProblem(
    DDEFunction(f_dde_DDETST_H2;
        analytic = fanalytic_dde_DDETST_H234),
    h_dde_DDETST_H234, (0.0, π);
    dependent_lags = ((u, p, t) -> t * (1 - u^2),),
    neutral = true)

# Problem H3
"""
    prob_dde_DDETST_H3

Same delay differential equation as [`prob_dde_DDETST_H2`](@ref) with ``L_3 = 0.3``.

# References

Hayashi, H. (1996). Numerical solution of retarded and neutral delay differential equations
using continuous Runge-Kutta methods, PhD thesis, Department of Computer Science, University
of Toronto, Toronto, Canada.
"""
prob_dde_DDETST_H3

let L₃ = 0.3
    global function f_dde_DDETST_H3(u, h, p, t)
        s, c = sincos(t)
        a = t * s^2

        z = t * u^2

        c * (1 + h(p, z)) + L₃ * u * h(p, z, Val{1}) + (1 - L₃) * s * cos(a) - sin(t + a)
    end
end

const prob_dde_DDETST_H3 = remake(prob_dde_DDETST_H2;
    f = DDEFunction(f_dde_DDETST_H3;
        analytic = fanalytic_dde_DDETST_H234))

# Problem H4
"""
    prob_dde_DDETST_H4

Same delay differential equation as [`prob_dde_DDETST_H2`](@ref) with ``L_3 = 0.5``.

# References

Hayashi, H. (1996). Numerical solution of retarded and neutral delay differential equations
using continuous Runge-Kutta methods, PhD thesis, Department of Computer Science, University
of Toronto, Toronto, Canada.
"""
prob_dde_DDETST_H4

let L₃ = 0.5
    global function f_dde_DDETST_H4(u, h, p, t)
        s, c = sincos(t)
        a = t * s^2

        z = t * u^2

        c * (1 + h(p, z)) + L₃ * u * h(p, z, Val{1}) + (1 - L₃) * s * cos(a) - sin(t + a)
    end
end

const prob_dde_DDETST_H4 = remake(prob_dde_DDETST_H2;
    f = DDEFunction(f_dde_DDETST_H4;
        analytic = fanalytic_dde_DDETST_H234))

#=
The following examples are taken from the source code of the delay differential equation
solver RADAR5 by Nicola Guglielmi and Ernst Hairer:
http://www.unige.ch/~hairer/radar5-v2.1.tar
=#

# OREGONATOR example
@doc raw"""
    prob_dde_RADAR5_oregonator

Delay differential equation model from chemical kinetics, given by

```math
\begin{align*}
  u_1'(t) &= - k_1 A u_2(t) - k_2 u_1(t) u_2(t - \tau) + k_3 B u_1(t) - 2 k_4 u_1(t)^2, \\
  u_2'(t) &= - k_1 A u_2(t) - k_2 u_1(t) u_2(t - \tau) + f k_3 B u_1(t),
\end{align*}
```

for ``t \in [0, 100.5]`` with history function

```math
\begin{align*}
  \phi_1(t) &= 1e-10, \\
  \phi_2(t) &= 1e-5,
\end{align*}
```

for ``t \leq 0``, where ``k_1 = 1.34``, ``k_2 = 1.6e9``, ``k_3 = 8000``, ``k_4 = 4e7``,
``k_5 = 1``, ``f = 1``, ``A = 0.06``, ``B = 0.06``, and ``\tau = 0.15``.

# References

Epstein, I. and Luo, Y. (1991). Differential delay equations in chemical kinetics. Nonlinear
models, Journal of Chemical Physics (95), pp. 244-254.
"""
const prob_dde_RADAR5_oregonator =
  let k₁ = 1.34, k₂ = 1.6e9, k₃ = 8_000, k₄ = 4e7, f = 1, A = 0.06, B = 0.06, τ = 0.15
    global function f_dde_RADAR5_oregonator!(du, u, h, p, t)
      # past value of the second component
      v = h(p, t - τ; idxs = 2)

      # precalculations
      a = - k₁ * A * u[2] - k₂ * u[1] * v
      b = k₃ * B * u[1]

      du[1] = a + b - 2 * k₄ * u[1]^2
      du[2] = a + f * b

      nothing
    end

    global function h_dde_RADAR5_oregonator(p, t; idxs::Union{Nothing,Int} = nothing)
      t ≤ 0 || error("history function is only implemented for t ≤ 0")

      if idxs === nothing
        [1e-10, 1e-5]
      elseif idxs == 2
        1e-5
      else
        error("history function is only implemented for the second component")
      end
    end

    DDEProblem(f_dde_RADAR5_oregonator!, h_dde_RADAR5_oregonator, (0.0, 100.5);
               constant_lags = [τ])
  end

# ROBERTSON example
@doc raw"""
    prob_dde_RADAR5_robertson

Delay differential equation model of a chemical reaction with steady state solution, given
by

```math
\begin{align*}
  u_1'(t) &= - a u_1(t) + b u_2(t - \tau) u_3(t), \\
  u_2'(t) &= a u_1(t) - b u_2(t - \tau) u_3(t) - c u_2(t)^2, \\
  u_3'(t) &= c u_2(t)^2,
\end{align*}
```

for ``t \in [0, 10e10]`` with history function ``\phi_1(0) = 1``, ``\phi_2(t) = 0`` for
``t \in [-\tau, 0]``, and ``\phi_3(0) = 0``, where ``a = 0.04``, ``b = 10_000``,
``c = 3e7``, and ``\tau = 0.01``.

# References

Guglielmi, N. and Hairer, E. (2001). Implementing Radau IIA methods for stiff delay
differential equations, Computing (67), pp. 1-12.
"""
const prob_dde_RADAR5_robertson =
  let a = 0.04, b = 10_000, c = 3e7, τ = 0.01
    global function f_dde_RADAR5_robertson!(du, u, h, p, t)
      # past value of the second component
      v = h(p, t - τ; idxs = 2)

      # precalculations
      x = a * u[1] - b * v * u[3]
      y = c * u[2]^2

      du[1] = - x
      du[2] = x - y
      du[3] = y

      nothing
    end

    global function h_dde_RADAR5_robertson(p, t; idxs::Union{Nothing,Int} = nothing)
      - τ ≤ t ≤ 0 || error("history function is only implemented for t ∈ [-τ, 0]")

      if idxs === nothing
        [1.0, 0.0, 0.0]
      elseif idxs == 2
        0.0
      else
        error("history function is only implemented for the second component")
      end
    end

    DDEProblem(f_dde_RADAR5_robertson!, h_dde_RADAR5_robertson, (0.0, 1e10);
               constant_lags = [τ])
  end

# WALTMAN example
@doc raw"""
    prob_dde_RADAR5_waltman

Delay differential equation model of antibody production, given by

```math
\begin{align*}
  u_1'(t) &= - r u_1(t) u_2(t) - s u_1(t) u_4(t), \\
  u_2'(t) &= - r u_1(t) u_2(t) + \alpha r u_1(u_5(t)) u_2(u_5(t)) [t \geq t_0], \\
  u_3'(t) &= r u_1(t) u_2(t), \\
  u_4'(t) &= - s u_1(t) u_4(t) - \gamma u_4(t) + \beta r u_1(u_6(t)) u_2(u_6(t)) [t > t_1], \\
  u_5'(t) &= [t \geq t_0] \frac{u_1(t) u_2(t) + u_3(t)}{u_1(u_5(t)) u_2(u_5(t)) + u_3(u_5(t))}, \\
  u_6'(t) &= [t \geq t_1] \frac{1e-12 + u_2(t) + u_3(t)}{1e-12 + u_2(u_6(t)) + u_3(u_6(t))},
\end{align*}
```

for ``t \in [0, 300]`` with history function

```math
\begin{align*}
  \phi_1(t) &= \phi_0, \\
  \phi_2(t) &= 1e-15, \\
  \phi_3(t) &= 0, \\
  \phi_4(t) &= 0, \\
  \phi_5(t) &= 0, \\
  \phi_6(t) &= 0,
\end{align*}
```

for ``t \leq 0``, where ``\alpha = 1.8``, ``\beta = 20``, ``\gamma = 0.002``, ``r = 5e4``,
``s = 1e5``, ``t_0 = 32``, ``t_1 = 119``, and ``\phi_0 = 0.75e-4``.

# References

Waltman, P. (1978). A threshold model of antigen-stimulated antibody production, Theoretical
Immunology (8), pp. 437-453.
"""
prob_dde_RADAR5_waltman

let α = 1.8, β = 20, γ = 0.002, r = 50_000, s = 100_000
  global function f_dde_RADAR5_waltman!(du, u, h, p, t)
    # time of switches
    t₀, t₁ = p.t₀, p.t₁

    # precalculations
    u₁u₂ = u[1] * u[2]
    a = r * u₁u₂
    b = s * u[1] * u[4]

    du[1] = - a - b
    du[3] = a

    if t < t₀
      du[2] = - a
      du[5] = 0
    else
      v = h(p, u[5])
      v₁v₂ = v[1] * v[2]
      du[2] = - a + α * r * v₁v₂
      du[5] = (u₁u₂ + u[3]) / (v₁v₂ + v[3])
    end

    if t < t₁
      du[4] = - b - γ * u[4]
      du[6] = 0
    else
      w = h(p, u[6])
      du[4] = - b - γ * u[4] + β * r * w[1] * w[2]
      du[6] = (1e-12 + u[2] + u[3]) / (1e-12 + w[2] + w[3])
    end

    nothing
  end
end

function h_dde_RADAR5_waltman(p, t)
  t ≤ 0 || error("history function is only implemented for t ≤ 0")

  [p.ϕ₀, 1e-15, 0.0, 0.0, 0.0, 0.0]
end

const prob_dde_RADAR5_waltman =
  DDEProblem(f_dde_RADAR5_waltman!, (p, t) -> h_dde_RADAR5_waltman(p, 0.0),
             h_dde_RADAR5_waltman, (0.0, 300.0), (ϕ₀ = 0.75e-4, t₀ = 32, t₁ = 119);
             dependent_lags = ((u, p, t) -> t - u[5], (u, p, t) -> t - u[6]),
             order_discontinuity_t0 = 1)
const prob_dde_RADAR5_waltman_1 = prob_dde_RADAR5_waltman

@doc raw"""
    prob_dde_RADAR5_waltman_2

Same delay differential equation as [`prob_dde_RADAR5_waltman`] with ``t_0 = 32``,
``t_1 = 111``, and ``\phi_0 = 0.5e-4``.

# References

Waltman, P. (1978). A threshold model of antigen-stimulated antibody production, Theoretical
Immunology (8), pp. 437-453.
"""
const prob_dde_RADAR5_waltman_2 =
  remake(prob_dde_RADAR5_waltman; p = (ϕ₀ = 0.5e-4, t₀ = 32, t₁ = 111))

@doc raw"""
    prob_dde_RADAR5_waltman_3

Same delay differential equation as [`prob_dde_RADAR5_waltman`] with ``t_0 = 33``,
``t_1 = 145``, and ``\phi_0 = 1e-5``.

# References

Waltman, P. (1978). A threshold model of antigen-stimulated antibody production, Theoretical
Immunology (8), pp. 437-453.
"""
const prob_dde_RADAR5_waltman_3 =
  remake(prob_dde_RADAR5_waltman; p = (ϕ₀ = 1e-5, t₀ = 33, t₁ = 145))

@doc raw"""
    prob_dde_RADAR5_waltman_4

Same delay differential equation as [`prob_dde_RADAR5_waltman`] with ``t_0 = 34``,
``t_1 = 163``, and ``\phi_0 = 0.75e-5``.

# References

Waltman, P. (1978). A threshold model of antigen-stimulated antibody production, Theoretical
Immunology (8), pp. 437-453.
"""
const prob_dde_RADAR5_waltman_4 =
  remake(prob_dde_RADAR5_waltman; p = (ϕ₀ = 0.75e-5, t₀ = 34, t₁ = 163))

@doc raw"""
    prob_dde_RADAR5_waltman_5

Same delay differential equation as [`prob_dde_RADAR5_waltman`] with ``t_0 = 35``,
``t_1 = 197``, and ``\phi_0 = 0.5e-5``.

# References

Waltman, P. (1978). A threshold model of antigen-stimulated antibody production, Theoretical
Immunology (8), pp. 437-453.
"""
const prob_dde_RADAR5_waltman_5 =
  remake(prob_dde_RADAR5_waltman; p = (ϕ₀ = 0.5e-5, t₀ = 35, t₁ = 197))

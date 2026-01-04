module DAEProblemLibrary

using DiffEqBase, Markdown

#DAE Example Problems
export prob_dae_resrob, prob_dae_transamp

### DAE Problems

f = function (r, yp, y, p, tres)
    r[1] = -0.04 * y[1] + 1.0e4 * y[2] * y[3]
    r[2] = -r[1] - 3.0e7 * y[2] * y[2] - yp[2]
    r[1] -= yp[1]
    return r[3] = y[1] + y[2] + y[3] - 1.0
end
u0 = [1.0, 0, 0]
du0 = [-0.04, 0.04, 0.0]

@doc doc"""
The Robertson biochemical reactions in DAE form

```math
\begin{align*}
\frac{dy₁}{dt} &= -k₁y₁+k₃y₂y₃ \\
\frac{dy₂}{dt} &=  k₁y₁-k₂y₂^2-k₃y₂y₃ \\
1 &= y₁ + y₂ + y₃
\end{align*}
```
where ``k₁=0.04``, ``k₂=3×10^7``, ``k₃=10^4``. For details, see:
Hairer Norsett Wanner Solving Ordinary Differential Equations I - Nonstiff Problems Page 129
Usually solved on ``[0,10^{11}]``
"""
prob_dae_resrob = DAEProblem(f, du0, u0, (0.0, 100000.0))

C = [k * 1.0e-6 for k in 1:5]
Ub = 6
UF = 0.026
α = 0.99
β = 1.0e-6
R0 = 1000
R = 9000
Ue(t) = 0.1 * sin(200 * π * t)
g(x) = 1.0e-6 * (exp(x / 0.026) - 1)

function transamp(out, du, u, p, t)
    y1, y2, y3, y4, y5, y6, y7, y8 = u
    out[1] = -Ue(t) / R0 + y1 / R0 + C[1] * du[1] - C[1] * du[2]
    out[2] = -Ub / R + y2 * 2 / R - (α - 1) * g(y2 - y3) - C[1] * du[1] + C[1] * du[2]
    out[3] = -g(y2 - y3) + y3 / R + C[2] * du[3]
    out[4] = -Ub / R + y4 / R + α * g(y2 - y3) + C[3] * du[4] - C[3] * du[5]
    out[5] = -Ub / R + y5 * 2 / R - (α - 1) * g(y5 - y6) - C[3] * du[4] + C[3] * du[5]
    out[6] = -g(y5 - y6) + y6 / R + C[4] * du[6]
    out[7] = -Ub / R + y7 / R + α * g(y5 - y6) + C[5] * du[7] - C[5] * du[8]
    return out[8] = y8 / R - C[5] * du[7] + C[5] * du[8]
end

u0 = [0, Ub / 2, Ub / 2, Ub, Ub / 2, Ub / 2, Ub, 0]
du0 = [
    51.338775,
    51.338775,
    -Ub / (2 * (C[2] * R)),
    -24.9757667,
    -24.9757667,
    -Ub / (2 * (C[4] * R)),
    -10.00564453,
    -10.00564453,
]

@doc doc"""
The Transistor Amplifier model

```math
M\frac{dy}{dt} = f(t,y), \quad y(0)=y_0,\quad y'(0)=y_0'
```

```math
M = \begin{pmatrix}
-C_1 &  C_1 &    0 &    0 &    0 &    0 &    0 &    0 \\
 C_1 & -C_1 &    0 &    0 &    0 &    0 &    0 &    0 \\
   0 &    0 & -C_2 &    0 &    0 &    0 &    0 &    0 \\
   0 &    0 &    0 & -C_3 &  C_3 &    0 &    0 &    0 \\
   0 &    0 &    0 &  C_3 & -C_3 &    0 &    0 &    0 \\
   0 &    0 &    0 &    0 &    0 & -C_4 &    0 &    0 \\
   0 &    0 &    0 &    0 &    0 &    0 & -C_5 &  C_5 \\
   0 &    0 &    0 &    0 &    0 &    0 &  C_5 & -C_5
\end{pmatrix}
```

```math
f(t, y)=\begin{pmatrix}
-\frac{U_e(t)}{R_0} + \frac{y_1}{R_0} \\
-\frac{U_b}{R_2} + y_2\left(\frac{1}{R_1}+\frac{1}{R_2}\right) - (α-1) g\left(y_2-y_3\right) \\
-g\left(y_2-y_3\right) + \frac{y_3}{R_3} \\
-\frac{U_b}{R_4} + \frac{y_4}{R_4} + α g\left(y_2-y_3\right) \\
-\frac{U_b}{R_6} + y_5\left(\frac{1}{R_5}+\frac{1}{R_6}\right) - (α-1) g\left(y_5-y_6\right) \\
-g\left(y_5-y_6\right) + \frac{y_6}{R_7} \\
-\frac{U_b}{R_8} + \frac{y_7}{R_8} + α g\left(y_5-y_6\right) \\
\frac{y_8}{R_9}
\end{pmatrix}
```

## Reference
DAE testset: <https://archimede.uniba.it/~testset/problems/transamp.php>
"""
prob_dae_transamp = DAEProblem(transamp, du0, u0, (0.0, 0.2))

end # module

module DAEProblemLibrary

using DiffEqBase, Markdown

#DAE Example Problems
export prob_dae_resrob, prob_dae_transamp

### DAE Problems

f = function (r, yp, y, p, tres)
    r[1] = -0.04 * y[1] + 1.0e4 * y[2] * y[3]
    r[2] = -r[1] - 3.0e7 * y[2] * y[2] - yp[2]
    r[1] -= yp[1]
    r[3] = y[1] + y[2] + y[3] - 1.0
end
u0 = [1.0, 0, 0]
du0 = [-0.04, 0.04, 0.0]

@doc doc"""
The Robertson biochemical reactions in DAE form

```math
\frac{dy₁}{dt} = -k₁y₁+k₃y₂y₃
```
```math
\frac{dy₂}{dt} =  k₁y₁-k₂y₂^2-k₃y₂y₃
```
```math
1 = y₁ + y₂ + y₃
```
where ``k₁=0.04``, ``k₂=3\times10^7``, ``k₃=10^4``. For details, see:
Hairer Norsett Wanner Solving Ordinary Differential Equations I - Nonstiff Problems Page 129
Usually solved on ``[0,1e11]``
"""
prob_dae_resrob = DAEProblem(f, du0, u0, (0.0, 100000.0))

C = [k * 1e-6 for k in 1:5]
Ub = 6
UF = 0.026
α = 0.99
β = 1e-6
R0 = 1000
R = 9000
Ue(t) = 0.1 * sin(200 * π * t)
g(x) = 1e-6 * (exp(x / 0.026) - 1)

function transamp(out, du, u, p, t)
    y1, y2, y3, y4, y5, y6, y7, y8 = u
    out[1] = -Ue(t) / R0 + y1 / R0 + C[1] * du[1] - C[1] * du[2]
    out[2] = -Ub / R + y2 * 2 / R - (α - 1) * g(y2 - y3) - C[1] * du[1] + C[1] * du[2]
    out[3] = -g(y2 - y3) + y3 / R + C[2] * du[3]
    out[4] = -Ub / R + y4 / R + α * g(y2 - y3) + C[3] * du[4] - C[3] * du[5]
    out[5] = -Ub / R + y5 * 2 / R - (α - 1) * g(y5 - y6) - C[3] * du[4] + C[3] * du[5]
    out[6] = -g(y5 - y6) + y6 / R + C[4] * du[6]
    out[7] = -Ub / R + y7 / R + α * g(y5 - y6) + C[5] * du[7] - C[5] * du[8]
    out[8] = y8 / R - C[5] * du[7] + C[5] * du[8]
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
    -10.00564453
]

@doc doc"""
The Transistor Amplifier model

M\frac{dy}{dt}=f(t,y),\quad y(0)=y_0,\quad y'(0)=y_0'

```math
M=\left(\begin{array}{cccccccc}
-C_{1} & C_{1} & 0 & 0 & 0 & 0 & 0 & 0 \\
C_{1} & -C_{1} & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & -C_{2} & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & -C_{3} & C_{3} & 0 & 0 & 0 \\
0 & 0 & 0 & C_{3} & -C_{3} & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & -C_{4} & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & -C_{5} & C_{5} \\
0 & 0 & 0 & 0 & 0 & 0 & C_{5} & -C_{5}
\end{array}\right)
```

```math
f(t, y)=\left(\begin{array}{c}
-\frac{U_{e}(t)}{R_{0}}+\frac{y_{1}}{R_{0}} \\
-\frac{U_{b}}{R_{2}}+y_{2}\left(\frac{1}{R_{1}}+\frac{1}{R_{2}}\right)-(\alpha-1) g\left(y_{2}-y_{3}\right) \\
-g\left(y_{2}-y_{3}\right)+\frac{y_{3}}{R_{3}} \\
-\frac{U_{b}}{R_{4}}+\frac{y_{4}}{R_{4}}+\alpha g\left(y_{2}-y_{3}\right) \\
-\frac{U_{b}}{R_{6}}+y_{5}\left(\frac{1}{R_{5}}+\frac{1}{R_{6}}\right)-(\alpha-1) g\left(y_{5}-y_{6}\right) \\
-g\left(y_{5}-y_{6}\right)+\frac{y_{6}}{R_{7}} \\
-\frac{U_{b}}{R_{8}}+\frac{y_{7}}{R_{8}}+\alpha g\left(y_{5}-y_{6}\right) \\
\frac{y_{8}}{R_{9}}
\end{array}\right)
```

## Reference
DAE testset: https://archimede.uniba.it/~testset/problems/transamp.php
"""
prob_dae_transamp = DAEProblem(transamp, du0, u0, (0.0, 0.2))

end # module

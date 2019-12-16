using DiffEqBase, Markdown

#DAE Example Problems
export prob_dae_resrob

### DAE Problems

f = function (r, yp, y, p, tres)
    r[1]  = -0.04*y[1] + 1.0e4*y[2]*y[3]
    r[2]  = -r[1] - 3.0e7*y[2]*y[2] - yp[2]
    r[1] -=  yp[1]
    r[3]  =  y[1] + y[2] + y[3] - 1.0
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
prob_dae_resrob = DAEProblem(f,du0,u0,(0.0,100000.0))

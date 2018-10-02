## Lotka-Volterra

lotka = @ode_def_all LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d

"""
Lotka-Voltera Equations (Non-stiff)

```math
\\frac{dx}{dt} = ax - bxy
\\frac{dy}{dt} = -cy + dxy
```

with initial condition ``x=y=1``
"""
prob_ode_lotkavoltera = ODEProblem(lotka,[1.0,1.0],(0.0,1.0),[1.5,1.0,3.0,1.0])

## Fitzhugh-Nagumo

fitz = @ode_def_all FitzhughNagumo begin
  dv = v - v^3/3 -w + l
  dw = τinv*(v +  a - b*w)
end a b τinv l
"""
Fitzhugh-Nagumo (Non-stiff)

```math
\\frac{dv}{dt} = v - \\frac{v^3}{3} - w + I_{est}
τ \\frac{dw}{dt} = v + a -bw
```

with initial condition ``v=w=1``
"""
prob_ode_fitzhughnagumo = ODEProblem(fitz,[1.0;1.0],(0.0,1.0),(0.7,0.8,1/12.5,0.5))

#Van der Pol Equations
van = @ode_def_all VanDerPol begin
  dy = μ*((1-x^2)*y - x)
  dx = 1*y
end μ

"""
Van der Pol Equations

```math
\\begin{align}
\\frac{dx}{dt} &= y \\\\
\\frac{dy}{dt} &= μ((1-x^2)y -x)
\\end{align}
```

with ``μ=1.0`` and ``u0=[0,\\sqrt{3}]``

Non-stiff parameters.
"""
prob_ode_vanderpol = ODEProblem(van,[0;sqrt(3)],(0.0,1.0),1.0)

"""Van der Pol Equations

```math
\\begin{align}
\\frac{dx}{dt} &= y \\\\
\\frac{dy}{dt} &= μ(1-x^2)y -x
\\end{align}
```

with ``μ=10^6`` and ``u0=[0,\\sqrt{3}]``

Stiff parameters.
"""
prob_ode_vanstiff = ODEProblem(van,[0;sqrt(3)],(0.0,1.0),1e6)

# ROBER

rober = @ode_def_all Rober begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁ k₂ k₃

"""
The Robertson biochemical reactions: (Stiff)

```math
\\begin{align}
\\frac{dy₁}{dt} &= -k₁y₁+k₃y₂y₃  \\\\
\\frac{dy₂}{dt} &=  k₁y₁-k₂y₂^2-k₃y₂y₃ \\\\
\\frac{dy₃}{dt} &=  k₂y₂^2
\\end{align}
```

where ``k₁=0.04``, ``k₂=3\\times10^7``, ``k₃=10^4``. For details, see:

Hairer Norsett Wanner Solving Ordinary Differential Equations I - Nonstiff Problems Page 129

Usually solved on `[0,1e11]`
"""
prob_ode_rober = ODEProblem(rober,[1.0;0.0;0.0],(0.0,1e11),(0.04,3e7,1e4))

# Three Body
const threebody_μ = big(0.012277471); const threebody_μ′ = 1 - threebody_μ

threebody = (du,u,p,t) -> begin
  # 1 = y₁
  # 2 = y₂
  # 3 = y₁'
  # 4 = y₂'
  D₁ = ((u[1]+threebody_μ)^2 + u[2]^2)^(3/2)
  D₂ = ((u[1]-threebody_μ′)^2 + u[2]^2)^(3/2)
  du[1] = u[3]
  du[2] = u[4]
  du[3] = u[1] + 2u[4] - threebody_μ′*(u[1]+threebody_μ)/D₁ - threebody_μ*(u[1]-threebody_μ′)/D₂
  du[4] = u[2] - 2u[3] - threebody_μ′*u[2]/D₁ - threebody_μ*u[2]/D₂
end
"""
The ThreeBody problem as written by Hairer: (Non-stiff)

```math
\\begin{align}
y₁′′ &= y₁ + 2y₂′ - μ′\\frac{y₁+μ}{D₁} - μ\\frac{y₁-μ′}{D₂} \\\\
y₂′′ &= y₂ - 2y₁′ - μ′\\frac{y₂}{D₁} - μ\\frac{y₂}{D₂} \\\\
D₁ &= ((y₁+μ)^2 + y₂^2)^{3/2} \\\\
D₂ &= ((y₁-μ′)^2+y₂^2)^{3/2} \\\\
μ &= 0.012277471 \\\\
μ′ &=1-μ
\\end{align}
```

From Hairer Norsett Wanner Solving Ordinary Differential Equations I - Nonstiff Problems Page 129

Usually solved on `t₀ = 0.0`; `T = "17.0652165601579625588917206249"`
Periodic with that setup.
"""
prob_ode_threebody = ODEProblem(threebody,[0.994, 0.0, 0.0, big(-2.00158510637908252240537862224)],(big(0.0),big(17.0652165601579625588917206249)))

# Rigid Body Equations

rigid = @ode_def_all RigidBody begin
  dy₁  = I₁*y₂*y₃
  dy₂  = I₂*y₁*y₃
  dy₃  = I₃*y₁*y₂
end I₁ I₂ I₃

"""
Rigid Body Equations (Non-stiff)

```math
\\begin{align}
\\frac{dy₁}{dt}  &= I₁y₂y₃ \\\\
\\frac{dy₂}{dt}  &= I₂y₁y₃ \\\\
\\frac{dy₃}{dt}  &= I₃y₁y₂
\\end{align}
```

with ``I₁=-2``, ``I₂=1.25``, and ``I₃=-1/2``.

The initial condition is ``y=[1.0;0.0;0.9]``.

From Solving Differential Equations in R by Karline Soetaert

or Hairer Norsett Wanner Solving Ordinary Differential Equations I - Nonstiff Problems Page 244

Usually solved from 0 to 20.
"""
prob_ode_rigidbody = ODEProblem(rigid,[1.0,0.0,0.9],(0.0,20.0),(-2.0,1.25,-0.5))

# Pleiades Problem

pleiades = (du,u,p,t) -> begin
  x = view(u,1:7)   # x
  y = view(u,8:14)  # y
  v = view(u,15:21) # x′
  w = view(u,22:28) # y′
  du[1:7] .= v
  du[8:14].= w
  for i in 14:28
    du[i] = zero(eltype(u))
  end
  for i=1:7,j=1:7
    if i != j
      r = ((x[i]-x[j])^2 + (y[i] - y[j])^2)^(3/2)
      du[14+i] += j*(x[j] - x[i])/r
      du[21+i] += j*(y[j] - y[i])/r
    end
  end
end
"""
Pleiades Problem (Non-stiff)

```math
\\begin{align}
xᵢ′′ &= \\sum_{j≠i} mⱼ(xⱼ-xᵢ)/rᵢⱼ \\\\
yᵢ′′ &= \\sum_{j≠i} mⱼ(yⱼ-yᵢ)/rᵢⱼ
\\end{align}
```

where

```math
rᵢⱼ = ((xᵢ-xⱼ)^2 + (yᵢ-yⱼ)^2)^{3/2}
```

and initial conditions are

```math
\\begin{align}
x₁(0)&=3  \\\\
x₂(0)&=3  \\\\
x₃(0)&=-1  \\\\
x₄(0)&=-3  \\\\
x₅(0)&=2  \\\\
x₆(0)&=-2  \\\\
x₇(0)&=2  \\\\
y₁(0)&=3  \\\\
y₂(0)&=-3  \\\\
y₃(0)&=2  \\\\
y₄(0)&=0  \\\\
y₅(0)&=0  \\\\
y₆(0)&=-4  \\\\
y₇(0)&=4
\\end{align}
```

and with ``xᵢ′(0)=yᵢ′(0)=0`` except for

```math
\\begin{align}
x₆′(0)&=1.75 \\\\
x₇′(0)&=-1.5 \\\\
y₄′(0)&=-1.25 \\\\
y₅′(0)&=1
\\end{align}
```

From Hairer Norsett Wanner Solving Ordinary Differential Equations I - Nonstiff Problems Page 244

Usually solved from 0 to 3.
"""
prob_ode_pleiades = ODEProblem(pleiades,[3.0,3.0,-1.0,-3.0,2.0,-2.0,2.0,3.0,-3.0,2.0,0,0,-4.0,4.0,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0],(0.0,3.0))



Random.seed!(100)
const mm_A = rand(4,4)
mm_linear = function (du,u,p,t)
  mul!(du,mm_A,u)
end
const MM_linear =Matrix(Diagonal(0.5ones(4)))
mm_f = ODEFunction(mm_linear;analytic = (u0,p,t) -> exp(inv(MM_linear)*mm_A*t)*u0, mass_matrix=MM_linear)
prob_ode_mm_linear = ODEProblem(mm_f,rand(4),(0.0,1.0))

"""
[Hires Problem](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/Hires.ipynb) (Stiff)

It is in the form of ``\\frac{dy}{dt}=f(y), \\quad y(0)=y0,`` with

```math
y \\in ℝ^8, \\quad 0 ≤ t ≤ 321.8122

where ``f`` is defined by

```math
f(y) = \\begin{pmatrix}
−1.71y_1 & +0.43y_2 & +8.32y_3 & +0.0007y_4 & \\\\
1.71y_1 & −8.75y_2 & & & \\\\
−10.03y_3 & +0.43y_4 & +0.035y_5 & & \\\\
8.32y_2 & +1.71y_3 & −1.12y_4 & & \\\\
−1.745y_5 & +0.43y_6 & +0.43y_7 & & \\\\
−280y_6y_8 & +0.69y_4 & +1.71y_5 & −0.43y_6 & +0.69y_7 \\\\
280y_6y_8 & −1.81y_7 & & & \\\\
−280y_6y_8 & +1.81y_7 & & &
\\end{pmatrix}
```

http://www.radford.edu/~thompson/vodef90web/problems/demosnodislin/Demos_Pitagora/DemoHires/demohires.pdf
"""
hires = @ode_def_all Hires begin
  dy1 = -p1*y1 + p2*y2 + p3*y3 + p4
  dy2 = p1*y1 - p5*y2
  dy3 = -p6*y3 + p2*y4 + p7*y5
  dy4 = p3*y2 + p1*y3 - p8*y4
  dy5 = -p9*y5 + p2*y6 + p2*y7
  dy6 = -p10*y6*y8 + p11*y4 + p1*y5 -
           p2*y6 + p11*y7
  dy7 = p10*y6*y8 - p12*y7
  dy8 = -p10*y6*y8 + p12*y7
end p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12
u0 = zeros(8)
u0[1] = 1
u0[8] = 0.0057
prob_ode_hires = ODEProblem(hires,u0,(0.0,321.8122), (1.71, 0.43, 8.32, 0.0007, 8.75,
                                                      10.03, 0.035, 1.12, 1.745, 280.0,
                                                      0.69, 1.81))

"""
[Orego Problem](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/Orego.ipynb) (Stiff)

It is in the form of ``\\frac{dy}{dt}=f(y), \\quad y(0)=y0,`` with

```math
y \\in ℝ^3, \\quad 0 ≤ t ≤ 360

where ``f`` is defined by

```math
f(y) = \\begin{pmatrix}
s(y_2 - y_1(1-qy_1-y_2)) \\\\
(y_3 - y_2(1+y_1))/s \\\\
w(y_1-y_3)
\\end{pmatrix}
```

where ``s=77.27``, ``w=0.161`` and ``q=8.375⋅10^{-6}``.

http://www.radford.edu/~thompson/vodef90web/problems/demosnodislin/Demos_Pitagora/DemoOrego/demoorego.pdf
"""
orego = @ode_def_all Orego begin
  dy1 = p1*(y2+y1*(1-p2*y1-y2))
  dy2 = (y3-(1+y1)*y2)/p1
  dy3 = p3*(y1-y3)
end p1 p2 p3
prob_ode_orego = ODEProblem(orego,[1.0,2.0,3.0],(0.0,30.0),[77.27,8.375e-6,0.161])

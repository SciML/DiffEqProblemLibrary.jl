srand(100)

### ODE Examples

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

## Lotka-Volterra


lotka = @ode_def_nohes LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d

"""
Lotka-Voltera Equations

```math
\\frac{dx}{dt} = ax - bxy
\\frac{dy}{dt} = -cy + dxy
```

with initial condition ``x=y=1``
"""
prob_ode_lotkavoltera = ODEProblem(lotka,[1.0,1.0],(0.0,1.0),[1.5,1.0,3.0,1.0])

## Fitzhugh-Nagumo

fitz = @ode_def_nohes FitzhughNagumo begin
  dv = v - v^3/3 -w + l
  dw = τinv*(v +  a - b*w)
end a b τinv l
"""
Fitzhugh-Nagumo

```math
\\frac{dv}{dt} = v - \\frac{v^3}{3} - w + I_{est}
τ \\frac{dw}{dt} = v + a -bw
```

with initial condition ``v=w=1``
"""
prob_ode_fitzhughnagumo = ODEProblem(fitz,[1.0;1.0],(0.0,1.0),(0.7,0.8,1/12.5,0.5))

#Van der Pol Equations
van = @ode_def_noinvhes VanDerPol begin
  dy = μ*(1-x^2)*y - x
  dx = 1*y
end μ

"""
Van der Pol Equations

```math
\\begin{align}
\\frac{dx}{dt} &= y \\\\
\\frac{dy}{dt} &= μ(1-x^2)y -x
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

rober = @ode_def_noinvjac Rober begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁ k₂ k₃

"""
The Robertson biochemical reactions:

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
const threebody_μ = parse(BigFloat,"0.012277471"); const threebody_μ′ = 1 - threebody_μ

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
The ThreeBody problem as written by Hairer:

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

Usually solved on `t₀ = 0.0`; `T = parse(BigFloat,"17.0652165601579625588917206249")`
Periodic with that setup.
"""
prob_ode_threebody = ODEProblem(threebody,[0.994, 0.0, 0.0, parse(BigFloat,"-2.00158510637908252240537862224")],(parse(BigFloat,"0.0"),parse(BigFloat,"17.0652165601579625588917206249")))

# Rigid Body Equations

rigid = @ode_def_noinvjac RigidBody begin
  dy₁  = I₁*y₂*y₃
  dy₂  = I₂*y₁*y₃
  dy₃  = I₃*y₁*y₂
end I₁ I₂ I₃

"""
Rigid Body Equations

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
  for i in 14:21
    du[i] = zero(u)
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
Pleiades Problem

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


srand(100)
const mm_A = rand(4,4)
mm_linear = function (du,u,p,t)
  A_mul_B!(du,mm_A,u)
end
const MM_linear =full(Diagonal(0.5ones(4)))
(::typeof(mm_linear))(::Type{Val{:analytic}},u0,p,t) = expm(inv(MM_linear)*mm_A*t)*u0
prob_ode_mm_linear = ODEProblem(mm_linear,rand(4),(0.0,1.0),mass_matrix=MM_linear)

function brusselator_2d_loop(du, u, p, t)
  @inbounds begin
    A, B, α, dx, N = p
    α = α/dx^2
    # Interior
    for i in 2:N-1, j in 2:N-1
      du[i,j,1] = α*(u[i-1,j,1] + u[i+1,j,1] + u[i,j+1,1] + u[i,j-1,1] - 4u[i,j,1]) +
      B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1]
    end
    for i in 2:N-1, j in 2:N-1
      du[i,j,2] = α*(u[i-1,j,2] + u[i+1,j,2] + u[i,j+1,2] + u[i,j-1,2] - 4u[i,j,2]) +
      A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
    end

    # Boundary @ edges
    for j in 2:N-1
      i = 1
      du[1,j,1] = α*(2u[i+1,j,1] + u[i,j+1,1] + u[i,j-1,1] - 4u[i,j,1]) +
      B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1]
    end
    for j in 2:N-1
      i = 1
      du[1,j,2] = α*(2u[i+1,j,2] + u[i,j+1,2] + u[i,j-1,2] - 4u[i,j,2]) +
      A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
    end
    for j in 2:N-1
      i = N
      du[end,j,1] = α*(2u[i-1,j,1] + u[i,j+1,1] + u[i,j-1,1] - 4u[i,j,1]) +
      B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1]
    end
    for j in 2:N-1
      i = N
      du[end,j,2] = α*(2u[i-1,j,2] + u[i,j+1,2] + u[i,j-1,2] - 4u[i,j,2]) +
      A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
    end
    for i in 2:N-1
      j = 1
      du[i,1,1] = α*(u[i-1,j,1] + u[i+1,j,1] + 2u[i,j+1,1] - 4u[i,j,1]) +
      B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1]
    end
    for i in 2:N-1
      j = 1
      du[i,1,2] = α*(u[i-1,j,2] + u[i+1,j,2] + 2u[i,j+1,2] - 4u[i,j,2]) +
      A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
    end
    for i in 2:N-1
      j = N
      du[i,end,1] = α*(u[i-1,j,1] + u[i+1,j,1] + 2u[i,j-1,1] - 4u[i,j,1]) +
      B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1]
    end
    for i in 2:N-1
      j = N
      du[i,end,2] = α*(u[i-1,j,2] + u[i+1,j,2] + 2u[i,j-1,2] - 4u[i,j,2]) +
      A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
    end

    # Boundary @ four vertexes
    i = 1; j = 1
    du[1,1,1] = α*(2u[i+1,j,1] + 2u[i,j+1,1] - 4u[i,j,1]) +
    B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1]
    du[1,1,2] = α*(2u[i+1,j,2] + 2u[i,j+1,2] - 4u[i,j,2]) +
    A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]

    i = 1; j = N
    du[1,N,1] = α*(2u[i+1,j,1] + 2u[i,j-1,1] - 4u[i,j,1]) +
    B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1]
    du[1,N,2] = α*(2u[i+1,j,2] + 2u[i,j-1,2] - 4u[i,j,2]) +
    A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]

    i = N; j = 1
    du[N,1,1] = α*(2u[i-1,j,1] + 2u[i,j+1,1] - 4u[i,j,1]) +
    B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1]
    du[N,1,2] = α*(2u[i-1,j,2] + 2u[i,j+1,2] - 4u[i,j,2]) +
    A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]

    i = N; j = N
    du[end,end,1] = α*(2u[i-1,j,1] + 2u[i,j-1,1] - 4u[i,j,1]) +
    B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1]
    du[end,end,2] = α*(2u[i-1,j,2] + 2u[i,j-1,2] - 4u[i,j,2]) +
    A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
  end
end
function init_brusselator_2d(xyd)
  M = length(xyd)
  u = zeros(M, M, 2)
  for I in CartesianRange((M, M))
    u[I,1] = (2 + 0.25xyd[I[2]])
    u[I,2] = (1 + 0.8xyd[I[1]])
  end
  u
end
const N_brusselator_2d = 128
prob_ode_brusselator_2d = ODEProblem(brusselator_2d_loop,
                                     init_brusselator_2d(linspace(0,1,N_brusselator_2d)),
                                     (0.,1),
                                     (3.4, 1., 0.002, 1/(N_brusselator_2d-1),
                                     N_brusselator_2d))

const N_brusselator_1d = 40
const D_brusselator_u = DerivativeOperator{Float64}(2,2,1/(N_brusselator_1d-1),
                                                    N_brusselator_1d,
                                                    :Dirichlet,:Dirichlet;
                                                    BC=(1.,1.))
const D_brusselator_v = DerivativeOperator{Float64}(2,2,1/(N_brusselator_1d-1),
                                                    N_brusselator_1d,
                                                    :Dirichlet,:Dirichlet;
                                                    BC=(3.,3.))
function brusselator_1d(du, u_, p, t)
    A, B, α, buffer = p
    u = @view(u_[:, 1])
    v = @view(u_[:, 2])
    A_mul_B!(buffer, D_brusselator_u, u)
    Du = buffer
    @. du[:, 1] = A + u^2*v - (B+1)*u + α*Du

    A_mul_B!(buffer, D_brusselator_v, v)
    Dv = buffer
    @. du[:, 2] = B*u - u^2*v + α*Dv
    nothing
end
function init_brusselator_1d(N)
  u = zeros(N, 2)
  x = linspace(0, 1, N)
  for i in 1:N
    u[i, 1] = 1 + sin(2pi*x[i])
    u[i, 2] = 3.
  end
  u
end
prob_ode_brusselator_1d = ODEProblem(brusselator_1d,
                                    init_brusselator_1d(N_brusselator_1d),
                                    (0.,10.),
                                    (1., 3., 1/50, zeros(N_brusselator_1d)))

orego = @ode_def Orego begin
  dy1 = p1*(y2+y1*(1-p2*y1-y2))
  dy2 = (y3-(1+y1)*y2)/p1
  dy3 = p3*(y1-y3)
end p1 p2 p3
prob_ode_orego = ODEProblem(orego,[1.0,2.0,3.0],(0.0,30.0),[77.27,8.375e-6,0.161])

hires = @ode_def Hires begin
  dy1 = -1.71*y1 + 0.43*y2 + 8.32*y3 + 0.0007
  dy2 = 1.71*y1 - 8.75*y2
  dy3 = -10.03*y3 + 0.43*y4 + 0.035*y5
  dy4 = 8.32*y2 + 1.71*y3 - 1.12*y4
  dy5 = -1.745*y5 + 0.43*y6 + 0.43*y7
  dy6 = -280.0*y6*y8 + 0.69*y4 + 1.71*y5 -
           0.43*y6 + 0.69*y7
  dy7 = 280.0*y6*y8 - 1.81*y7
  dy8 = -280.0*y6*y8 + 1.81*y7
end
u0 = zeros(8)
u0[1] = 1
u0[8] = 0.0057
prob_ode_hires = ODEProblem(hires,u0,(0.0,321.8122))

const k1=.35e0
const k2=.266e2
const k3=.123e5
const k4=.86e-3
const k5=.82e-3
const k6=.15e5
const k7=.13e-3
const k8=.24e5
const k9=.165e5
const k10=.9e4
const k11=.22e-1
const k12=.12e5
const k13=.188e1
const k14=.163e5
const k15=.48e7
const k16=.35e-3
const k17=.175e-1
const k18=.1e9
const k19=.444e12
const k20=.124e4
const k21=.21e1
const k22=.578e1
const k23=.474e-1
const k24=.178e4
const k25=.312e1

function pollution(dy,y,p,t)
 r1  = k1 *y[1]
 r2  = k2 *y[2]*y[4]
 r3  = k3 *y[5]*y[2]
 r4  = k4 *y[7]
 r5  = k5 *y[7]
 r6  = k6 *y[7]*y[6]
 r7  = k7 *y[9]
 r8  = k8 *y[9]*y[6]
 r9  = k9 *y[11]*y[2]
 r10 = k10*y[11]*y[1]
 r11 = k11*y[13]
 r12 = k12*y[10]*y[2]
 r13 = k13*y[14]
 r14 = k14*y[1]*y[6]
 r15 = k15*y[3]
 r16 = k16*y[4]
 r17 = k17*y[4]
 r18 = k18*y[16]
 r19 = k19*y[16]
 r20 = k20*y[17]*y[6]
 r21 = k21*y[19]
 r22 = k22*y[19]
 r23 = k23*y[1]*y[4]
 r24 = k24*y[19]*y[1]
 r25 = k25*y[20]

 dy[1]  = -r1-r10-r14-r23-r24+
          r2+r3+r9+r11+r12+r22+r25
 dy[2]  = -r2-r3-r9-r12+r1+r21
 dy[3]  = -r15+r1+r17+r19+r22
 dy[4]  = -r2-r16-r17-r23+r15
 dy[5]  = -r3+r4+r4+r6+r7+r13+r20
 dy[6]  = -r6-r8-r14-r20+r3+r18+r18
 dy[7]  = -r4-r5-r6+r13
 dy[8]  = r4+r5+r6+r7
 dy[9]  = -r7-r8
 dy[10] = -r12+r7+r9
 dy[11] = -r9-r10+r8+r11
 dy[12] = r9
 dy[13] = -r11+r10
 dy[14] = -r13+r12
 dy[15] = r14
 dy[16] = -r18-r19+r16
 dy[17] = -r20
 dy[18] = r20
 dy[19] = -r21-r22-r24+r23+r25
 dy[20] = -r25+r24
end

function pollution(::Type{Val{:jac}},J,y,p,t)
  J .= 0.0
  J[1,1]   = -k1-k10*y[11]-k14*y[6]-k23*y[4]-k24*y[19]
  J[1,11]  = -k10*y[1]+k9*y[2]
  J[1,6]   = -k14*y[1]
  J[1,4]   = -k23*y[1]+k2*y[2]
  J[1,19]  = -k24*y[1]+k22
  J[1,2]   = k2*y[4]+k9*y[11]+k3*y[5]+k12*y[10]
  J[1,13]  = k11
  J[1,20]  = k25
  J[1,5]   = k3*y[2]
  J[1,10]  = k12*y[2]

  J[2,4]   = -k2*y[2]
  J[2,5]   = -k3*y[2]
  J[2,11]  = -k9*y[2]
  J[2,10]  = -k12*y[2]
  J[2,19]  = k21
  J[2,1]   = k1
  J[2,2]   = -k2*y[4]-k3*y[5]-k9*y[11]-k12*y[10]

  J[3,1]   = k1
  J[3,4]   = k17
  J[3,16]  = k19
  J[3,19]  = k22
  J[3,3]   = -k15

  J[4,4]   = -k2*y[2]-k16-k17-k23*y[1]
  J[4,2]   = -k2*y[4]
  J[4,1]   = -k23*y[4]
  J[4,3]   = k15

  J[5,5]   = -k3*y[2]
  J[5,2]   = -k3*y[5]
  J[5,7]   = 2k4+k6*y[6]
  J[5,6]   = k6*y[7]+k20*y[17]
  J[5,9]   = k7
  J[5,14]  = k13
  J[5,17]  = k20*y[6]

  J[6,6]   = -k6*y[7]-k8*y[9]-k14*y[1]-k20*y[17]
  J[6,7]   = -k6*y[6]
  J[6,9]   = -k8*y[6]
  J[6,1]   = -k14*y[6]
  J[6,17]  = -k20*y[6]
  J[6,2]   = k3*y[5]
  J[6,5]   = k3*y[2]
  J[6,16]  = 2k18

  J[7,7]   = -k4-k5-k6*y[6]
  J[7,6]   = -k6*y[7]
  J[7,14]  = k13

  J[8,7]   = k4+k5+k6*y[6]
  J[8,6]   = k6*y[7]
  J[8,9]   = k7

  J[9,9]   = -k7-k8*y[6]
  J[9,6]   = -k8*y[9]

  J[10,10] = -k12*y[2]
  J[10,2]  = -k12*y[10]+k9*y[11]
  J[10,9]  = k7
  J[10,11] = k9*y[2]

  J[11,11] = -k9*y[2]-k10*y[1]
  J[11,2]  = -k9*y[11]
  J[11,1]  = -k10*y[11]
  J[11,9]  = k8*y[6]
  J[11,6]  = k8*y[9]
  J[11,13] = k11

  J[12,11] = k9*y[2]
  J[12,2]  = k9*y[11]

  J[13,13] = -k11
  J[13,11] = k10*y[1]
  J[13,1]  = k10*y[11]

  J[14,14] = -k13
  J[14,10] = k12*y[2]
  J[14,2]  = k12*y[10]

  J[15,1]  = k14*y[6]
  J[15,6]  = k14*y[1]

  J[16,16] = -k18-k19
  J[16,4]  = k16

  J[17,17] = -k20*y[6]
  J[17,6]  = -k20*y[17]

  J[18,17] = k20*y[6]
  J[18,6]  = k20*y[17]

  J[19,19] = -k21-k22-k24*y[1]
  J[19,1]  = -k24*y[19]+k23*y[4]
  J[19,4]  = k23*y[1]
  J[19,20] = k25

  J[20,20] = -k25
  J[20,1]  = k24*y[19]
  J[20,19] = k24*y[1]

  return nothing
end
u0 = zeros(20)
u0[2]  = 0.2
u0[4]  = 0.04
u0[7]  = 0.1
u0[8]  = 0.3
u0[9]  = 0.01
u0[17] = 0.007
prob_ode_pollution = ODEProblem(pollution,u0,(0.0,60.0))

######################################################
# Filament
######################################################

const T = Float64
abstract type AbstractFilamentCache end
abstract type AbstractMagneticForce end
abstract type AbstractInextensibilityCache end
abstract type AbstractSolver end
abstract type AbstractSolverCache end
struct FerromagneticContinuous <: AbstractMagneticForce
    ω :: T
    F :: Vector{T}
end

mutable struct FilamentCache{
        MagneticForce        <: AbstractMagneticForce,
        InextensibilityCache <: AbstractInextensibilityCache,
        SolverCache          <: AbstractSolverCache
            } <: AbstractFilamentCache
    N  :: Int
    μ  :: T
    Cm :: T
    x  :: SubArray{T,1,Vector{T},Tuple{StepRange{Int,Int}},true}
    y  :: SubArray{T,1,Vector{T},Tuple{StepRange{Int,Int}},true}
    z  :: SubArray{T,1,Vector{T},Tuple{StepRange{Int,Int}},true}
    A  :: Matrix{T}
    P  :: InextensibilityCache
    F  :: MagneticForce
    Sc :: SolverCache
end
struct NoHydroProjectionCache <: AbstractInextensibilityCache
    J         :: Matrix{T}
    P         :: Matrix{T}
    J_JT      :: Matrix{T}
    J_JT_LDLT :: Base.LinAlg.LDLt{T, SymTridiagonal{T}}
    P0        :: Matrix{T}

    NoHydroProjectionCache(N::Int) = new(
        zeros(N, 3*(N+1)),          # J
        zeros(3*(N+1), 3*(N+1)),    # P
        zeros(N,N),                 # J_JT
        Base.LinAlg.LDLt{T,SymTridiagonal{T}}(SymTridiagonal(zeros(N), zeros(N-1))),
        zeros(N, 3*(N+1))
    )
end
struct DiffEqSolverCache <: AbstractSolverCache
    S1 :: Vector{T}
    S2 :: Vector{T}

    DiffEqSolverCache(N::Integer) = new(zeros(T,3*(N+1)), zeros(T,3*(N+1)))
end
function FilamentCache(N=20; Cm=32, ω=200, Solver=SolverDiffEq)
    InextensibilityCache = NoHydroProjectionCache
    SolverCache = DiffEqSolverCache
    tmp = zeros(3*(N+1))
    FilamentCache{FerromagneticContinuous, InextensibilityCache, SolverCache}(
        N, N+1, Cm, view(tmp,1:3:3*(N+1)), view(tmp,2:3:3*(N+1)), view(tmp,3:3:3*(N+1)),
        zeros(3*(N+1), 3*(N+1)), # A
        InextensibilityCache(N), # P
        FerromagneticContinuous(ω, zeros(3*(N+1))),
        SolverCache(N)
    )
end
function stiffness_matrix!(f::AbstractFilamentCache)
    N, μ, A = f.N, f.μ, f.A
    A[:] = eye(3*(N+1))
    for i in 1 : 3
        A[i,i] =    1
        A[i,3+i] = -2
        A[i,6+i] =  1

        A[3+i,i]   = -2
        A[3+i,3+i] =  5
        A[3+i,6+i] = -4
        A[3+i,9+i] =  1

        A[3*(N-1)+i,3*(N-3)+i] =  1
        A[3*(N-1)+i,3*(N-2)+i] = -4
        A[3*(N-1)+i,3*(N-1)+i] =  5
        A[3*(N-1)+i,3*N+i]     = -2

        A[3*N+i,3*(N-2)+i]     =  1
        A[3*N+i,3*(N-1)+i]     = -2
        A[3*N+i,3*N+i]         =  1

        for j in 2 : N-2
            A[3*j+i,3*j+i]     =  6
            A[3*j+i,3*(j-1)+i] = -4
            A[3*j+i,3*(j+1)+i] = -4
            A[3*j+i,3*(j-2)+i] =  1
            A[3*j+i,3*(j+2)+i] =  1
        end
    end
    scale!(A, -μ^4)
    nothing
end
function update_separate_coordinates!(f::AbstractFilamentCache, r)
    N, x, y, z = f.N, f.x, f.y, f.z
    @inbounds for i in 1 : length(x)
        x[i] = r[3*i-2]
        y[i] = r[3*i-1]
        z[i] = r[3*i]
    end
    nothing
end

function update_united_coordinates!(f::AbstractFilamentCache, r)
    N, x, y, z = f.N, f.x, f.y, f.z
    @inbounds for i in 1 : length(x)
        r[3*i-2] = x[i]
        r[3*i-1] = y[i]
        r[3*i]   = z[i]
    end
    nothing
end

function update_united_coordinates(f::AbstractFilamentCache)
    r = zeros(T, 3*length(f.x))
    update_united_coordinates!(f, r)
    r
end

function initialize!(initial_conf_type::Symbol, f::AbstractFilamentCache)
    N, x, y, z = f.N, f.x, f.y, f.z
    if initial_conf_type == :StraightX
        x[:] = linspace(0, 1, N+1)
        y[:] = 0 .* x
        z[:] = 0 .* x
    else
        error("Unknown initial configuration requested.")
    end
    update_united_coordinates(f)
end

function magnetic_force!(::FerromagneticContinuous, f::AbstractFilamentCache, t)
    # TODO: generalize this for different magnetic fields as well
    N, μ, Cm, ω, F = f.N, f.μ, f.Cm, f.F.ω, f.F.F
    F[1]         = -μ * Cm * cos(ω*t)
    F[2]         = -μ * Cm * sin(ω*t)
    F[3*(N+1)-2] =  μ * Cm * cos(ω*t)
    F[3*(N+1)-1] =  μ * Cm * sin(ω*t)
    nothing
end

struct SolverDiffEq <: AbstractSolver end

function (f::FilamentCache)(dr, r, p, t)
    @views f.x, f.y, f.z = r[1:3:end], r[2:3:end], r[3:3:end]
    jacobian!(f)
    projection!(f)
    magnetic_force!(f.F, f, t)
    A, P, F, S1, S2 = f.A, f.P.P, f.F.F, f.Sc.S1, f.Sc.S2

    # implement dr = P * (A*r + F) in an optimized way to avoid temporaries
    A_mul_B!(S1, A, r)
    S1 .+= F
    A_mul_B!(S2, P, S1)
    copy!(dr, S2)
    return dr
end

function (f::FilamentCache)(::Type{Val{:jac}}, J, r, p, t)
    A_mul_B!(J, f.P.P, f.A)
    nothing
end

function jacobian!(f::FilamentCache)
    N, x, y, z, J = f.N, f.x, f.y, f.z, f.P.J
    @inbounds for i in 1 : N
        J[i, 3*i-2]     = -2 * (x[i+1]-x[i])
        J[i, 3*i-1]     = -2 * (y[i+1]-y[i])
        J[i, 3*i]       = -2 * (z[i+1]-z[i])
        J[i, 3*(i+1)-2] =  2 * (x[i+1]-x[i])
        J[i, 3*(i+1)-1] =  2 * (y[i+1]-y[i])
        J[i, 3*(i+1)]   =  2 * (z[i+1]-z[i])
    end
    nothing
end

function projection!(f::FilamentCache)
    # implement P[:] = I - J'/(J*J')*J in an optimized way to avoid temporaries
    J, P, J_JT, J_JT_LDLT, P0 = f.P.J, f.P.P, f.P.J_JT, f.P.J_JT_LDLT, f.P.P0
    A_mul_Bt!(J_JT, J, J)
    LDLt_inplace!(J_JT_LDLT, J_JT)
    A_ldiv_B!(P0, J_JT_LDLT, J)
    At_mul_B!(P, P0, J)
    subtract_from_identity!(P)
    nothing
end

function subtract_from_identity!(A)
    scale!(-1, A)
    @inbounds for i in 1 : size(A,1)
        A[i,i] += 1
    end
    nothing
end

function LDLt_inplace!{T<:Real}(L::Base.LinAlg.LDLt{T,SymTridiagonal{T}}, A::Matrix{T})
    n = size(A,1)
    dv, ev = L.data.dv, L.data.ev
    @inbounds for (i,d) in enumerate(diagind(A))
        dv[i] = A[d]
    end
    @inbounds for (i,d) in enumerate(diagind(A,-1))
        ev[i] = A[d]
    end
    @inbounds @simd for i in 1 : n-1
        ev[i]   /= dv[i]
        dv[i+1] -= abs2(ev[i]) * dv[i]
    end
    L
end

function filament_prob(::SolverDiffEq; N=20, Cm=32, ω=200, time_end=1.)
    f = FilamentCache(N, Solver=SolverDiffEq, Cm=Cm, ω=ω)
    r0 = initialize!(:StraightX, f)
    stiffness_matrix!(f)
    prob = ODEProblem(f, r0, (0., time_end))
end
prob_ode_filament = filament(SolverDiffEq())

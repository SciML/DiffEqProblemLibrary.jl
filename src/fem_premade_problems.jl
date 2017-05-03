### Finite Element Examples

analytic_moving(t,x) = 0.1*(1-exp.(-100*(t-0.5).^2)).*exp.(-25((x[:,1]-t+0.5).^2 + (x[:,2]-t+0.5).^2))
Du = (t,x) -> -50[analytic_moving(t,x).*(0.5-t+x[:,1])  analytic_moving(t,x).*(0.5-t+x[:,2])]
f = (t,x) -> (-5).*exp.((-25).*((3/2)+6.*t.^2+x[:,1]+x[:,1].^2+x[:,2]+x[:,2].^2+(-2).*t.*(3+x[:,1]+
  x[:,2]))).*((-20)+(-100).*t.^2+(-49).*x[:,1]+(-50).*x[:,1].^2+(-49).*x[:,2]+(-50).*
  x[:,2].^2+2.*t.*(47+50.*x[:,1]+50.*x[:,2])+exp.(25.*(1+(-2).*t).^2).*(22+
  100.*t.^2+49.*x[:,1]+50.*x[:,1].^2+49.*x[:,2]+50.*x[:,2].^2+(-2).*t.*(49+50.*x[:,1]+50.*x[:,2])))
tspan = (0.0,2.0)
dx = 1//2^(3)
dt = 1//2^(9)
mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:dirichlet)
"""
Example problem defined by the solution:
```math
u(x,y,t)=\\frac{1}{10}(1-\\exp(-100(t-\\frac{1}{2})^2))\\exp(-25((x-t+0.5)^2 + (y-t+0.5)^2))
```

This will have a mound which moves across the screen. Good animation test.
"""
prob_femheat_moving = HeatProblem(analytic_moving,Du,f,mesh)

N = 2 #Number of different dt to solve at, 2 for test speed
topdt = 6 # 1//2^(topdt-1) is the max dt. Small for test speed
dts = 1.//2.^(topdt-1:-1:N)
dxs = 1//2^(5) * ones(dts) #Run at 2^-7 for best plot
probs = [HeatProblem(analytic_moving,Du,f,parabolic_squaremesh([0 1 0 1],dxs[i],dts[i],(0.0,1.0),:dirichlet)) for i in eachindex(dts)]
cs_femheat_moving_dt = ConvergenceSetup(probs,dts)
dxs = 1//2^(4) * ones(dts) #Run at 2^-7 for best plot
probs = [HeatProblem(analytic_moving,Du,f,parabolic_squaremesh([0 1 0 1],dxs[i],dts[i],(0.0,1.0),:dirichlet)) for i in eachindex(dts)]
cs_femheat_moving_faster_dt = ConvergenceSetup(probs,dts)

#Not good plots, but quick for unit tests
dxs = 1.//2.^(2:-1:1)
dts = 1//2^(6) * ones(dxs) #Run at 2^-7 for best plot
probs = [HeatProblem(analytic_moving,Du,f,parabolic_squaremesh([0 1 0 1],dxs[i],dts[i],(0.0,1.0),:dirichlet)) for i in eachindex(dts)]
cs_femheat_moving_dx = ConvergenceSetup(probs,dxs)

tspan = (0.0,1.0)
dx = 1//2^(3)
dt = 1//2^(7)
mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:dirichlet)
"""
Example problem defined by the solution:
```math
u(x,y,t)=\\frac{1}{10}(1-\\exp(-100(t-\\frac{1}{2})^2))\\exp(-25((x-t+0.5)^2 + (y-t+0.5)^2))
```

This will have a mound which moves across the screen. Good animation test.
"""
prob_femheat_moving7 = HeatProblem(analytic_moving,Du,f,mesh)

analytic_diffuse(t,x) = exp.(-10((x[:,1]-.5).^2 + (x[:,2]-.5).^2 )-t)
f = (t,x) -> exp.(-t-5*(1-2x[:,1]+2x[:,1].^2 - 2x[:,2] +2x[:,2].^2)).*(-161 + 400*(x[:,1] - x[:,1].^2 + x[:,2] - x[:,2].^2))
Du = (t,x) -> -20[analytic_diffuse(t,x).*(x[:,1]-.5) analytic_diffuse(t,x).*(x[:,2]-.5)]
tspan = (0.0,1.0)
dx = 1//2^(3)
dt = 1//2^(7)
mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:dirichlet)
"""
Example problem defined by the solution:

```math
u(x,y,t)=\\exp(-10((x-\\frac{1}{2})^2 + (y-\\frac{1}{2})^2 )-t)
```

This is a Gaussian centered at ``(\\frac{1}{2},\\frac{1}{2})`` which diffuses over time.
"""
prob_femheat_diffuse = HeatProblem(analytic_diffuse,Du,f,mesh)


f = (t,x)  -> zeros(size(x,1))
if VERSION < v"0.6-"
  u0_func = (x) -> float((abs.(x[:,1]-.5) .< 1e-6) & (abs.(x[:,2]-.5) .< 1e-6)) #Only mass at middle of (0.0,1.0)^2
else
  u0_func = (x) -> float((abs.(x[:,1]-.5) .< 1e-6) .& (abs.(x[:,2]-.5) .< 1e-6)) #Only mass at middle of (0.0,1.0)^2
end
tspan = (0.0,1/2^(5))
dx = 1//2^(3)
dt = 1//2^(9)
mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:dirichlet)
u0 = u0_func(mesh.node)
"""
Example problem which starts with a Dirac δ cenetered at (0.5,0.5) and solves with ``f=gD=0``.
This gives the Green's function solution.
"""
prob_femheat_pure = HeatProblem(u0,f,mesh)

tspan =  (0.0,1/2^(5))
dt = 1//2^(11)
mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:dirichlet)
u0 = u0_func(mesh.node)
"""
Example problem which starts with a Dirac δ cenetered at (0.5,0.5) and solves with ``f=gD=0``.
This gives the Green's function solution.
"""
prob_femheat_pure11 = HeatProblem(u0,f,mesh)

f = (t,x,u) -> ones(size(x,1)) - .5u
u0_func = (x) -> zeros(size(x,1))
tspan = (0.0,1.0)
dx = 1//2^(3)
dt = 1//2^(7)
mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:neumann)
u0 = u0_func(mesh.node)
"""
Homogenous reaction-diffusion problem which starts with 0 and solves with ``f(u)=1-u/2``
"""
prob_femheat_birthdeath = HeatProblem(u0,f,mesh)


f = (t,x,u)  -> [ones(size(x,1))-.5u[:,1]   ones(size(x,1))-u[:,2]]
u0_func = (x) -> ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
tspan = (0.0,5.0)
dx = 1/2^(1)
dt = 1/2^(7)
mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:neumann)
u0 = u0_func(mesh.node)
"""
Homogenous reaction-diffusion which starts at 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=1-v``
"""
prob_femheat_birthdeathsystem = HeatProblem(u0,f,mesh)

f  = (t,x,u)  -> [ones(size(x,1))-.5u[:,1]     .5u[:,1]-u[:,2]]
u0_func = (x) -> ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
tspan = (0.0,5.0)
dx = 1/2^(1)
dt = 1/2^(7)
mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:neumann)
u0 = u0_func(mesh.node)
"""
Homogenous reaction-diffusion which starts with 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=.5u-v``
"""
prob_femheat_birthdeathinteractingsystem = HeatProblem(u0,f,mesh)

#=
f = (t,x,u)  -> [zeros(size(x,1))    zeros(size(x,1))]
u0 = (x) -> [float((abs.(x[:,1]-.5) .< 1e-6) .& (abs.(x[:,2]-.5) .< 1e-6)) float((abs.(x[:,1]-.5) .< 1e-6) .& (abs.(x[:,2]-.5) .< 1e-6))]  # size (x,2), 2 meaning 2 variables
"""
Example problem which solves the homogeneous Heat equation with all mass starting at (1/2,1/2) with two different diffusion constants,
``D₁=0.01`` and ``D₂=0.001``. Good animation test.
"""
prob_femheat_diffusionconstants = HeatProblem(u0,f,mesh,D=[.01 .001])
=#

#=
"""
`heatProblemExample_grayscott(;ρ=.03,k=.062,D=[1e-3 .5e-3])`

The Gray-Scott equations with quasi-random initial conditions. The reaction equations
are given by:

```math
\\begin{align}
u_t &= uv^2 + ρ(1-v) \\\\
v_t &= uv^2 - (ρ+k)v
\\end{align}
```
"""
function heatProblemExample_grayscott(;ρ=.03,k=.062,D=[1e-3 .5e-3])
  f₁(t,x,u)  = u[:,1].*u[:,2].*u[:,2] + ρ*(1-u[:,2])
  f₂(t,x,u)  = u[:,1].*u[:,2].*u[:,2] -(ρ+k).*u[:,2]
  f(t,x,u) = [f₁(t,x,u) f₂(t,x,u)]
  u0(x) = [ones(size(x,1))+rand(size(x,1)) .25.*float(((.2.<x[:,1].<.6) .&
          (.2.<x[:,2].<.6)) | ((.85.<x[:,1]) .& (.85.<x[:,2])))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u0,f,D=D))
end

"""
`heatProblemExample_gierermeinhardt(;a=1,α=1,D=[0.01 1.0],ubar=1,vbar=0,β=10,startNoise=0.01)`

The Gierer-Meinhardt equations wtih quasi-random initial perturbations.

```math
\\begin{align}
u_t &= \\frac{au}{v^2} + \\bar{u} -αu \\\\
v_t &= au^2 + \\bar{v} - βv
\\end{align}
```

The equation starts at the steady state

```math
\\begin{align}
u_{ss} &= \\frac{\\bar{u}+β}{α} \\\\
v_{ss} &= \\frac{α}{β} u_{ss}^2
\\end{align}
```

with a bit of noise.

"""
function heatProblemExample_gierermeinhardt(;a=1,α=1,D=[0.01 1.0],ubar=1,vbar=0,β=10,startNoise=0.01)
  f₁(t,x,u)  = a*u[:,1].*u[:,1]./u[:,2] + ubar - α*u[:,1]
  f₂(t,x,u)  = a*u[:,1].*u[:,1] + vbar -β.*u[:,2]
  f(t,x,u) = [f₁(t,x,u) f₂(t,x,u)]
  uss = (ubar +β)/α
  vss = (α/β)*uss.^2
  u0(x) = [uss*ones(size(x,1))+startNoise*rand(size(x,1)) vss*ones(size(x,1))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u0,f,D=D))
end
=#

f = (t,x,u)  -> ones(size(x,1)) - .5u
u0_func = (x) -> zeros(size(x,1))
σ = (t,x,u) -> 1u.^2
tspan =  (0.0,5.0)
dx = 1//2^(3)
dt = 1//2^(5)
mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:neumann)
u0 = u0_func(mesh.node)
"""
Homogenous stochastic reaction-diffusion problem which starts with 0
and solves with ``f(u)=1-u/2`` with noise ``σ(u)=10u^2``
"""
prob_femheat_stochasticbirthdeath = HeatProblem(u0,f,mesh,σ=σ)

dx = 1//2^(1)
dt = 1//2^(1)
mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:neumann)
u0 = u0_func(mesh.node)
prob_femheat_stochasticbirthdeath_fast = HeatProblem(u0,f,mesh,σ=σ)

## Poisson

f = (x) -> sin.(2π.*x[:,1]).*cos.(2π.*x[:,2])
analytic = (x) -> sin.(2π.*x[:,1]).*cos.(2π.*x[:,2])/(8π*π)
Du = (x) -> [cos.(2*pi.*x[:,1]).*cos.(2*pi.*x[:,2])./(4*pi) -sin.(2π.*x[:,1]).*sin.(2π.*x[:,2])./(4π)]
dx = 1//2^(5)
mesh = notime_squaremesh([0 1 0 1],dx,:dirichlet)
"""
Problem defined by the solution: ``u(x,y)= \\sin(2πx)\\cos(2πy)/(8π^2)``
"""
prob_poisson_wave = PoissonProblem(f,analytic,Du,mesh)

dxs = 1.//2.^(4:-1:2) # 4 for testing, use 7 for good graph
probs = [PoissonProblem(f,analytic,Du,notime_squaremesh([0 1 0 1],dx,:dirichlet)) for dx in dxs]
cs_fempoisson_wave = ConvergenceSetup(probs,dts)

σ = (x) -> 5 #Additive noise
dx = 1//2^(5)
mesh = notime_squaremesh([0 1 0 1],dx,:dirichlet)
"""
Problem with deterministic solution: ``u(x,y)= \\sin(2πx)\\cos(2πy)/(8π^2)``
and additive noise ``σ(x,y)=5``
"""
prob_poisson_noisywave = PoissonProblem(f,analytic,Du,mesh,σ=σ)

f = (x,u) -> ones(size(x,1)) - .5u
dx = 1//2^(3)
mesh = notime_squaremesh([0 1 0 1],dx,:neumann)
"""
Nonlinear Poisson equation with ``f(u)=1-u/2``.
Corresponds to the steady state of a humogenous reaction-diffusion equation
with the same ``f``.
"""
prob_poisson_birthdeath = PoissonProblem(f,mesh,numvars=1)

f  = (x,u) -> [ones(size(x,1))-.5u[:,1]     ones(size(x,1))-u[:,2]]
u0_func = (x) -> .5*ones(size(x,1),2) # size (x,2), 2 meaning 2 variables
dx = 1//2^(1)
mesh = notime_squaremesh([0 1 0 1],dx,:neumann)
u0 = u0_func(mesh.node)
"""
Nonlinear Poisson equation with ``f(u)=1-u/2`` and ``f(v)=1-v`` and initial
condition homogenous 1/2. Corresponds to the steady state of a humogenous
reaction-diffusion equation with the same ``f``.
"""
prob_poisson_birthdeathsystem = PoissonProblem(f,mesh,u0=u0)

f  = (x,u) -> [ones(size(x,1))-.5u[:,1]     .5u[:,1]-u[:,2]]
u0_func = (x) -> ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
dx = 1//2^(1)
mesh = notime_squaremesh([0 1 0 1],dx,:neumann)
u0 = u0_func(mesh.node)
"""
Nonlinear Poisson equation with ``f(u)=1-u/2`` and ``f(v)=.5u-v`` and initial
condition homogenous 1/2. Corresponds to the steady state of a humogenous
reaction-diffusion equation with the same ``f``.
"""
prob_poisson_birthdeathinteractingsystem = PoissonProblem(f,mesh,u0=u0)

### Finite Element Examples

analytic_moving(t,x) = 0.1*(1-exp.(-100*(t-0.5).^2)).*exp.(-25((x[:,1]-t+0.5).^2 + (x[:,2]-t+0.5).^2))
Du = (t,x) -> -50[analytic_moving(t,x).*(0.5-t+x[:,1])  analytic_moving(t,x).*(0.5-t+x[:,2])]
f = (t,x) -> (-5).*exp.((-25).*((3/2)+6.*t.^2+x[:,1]+x[:,1].^2+x[:,2]+x[:,2].^2+(-2).*t.*(3+x[:,1]+
  x[:,2]))).*((-20)+(-100).*t.^2+(-49).*x[:,1]+(-50).*x[:,1].^2+(-49).*x[:,2]+(-50).*
  x[:,2].^2+2.*t.*(47+50.*x[:,1]+50.*x[:,2])+exp.(25.*(1+(-2).*t).^2).*(22+
  100.*t.^2+49.*x[:,1]+50.*x[:,1].^2+49.*x[:,2]+50.*x[:,2].^2+(-2).*t.*(49+50.*x[:,1]+50.*x[:,2])))
"""
Example problem defined by the solution:
```math
u(x,y,t)=\\frac{1}{10}(1-\\exp(-100(t-\\frac{1}{2})^2))\\exp(-25((x-t+0.5)^2 + (y-t+0.5)^2))
```

This will have a mound which moves across the screen. Good animation test.
"""
prob_femheat_moving = HeatProblem(analytic_moving,Du,f)



analytic_diffuse(t,x) = exp.(-10((x[:,1]-.5).^2 + (x[:,2]-.5).^2 )-t)
f = (t,x) -> exp.(-t-5*(1-2x[:,1]+2x[:,1].^2 - 2x[:,2] +2x[:,2].^2)).*(-161 + 400*(x[:,1] - x[:,1].^2 + x[:,2] - x[:,2].^2))
Du = (t,x) -> -20[analytic_diffuse(t,x).*(x[:,1]-.5) analytic_diffuse(t,x).*(x[:,2]-.5)]
"""
Example problem defined by the solution:

```math
u(x,y,t)=\\exp(-10((x-\\frac{1}{2})^2 + (y-\\frac{1}{2})^2 )-t)
```

This is a Gaussian centered at ``(\\frac{1}{2},\\frac{1}{2})`` which diffuses over time.
"""
prob_femheat_diffuse = HeatProblem(analytic_diffuse,Du,f)


f = (t,x)  -> zeros(size(x,1))
u₀ = (x) -> float((abs.(x[:,1]-.5) .< 1e-6) & (abs.(x[:,2]-.5) .< 1e-6)) #Only mass at middle of (0,1)^2
"""
Example problem which starts with a Dirac δ cenetered at (0.5,0.5) and solves with ``f=gD=0``.
This gives the Green's function solution.
"""
prob_femheat_pure = HeatProblem(u₀,f)


f = (t,x,u) -> ones(size(x,1)) - .5u
u₀ = (x) -> zeros(size(x,1))
"""
Homogenous reaction-diffusion problem which starts with 0 and solves with ``f(u)=1-u/2``
"""
prob_femheat_birthdeath = HeatProblem(u₀,f)


f = (t,x,u)  -> [ones(size(x,1))-.5u[:,1]   ones(size(x,1))-u[:,2]]
u₀ = (x) -> ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
"""
Homogenous reaction-diffusion which starts at 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=1-v``
"""
prob_femheat_birthdeathsystem = HeatProblem(u₀,f)

f = (t,x,u)  -> [zeros(size(x,1))    zeros(size(x,1))]
u₀ = (x) -> [float((abs.(x[:,1]-.5) .< 1e-6) & (abs.(x[:,2]-.5) .< 1e-6)) float((abs.(x[:,1]-.5) .< 1e-6) & (abs.(x[:,2]-.5) .< 1e-6))]  # size (x,2), 2 meaning 2 variables
"""
Example problem which solves the homogeneous Heat equation with all mass starting at (1/2,1/2) with two different diffusion constants,
``D₁=0.01`` and ``D₂=0.001``. Good animation test.
"""
prob_femheat_diffusionconstants = HeatProblem(u₀,f,D=[.01 .001])

f  = (t,x,u)  -> [ones(size(x,1))-.5u[:,1]     .5u[:,1]-u[:,2]]
u₀ = (x) -> ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
"""
Homogenous reaction-diffusion which starts with 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=.5u-v``
"""
prob_femheat_birthdeathinteractingsystem = HeatProblem(u₀,f)

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
  u₀(x) = [ones(size(x,1))+rand(size(x,1)) .25.*float(((.2.<x[:,1].<.6) &
          (.2.<x[:,2].<.6)) | ((.85.<x[:,1]) & (.85.<x[:,2])))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
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
  u₀(x) = [uss*ones(size(x,1))+startNoise*rand(size(x,1)) vss*ones(size(x,1))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
end

f = (t,x,u)  -> ones(size(x,1)) - .5u
u₀ = (x) -> zeros(size(x,1))
σ = (t,x,u) -> 1u.^2
"""
Homogenous stochastic reaction-diffusion problem which starts with 0
and solves with ``f(u)=1-u/2`` with noise ``σ(u)=10u^2``
"""
prob_femheat_stochasticbirthdeath = HeatProblem(u₀,f,σ=σ)

## Poisson

f = (x) -> sin.(2π.*x[:,1]).*cos.(2π.*x[:,2])
analytic = (x) -> sin.(2π.*x[:,1]).*cos.(2π.*x[:,2])/(8π*π)
Du = (x) -> [cos.(2*pi.*x[:,1]).*cos.(2*pi.*x[:,2])./(4*pi) -sin.(2π.*x[:,1]).*sin.(2π.*x[:,2])./(4π)]
"""
Problem defined by the solution: ``u(x,y)= \\sin(2πx)\\cos(2πy)/(8π^2)``
"""
prob_poisson_wave = PoissonProblem(f,analytic,Du)

σ = (x) -> 5 #Additive noise
"""
Problem with deterministic solution: ``u(x,y)= \\sin(2πx)\\cos(2πy)/(8π^2)``
and additive noise ``σ(x,y)=5``
"""
prob_poisson_noisywave = PoissonProblem(f,analytic,Du,σ=σ)

f = (x,u) -> ones(size(x,1)) - .5u
"""
Nonlinear Poisson equation with ``f(u)=1-u/2``.
Corresponds to the steady state of a humogenous reaction-diffusion equation
with the same ``f``.
"""
prob_poisson_birthdeath = PoissonProblem(f,numvars=1)

f  = (x,u) -> [ones(size(x,1))-.5u[:,1]     ones(size(x,1))-u[:,2]]
u₀ = (x) -> .5*ones(size(x,1),2) # size (x,2), 2 meaning 2 variables

"""
Nonlinear Poisson equation with ``f(u)=1-u/2`` and ``f(v)=1-v`` and initial
condition homogenous 1/2. Corresponds to the steady state of a humogenous
reaction-diffusion equation with the same ``f``.
"""
prob_poisson_birthdeathsystem = PoissonProblem(f,u₀=u₀)

f  = (x,u) -> [ones(size(x,1))-.5u[:,1]     .5u[:,1]-u[:,2]]
u₀ = (x) -> ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables

"""
Nonlinear Poisson equation with ``f(u)=1-u/2`` and ``f(v)=.5u-v`` and initial
condition homogenous 1/2. Corresponds to the steady state of a humogenous
reaction-diffusion equation with the same ``f``.
"""
prob_poisson_birthdeathinteractingsystem = PoissonProblem(f,u₀=u₀)

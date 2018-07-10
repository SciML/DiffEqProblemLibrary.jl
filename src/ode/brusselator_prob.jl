brusselator_f(x, y, t) = ifelse((((x-0.3)^2 + (y-0.6)^2) <= 0.1^2) &&
                                (t >= 1.1), 5., 0.)
function limit(a, N)
  if a == N+1
    return 1
  elseif a == 0
    return N
  else
    return a
  end
end
function brusselator_2d_loop(du, u, p, t)
  @inbounds begin
    A, B, α, xyd, dx, N = p
    α = α/dx^2
    for I in CartesianIndices((N, N))
      x = xyd[I[1]]
      y = xyd[I[2]]
      i = I[1]
      j = I[2]
      ip1 = limit(i+1, N)
      im1 = limit(i-1, N)
      jp1 = limit(j+1, N)
      jm1 = limit(j-1, N)
      du[i,j,1] = α*(u[im1,j,1] + u[ip1,j,1] + u[i,jp1,1] + u[i,jm1,1] - 4u[i,j,1]) +
      B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1] + brusselator_f(x, y, t)
    end
    for I in CartesianIndices((N, N))
      i = I[1]
      j = I[2]
      ip1 = limit(i+1, N)
      im1 = limit(i-1, N)
      jp1 = limit(j+1, N)
      jm1 = limit(j-1, N)
      du[i,j,2] = α*(u[im1,j,2] + u[ip1,j,2] + u[i,jp1,2] + u[i,jm1,2] - 4u[i,j,2]) +
      A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
    end
  end
end
function init_brusselator_2d(xyd)
  N = length(xyd)
  u = zeros(N, N, 2)
  for I in CartesianIndices((N, N))
    x = xyd[I[1]]
    y = xyd[I[2]]
    u[I,1] = 22*(y*(1-y))^(3/2)
    u[I,2] = 27*(x*(1-x))^(3/2)
  end
  u
end
xyd_brusselator = range(0,stop=1,length=32)

@doc doc"""
2D Brusselator

```math
\\begin{align}
\\frac{\\partial u}{\\partial t} &= 1 + u^2v - 4.4u + \\alpha(\frac{\\partial^2 u}{\\partial x^2} + \frac{\\partial^2 u}{\\partial y^2}) + f(x, y, t)
\\frac{\\partial v}{\\partial t} &= 3.4u - u^2v + \\alpha(\frac{\\partial^2 u}{\\partial x^2} + \frac{\\partial^2 u}{\\partial y^2})
\\end{align}
```

where

```math
f(x, y, t) = \\begin{cases}
5 & \\quad \\text{if } (x-0.3)^2+(y-0.6)^2 ≤ 0.1^2 \\text{ and } t ≥ 1.1 \\\\
0 & \\quad \\text{else}
\\end{cases}
```

and the initial conditions are

```math
\\begin{align}
u(x, y, 0) &= 22\\cdot y(1-y)^{3/2} \\\\
v(x, y, 0) &= 27\\cdot x(1-x)^{3/2}
\\end{align}
```

with the periodic boundary condition

```math
\\begin{align}
u(x+1,y,t) &= u(x,y,t) \\\\
u(x,y+1,t) &= u(x,y,t)
\\end{align}
```

From Hairer Norsett Wanner Solving Ordinary Differential Equations II - Stiff and Differential-Algebraic Problems Page 152
"""
prob_ode_brusselator_2d = ODEProblem(brusselator_2d_loop,
                                     init_brusselator_2d(xyd_brusselator),
                                     (0.,11.5),
                                     (3.4, 1., 10.,
                                      xyd_brusselator, step(xyd_brusselator),
                                      length(xyd_brusselator)))

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
    mul!(buffer, D_brusselator_u, u)
    Du = buffer
    @. du[:, 1] = A + u^2*v - (B+1)*u + α*Du

    mul!(buffer, D_brusselator_v, v)
    Dv = buffer
    @. du[:, 2] = B*u - u^2*v + α*Dv
    nothing
end
function init_brusselator_1d(N)
  u = zeros(N, 2)
  x = range(0, stop=1, length=N)
  for i in 1:N
    u[i, 1] = 1 + sin(2pi*x[i])
    u[i, 2] = 3.
  end
  u
end

@doc doc"""
1D Brusselator

```math
\\begin{align}
\\frac{\\partial u}{\\partial t} &= A + u^2v - (B+1)u + \\alpha\frac{\\partial^2 u}{\\partial x^2}
\\frac{\\partial v}{\\partial t} &= Bu - u^2v + \\alpha\frac{\\partial^2 u}{\\partial x^2}
\\end{align}
```

and the initial conditions are

```math
\\begin{align}
u(x,0) &= 1+\\sin(2π x) \\\\
v(x,0) &= 3
\\end{align}
```

with the boundary condition

```math
\\begin{align}
u(0,t) &= u(1,t) = 1 \\\\
v(0,t) &= v(1,t) = 3
\\end{align}
```

From Hairer Norsett Wanner Solving Ordinary Differential Equations II - Stiff and Differential-Algebraic Problems Page 6
"""
prob_ode_brusselator_1d = ODEProblem(brusselator_1d,
                                    init_brusselator_1d(N_brusselator_1d),
                                    (0.,10.),
                                    (1., 3., 1/50, zeros(N_brusselator_1d)))

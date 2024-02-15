module BVProblemLibrary

using DiffEqBase, Markdown, SpecialFunctions

include("linear.jl")
include("nonlinear.jl")

################### flat_moon ############################
function flat_moon_f!(du, u, p, t)
    g = 1.62
    A = 3 * 1.62
    du[1] = u[3]
    du[2] = u[4]
    du[3] = A * cos(u[5])
    du[4] = (A * sin(u[5]) - g)
    du[5] = -u[6] * cos(u[5])
    du[6] = (u[6])^2 * sin(u[5])
    du[7] = 0
end
function flat_moon_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
    res_a[2] = u_a[2]
    res_a[3] = u_a[3]
    res_a[4] = u_a[4]
end
function flat_moon_bcb!(res_b, u_b, p)
    h = 185.2
    Vc = 1627
    res_b[1] = u_b[5] - h
    res_b[2] = u_b[6] - Vc
    res_b[3] = u_b[7]
end
flat_moon_tspan = (0, 700)
flat_moon_function = BVPFunction(flat_moon_f!, (flat_moon_bca!, flat_moon_bcb!),
    bcresid_prototype = (zeros(4), zeros(3)), twopoint = Val(true))
@doc raw"""
    flat_moon

This test problem is about the optimal-time launching of a satellite into orbit from flat moon without atmospheric drag.

Given by

```math
\frac{dz_1}{dt}=z_3t_f
```
```math
\frac{dz_2}{dt}=z_4t_f
```
```math
\frac{dz_3}{dt}=A\cos(z_5)t_f
```
```math
\frac{dz_4}{dt}=(A\sin(z_5)-g)t_f
```
```math
\frac{dz_5}{dt}=-z_6\cos(z_5)t_F
```
```math
\frac{dz_6}{dt}=z_6^2\sin(z_5)t_f
```
```math
\frac{dz_7}{dt}=0
```

with boundary condition

```math
z_1(0)=0
```
```math
z_2(0)=0
```
```math
z_3(0)=0
```
```math
z_4(0)=0
```
```math
z_5(1)=h
```
```math
z_6(1)=V_c
```
```math
z_7(1)=0
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=534)
"""
flat_moon = BVProblem(flat_moon_function, [0, 0, 0, 0, 0, 0, 0], flat_moon_tspan)

################### flat_earth ############################
function flat_earth_f!(du, u, p, t)
    Vc = sqrt(398600.4 / (6378.14 + 300)) * 1000
    h = 300000
    g = 9.80665
    Thrust2Weight = 3
    acc = Thrust2Weight * g

    du[1] = (u[3] * (Vc / h))
    du[2] = (u[4] * (Vc / h))
    du[3] = (acc * (1 / (abs(Vc) * sqrt(1 + u[6]^2))))
    du[4] = (acc * (u[6] / (abs(Vc) * sqrt(1 + u[6]^2))) - (g / Vc))
    du[5] = 0
    du[6] = (-u[5] * (Vc / h))
    du[7] = 0
end
function flat_earth_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
    res_a[2] = u_a[2]
    res_a[3] = u_a[3]
    res_a[4] = u_a[4]
end
function flat_earth_bcb!(res_b, u_b, p)
    res_b[1] = u_b[2] - 1
    res_b[2] = u_b[3] - 1
    res_b[3] = u_b[4]
end
flat_earth_function = BVPFunction(flat_earth_f!, (flat_earth_bca!, flat_earth_bcb!),
    bcresid_prototype = (zeros(4), zeros(3)), twopoint = Val(true))
flat_earth_tspan = (0, 700)
@doc raw"""
    flat_earth

Launch of a satellite into circular orbit from a flat Earth where we assume a uniform gravitational field ``g``.

Given by

```math
\frac{dz_1}{dt}=z_3\frac{V_c}{h}
```
```math
\frac{dz_2}{dt}=z_4\frac{V_c}{h}
```
```math
\frac{dz_3}{dt}=acc\frac{1}{|V_c|\sqrt{1+z_6^2}}
```
```math
\frac{dz_4}{dt}=acc\frac{1}{|V_c|\sqrt{1+z_6^2}}-frac{g}{V_c}
```
```math
\frac{dz_5}{dt}=0
```
```math
\frac{dz_6}{dt}=-z_5\frac{V_c}{h}
```
```math
\frac{dz_7}{dt}=0
```

with boundary condition

```math
z_1(0)=0
```
```math
z_2(0)=0
```
```math
z_3(0)=0
```
```math
z_4(0)=0
```
```math
z_5(1)=h
```
```math
z_6(1)=V_c
```
```math
z_7(1)=0
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=538)
"""
flat_earth = BVProblem(flat_earth_function, [0, 0, 0, 0, 0, 0, 0], flat_earth_tspan)

################### flat_earth_drag ############################
function flat_earth_drag_f!(du, u, p, t)
    fr = 2100000
    h = 180000
    m = 60880
    g_accel = 9.80665
    vc = 1000 * sqrt((398600.4) / (6378.14 + (h / 1000.0)))
    beta = 180000 / 840
    eta = 1.225 * 0.5 * 7.069 / 2
    xbardot = u[3] * (vc / h)
    ybardot = u[4] * (vc / h)
    Vxbardot = (fr / vc * (-u[6] / sqrt(u[6]^2.0 + u[7]^2.0)) -
                eta * exp(-u[2] * beta) * u[3] * sqrt(u[3]^2.0 + u[4]^2.0) * vc) / m
    Vybardot = (fr / vc * (-u[7] / sqrt(u[6]^2.0 + u[7]^2.0)) -
                eta * exp(-u[2] * beta) * u[4] * sqrt(u[3]^2.0 + u[4]^2.0) * vc) / m -
               g_accel / vc
    if sqrt(u[3]^2.0 + u[4]^2.0) == 0.0
        lambda_2_bar = 0.0
        lambda_3_bar = 0.0
        lambda_4_bar = -u[5] * (vc / h)
    else
        lambda_2_bar = -(u[6] * u[3] + u[7] * u[4]) * eta * beta *
                       sqrt(u[3]^2.0 + u[4]^2.0) * exp(-u[2] * beta) * vc / m
        lambda_3_bar = eta * exp(-u[2] * beta) * vc *
                       (u[6] * (2 * u[3]^2.0 + u[4]^2.0) + u[7] * u[3] * u[4]) /
                       sqrt(u[3]^2.0 + u[4]^2.0) / m
        lambda_4_bar = eta * exp(-u[2] * beta) * vc *
                       (u[7] * (u[3]^2.0 + 2.0 * u[4]^2.0) + u[6] * u[3] * u[4]) /
                       sqrt(u[3]^2.0 + u[4]^2.0) / m
    end
    du[1] = xbardot
    du[2] = ybardot
    du[3] = Vxbardot
    du[4] = Vybardot
    du[5] = lambda_2_bar
    du[6] = lambda_3_bar
    du[7] = lambda_4_bar
    du[8] = 0
end
function flat_earth_drag_bca!(res_a, u_a, p)
    res_a[1] = u_a[1]
    res_a[2] = u_a[2]
    res_a[3] = u_a[3]
    res_a[4] = u_a[4]
end
function flat_earth_drag_bcb!(res_b, u_b, p)
    fr = 2100000
    h = 180000
    m = 60880
    g_accel = 9.80665
    vc = 1000 * sqrt((398600.4) / (6378.14 + (h / 1000.0)))
    beta = 180000 / 840
    eta = 1.225 * 0.5 * 7.069 / 2
    res_b[1] = u_b[2] - 1
    res_b[2] = u_b[3] - 1
    res_b[3] = u_b[4]
    res_b[4] = (-sqrt(u_b[6]^2.0 + u_b[7]^2.0) * fr / m / vc -
                (u_b[6] * u_b[3]) * eta * exp(-beta) * sqrt(u_b[3]^2.0) * vc / m -
                u_b[7] * g_accel / vc) * u_b[8] + 1.0
end
flat_earth_drag_function = BVPFunction(
    flat_earth_drag_f!, (flat_earth_drag_bca!, flat_earth_drag_bcb!),
    bcresid_prototype = (zeros(4), zeros(4)), twopoint = Val(true))
flat_earth_drag_tspan = (0, 100)
@doc raw"""
    flat_earth_drag

Launch into circular orbit from a flat Earth including athmosferic drag.

Given by

```math
\frac{dz_1}{dt}=z_3\frac{V_c}{h}
```
```math
\frac{dz_2}{dt}=z_4\frac{V_c}{h}
```
```math
\frac{dz_3}{dt}=\frac{f}{V_c}(-\frac{z_6}{z_6^2+z_7^2}-V_c\eta\exp(-z_2\beta)z_3\sqrt{z_3^3+z_4^2})/m
```
```math
\frac{dz_4}{dt}=\frac{f}{V_c}(-\frac{z_7}{z_6^2+z_7^2}-V_c\eta\exp(-z_2\beta)z_4\sqrt{z_3^3+z_4^2})/m - g_{accel}/V_c
```
```math
\frac{dz_5}{dt}=-\eta\beta\exp(-z_2\beta)(z_6z_3+z_7z_4)\sqrt{z_3^3+z_4^2}\frac{V_c}{m}
```
```math
\frac{dz_6}{dt}=\eta\exp(-z_2\beta)(z_6(2z_3^2+z_4^2)+z_7z_3z_4)V_c/\sqrt{z_3^2+z_4^2}/m
```
```math
\frac{dz_7}{dt}=\eta\exp(-z_2\beta)(z_7(z_3^2+2z_4^2)+z_6z_3z_4)V_c/\sqrt{z_3^2+z_4^2}/m
```

with boundary condition

```math
z_1(0)=0
```
```math
z_2(0)=0
```
```math
z_3(0)=0
```
```math
z_4(0)=0
```
```math
z_5(1)=h
```
```math
z_6(1)=V_c
```
```math
z_7(1)=0
```

# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=544)
"""
flat_earth_drag = BVProblem(flat_earth_drag_function,
    [0, 0, 0, 0, 0, 0, 0, 0],
    flat_earth_drag_tspan)

################### measles ############################
function measles_f!(du, u, p, t)
    mu = 0.02
    lambda = 0.0279
    eta = 0.01
    B0 = 1575
    bt = B0 * (1.0 + (cos(2 * π * t)))

    du[1] = mu - (bt * u[1] * u[3])
    du[2] = (bt * u[1] * u[3]) - (u[2] / lambda)
    du[3] = (u[2] / lambda) - (u[3] / eta)
end

function measles_bc!(res, u, p, t)
    res[1] = u[1][1] - u[end][1]
    res[2] = u[1][2] - u[end][2]
    res[3] = u[1][3] - u[end][3]
end
measles_function = BVPFunction(measles_f!, measles_bc!)
measles_tspan = (0, 1)
@doc raw"""
    measles

This is an epidemiology model, about the spread of diseases.

Given by

```math
\frac{dy_1}{dt}=\mu-\beta(t)y_1y_3
```
```math
\frac{dy_2}{dt}=\beta(t)y_1y_3-y_2/\lambda
```
```math
\frac{dy_3}{dt}=y_2/\lambda-y_3/\eta
```

with boundary condition

```math
y(0)=y(1)
```


# Solution

No analytical solution

# References

[Reference](https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=555)
"""
measles = BVProblem(measles_function, [0, 0, 0], measles_tspan)

# Linear BVP Example Problems
export prob_bvp_linear_1, prob_bvp_linear_2, prob_bvp_linear_3, prob_bvp_linear_4,
       prob_bvp_linear_5,
       prob_bvp_linear_6, prob_bvp_linear_7, prob_bvp_linear_8, prob_bvp_linear_9,
       prob_bvp_linear_10,
       prob_bvp_linear_11, prob_bvp_linear_12, prob_bvp_linear_13, prob_bvp_linear_14,
       prob_bvp_linear_15,
       prob_bvp_linear_16, prob_bvp_linear_17, prob_bvp_linear_18

# Nonlinear BVP Example Problems
export prob_bvp_nonlinear_1, prob_bvp_nonlinear_2, prob_bvp_nonlinear_3,
       prob_bvp_nonlinear_4, prob_bvp_nonlinear_5,
       prob_bvp_nonlinear_6, prob_bvp_nonlinear_7, prob_bvp_nonlinear_8,
       prob_bvp_nonlinear_9,
       prob_bvp_nonlinear_10,
       prob_bvp_nonlinear_11, prob_bvp_nonlinear_12, prob_bvp_nonlinear_13,
       prob_bvp_nonlinear_14, prob_bvp_nonlinear_15

export flat_moon, flat_earth, flat_earth_drag, measles

end # module BVProblemLibrary

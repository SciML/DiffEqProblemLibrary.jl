__precompile__(false)

module SDEProblemLibrary

using DiffEqBase, Catalyst, Markdown

import RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

#SDE Example Problems
export prob_sde_wave, prob_sde_linear, prob_sde_cubic, prob_sde_2Dlinear,
       prob_sde_lorenz, prob_sde_2Dlinear, prob_sde_additive,
       prob_sde_additivesystem, oval2ModelExample, prob_sde_bistable,
       prob_sde_bruss, prob_sde_oscilreact
#SDE Stratonovich Example Problems
export prob_sde_linear_stratonovich, prob_sde_2Dlinear_stratonovich

### SDE Examples

f_linear(u, p, t) = 1.01u
σ_linear(u, p, t) = 0.87u
linear_analytic(u0, p, t, W) = @.(u0*exp(0.63155t + 0.87W))

@doc doc"""
```math
du_t = αudt + βudW_t
```
where ``α=1.01``, ``β=0.87``, and initial condition ``u_0=1/2``, with solution

```math
u(u_0,p,t,W_t)=u_0\exp((α-\frac{β^2}{2})t+βW_t)
```

"""
prob_sde_linear = SDEProblem(
    SDEFunction(f_linear, σ_linear,
        analytic = linear_analytic), 1 / 2,
    (0.0, 1.0))

linear_analytic_strat(u0, p, t, W) = @.(u0*exp(1.01t + 0.87W))
prob_sde_linear_stratonovich = SDEProblem(
    SDEFunction(f_linear, σ_linear,
        analytic = linear_analytic_strat), 1 / 2, (0.0, 1.0))
f_linear_iip(du, u, p, t) = @.(du=1.01 * u)
σ_linear_iip(du, u, p, t) = @.(du=0.87 * u)
@doc doc"""
8 linear SDEs (as a 4x2 matrix):

```math
du_t = αudt + βudW_t
```
where ``α=1.01``, ``β=0.87``, and initial condition ``u_0=\frac{1}{2}`` with solution

```math
u(u_0,p,t,W_t)=u_0\exp((α-\frac{β^2}{2})t+βW_t)
```
"""
prob_sde_2Dlinear = SDEProblem(
    SDEFunction(f_linear_iip, σ_linear_iip,
        analytic = linear_analytic),
    ones(4, 2) / 2, (0.0, 1.0))
prob_sde_2Dlinear_stratonovich = SDEProblem(
    SDEFunction(f_linear_iip, σ_linear_iip,
        analytic = linear_analytic_strat),
    ones(4, 2) / 2, (0.0, 1.0))

f_cubic(u, p, t) = -0.25 * u * (1 - u^2)
σ_cubic(u, p, t) = 0.5 * (1 - u^2)
cubic_analytic(u0, p, t, W) = @. ((1 + u0) * exp(W) + u0 - 1) / ((1 + u0) * exp(W) + 1 - u0)
ff_cubic = SDEFunction(f_cubic, σ_cubic, analytic = cubic_analytic)
@doc doc"""
```math
du_t = \frac{1}{4}u(1-u^2)dt + \frac{1}{2}(1-u^2)dW_t
```

and initial condition ``u_0=\frac{1}{2}``, with solution

```math
u(u0,p,t,W_t)=\frac{(1+u_0)\exp(W_t)+u)0-1}{(1+u_0)\exp(W_t)+1-u_0}
```
"""
prob_sde_cubic = SDEProblem(ff_cubic, 1 / 2, (0.0, 1.0))

f_wave(u, p, t) = @. -0.01 * sin(u) * cos(u)^3
σ_wave(u, p, t) = @. 0.1 * cos(u)^2
wave_analytic(u0, p, t, W) = @. atan(0.1 * W + tan(u0))
ff_wave = SDEFunction(f_wave, σ_wave, analytic = wave_analytic)
@doc doc"""
```math
du_t = -\frac{1}{100}\sin(u)\cos^3(u)dt + \frac{1}{10}\cos^{2}(u_t) dW_t
```

and initial condition ``u_0=1`` with solution

```math
u(u_0,p,t,W_t)=\arctan(\frac{W_t}{10} + \tan(u_0))
```
"""
prob_sde_wave = SDEProblem(ff_wave, 1.0, (0.0, 1.0))

f_additive(u, p, t) = @. p[2] / sqrt(1 + t) - u / (2 * (1 + t))
σ_additive(u, p, t) = @. p[1] * p[2] / sqrt(1 + t)
p = (0.1, 0.05)
additive_analytic(u0, p, t, W) = @. u0 / sqrt(1 + t) + p[2] * (t + p[1] * W) / sqrt(1 + t)
ff_additive = SDEFunction(f_additive, σ_additive, analytic = additive_analytic)
@doc doc"""
Additive noise problem

```math
u_t = (\frac{β}{\sqrt{1+t}}-\frac{1}{2(1+t)}u_t)dt + \frac{αβ}{\sqrt{1+t}}dW_t
```

and initial condition ``u_0=1`` with ``α=0.1`` and ``β=0.05``, with solution

```math
u(u_0,p,t,W_t)=\frac{u_0}{\sqrt{1+t}} + \frac{β(t+αW_t)}{\sqrt{1+t}}
```
"""
prob_sde_additive = SDEProblem(ff_additive, 1.0, (0.0, 1.0), p)

f_additive_iip(du, u, p, t) = @.(du=p[2] / sqrt(1 + t) - u / (2 * (1 + t)))
σ_additive_iip(du, u, p, t) = @.(du=p[1] * p[2] / sqrt(1 + t))
ff_additive_iip = SDEFunction(f_additive_iip, σ_additive_iip, analytic = additive_analytic)
p = ([0.1; 0.1; 0.1; 0.1], [0.5; 0.25; 0.125; 0.1115])
@doc doc"""
A multiple dimension extension of `additiveSDEExample`

"""
prob_sde_additivesystem = SDEProblem(ff_additive_iip, [1.0; 1.0; 1.0; 1.0],
    (0.0, 1.0), p)

function f_lorenz(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
end
σ_lorenz(du, u, p, t) = @.(du=3.0)
@doc doc"""
Lorenz Attractor with additive noise

```math
dx = σ(y-x)dt + αdW_t
```
```math
dy = (x(ρ-z) - y)dt + αdW_t
```
```math
dz = (xy - βz)dt + αdW_t
```

with ``σ=10``, ``ρ=28``, ``β=8/3``, ``α=3.0`` and initial condition ``u_0=[1;1;1]``.
"""
prob_sde_lorenz = SDEProblem(f_lorenz, σ_lorenz, ones(3), (0.0, 10.0), (10.0, 28.0, 2.66))

f_nltest(u, p, t) = (1 / 3) * u^(1 / 3) + 6 * u^(2 / 3)
σ_nltest(u, p, t) = u^(2 / 3)
analytic_nltest(u0, p, t, W) = (2t + 1 + W / 3)^3
ff_nltest = SDEFunction(f_nltest, σ_nltest, analytic = analytic_nltest)
@doc doc"""
Runge–Kutta methods for numerical solution of stochastic differential equations
Tocino and Ardanuy
"""
prob_sde_nltest = SDEProblem(ff_nltest, 1.0, (0.0, 10.0))

@doc doc"""
oval2ModelExample(;largeFluctuations=false,useBigs=false,noiseLevel=1)

A function which generates the Oval2 Epithelial-Mesenchymal Transition model
from:

Rackauckas, C., & Nie, Q. (2017). Adaptive methods for stochastic differential equations
via natural embeddings and rejection sampling with memory. Discrete and continuous
dynamical systems. Series B, 22(7), 2731.

19 SDEs which are only stiff during transitions between biological states.
"""
function oval2ModelExample(; largeFluctuations = false, useBigs = false, noiseLevel = 1)
    #Parameters
    J1_200 = 3.0
    J1_34 = 0.15
    J2_200 = 0.2
    J2_34 = 0.35
    J_2z = 0.9
    J_O = 0.918
    J_SO = 0.5
    # J_ZEB=0.06
    J_ecad1 = 0.1
    J_ecad2 = 0.3
    J_ncad1 = 0.4
    J_ncad2 = 0.4
    J_ncad3 = 2
    J_snail0 = 0.6
    J_snail1 = 1.8
    J_zeb = 3.0
    K1 = 1.0
    K2 = 1.0
    K3 = 1.0
    K4 = 1.0
    K5 = 1.0
    KTGF = 20.0
    Ks = 100.0
    # TGF0=0
    TGF_flg = 0.0
    Timescale = 1000.0
    dk_ZR1 = 0.5
    dk_ZR2 = 0.5
    dk_ZR3 = 0.5
    dk_ZR4 = 0.5
    dk_ZR5 = 0.5
    k0O = 0.35
    k0_200 = 0.0002
    k0_34 = 0.001
    k0_snail = 0.0005
    k0_zeb = 0.003
    kO = 1.2
    kOp = 10.0
    k_200 = 0.02
    k_34 = 0.019
    k_OT = 1.1
    k_SNAIL = 16.0
    k_TGF = 1.5
    k_ZEB = 16.0
    k_ecad0 = 5.0
    k_ecad1 = 15.0
    k_ecad2 = 5.0
    k_ncad0 = 5.0
    k_ncad1 = 2.0
    k_ncad2 = 5.0
    k_snail = 0.05
    k_tgf = 0.05
    k_zeb = 0.06
    kdO = 1.0
    kd_200 = 0.035
    kd_34 = 0.035
    kd_SNAIL = 1.6
    kd_SR1 = 0.9
    kd_TGF = 0.9
    kd_ZEB = 1.66
    kd_ecad = 0.05
    kd_ncad = 0.05
    kd_snail = 0.09
    kd_tgf = 0.1
    kd_tgfR = 1.0
    kd_zeb = 0.1
    kd_Op = 10.0
    lamda1 = 0.5
    lamda2 = 0.5
    lamda3 = 0.5
    lamda4 = 0.5
    lamda5 = 0.5
    lamdas = 0.5
    lamdatgfR = 0.8
    nO = 6.0
    nSO = 2.0
    nzo = 2.0
    GE = 1.0
    f = function (dy, y, p, t)
        # y(1) = snailt
        # y(2) = SNAIL
        # y(3) = miR34t
        # y(4) = SR1 # abundance of SNAIL/miR34 complex
        # y(5) = zebt
        # y(6) = ZEB
        # y(7) = miR200t
        # y(8) = ZR1 # abundance of ZEB/miR200 complex with i copies of miR200 bound on the sequence of ZEB1
        # y(9) = ZR2
        # y(10) = ZR3
        # y(11) = ZR4
        # y(12) = ZR5
        # y(13) = tgft
        # y(14) = TGF
        # y(15) = tgfR # abundance of TGF/miR200 complex
        # y(16) = Ecad
        # y(17) = Ncad
        # y(18) = Ovol2
        TGF0 = 0.5(t > 100)
        #ODEs
        dy[1] = k0_snail +
                k_snail * (((y[14] + TGF0) / J_snail0))^2 /
                (1 + (((y[14] + TGF0) / J_snail0))^2 + (y[19] / J_SO)^nSO) /
                (1 + y[2] / J_snail1) - kd_snail * (y[1] - y[4]) - kd_SR1 * y[4]
        dy[2] = k_SNAIL * (y[1] - y[4]) - kd_SNAIL * y[2]
        dy[3] = k0_34 + k_34 / (1 + ((y[2] / J1_34))^2 + ((y[6] / J2_34))^2) -
                kd_34 * (y[3] - y[4]) - kd_SR1 * y[4] + lamdas * kd_SR1 * y[4]
        dy[4] = Timescale * (Ks * (y[1] - y[4]) * (y[3] - y[4]) - y[4])
        dy[5] = k0_zeb +
                k_zeb * ((y[2] / J_zeb))^2 /
                (1 + ((y[2] / J_zeb))^2 + ((y[19] / J_2z))^nO) -
                kd_zeb * (y[5] - (5 * y[8] + 10 * y[9] + 10 * y[10] + 5 * y[11] + y[12])) -
                dk_ZR1 * 5 * y[8] - dk_ZR2 * 10 * y[9] - dk_ZR3 * 10 * y[10] -
                dk_ZR4 * 5 * y[11] - dk_ZR5 * y[12]
        dy[6] = k_ZEB * (y[5] - (5 * y[8] + 10 * y[9] + 10 * y[10] + 5 * y[11] + y[12])) -
                kd_ZEB * y[6]
        dy[7] = k0_200 + k_200 / (1 + ((y[2] / J1_200))^3 + ((y[6] / J2_200))^2) -
                kd_200 * (y[7] -
                 (5 * y[8] + 2 * 10 * y[9] + 3 * 10 * y[10] + 4 * 5 * y[11] + 5 * y[12]) -
                 y[15]) - dk_ZR1 * 5 * y[8] - dk_ZR2 * 2 * 10 * y[9] -
                dk_ZR3 * 3 * 10 * y[10] - dk_ZR4 * 4 * 5 * y[11] - dk_ZR5 * 5 * y[12] +
                lamda1 * dk_ZR1 * 5 * y[8] + lamda2 * dk_ZR2 * 2 * 10 * y[9] +
                lamda3 * dk_ZR3 * 3 * 10 * y[10] + lamda4 * dk_ZR4 * 4 * 5 * y[11] +
                lamda5 * dk_ZR5 * 5 * y[12] - kd_tgfR * y[15] + lamdatgfR * kd_tgfR * y[15]
        dy[8] = Timescale * (K1 *
                 (y[7] -
                  (5 * y[8] + 2 * 10 * y[9] + 3 * 10 * y[10] + 4 * 5 * y[11] + 5 * y[12]) -
                  y[15]) *
                 (y[5] - (5 * y[8] + 10 * y[9] + 10 * y[10] + 5 * y[11] + y[12])) - y[8])
        dy[9] = Timescale * (K2 *
                 (y[7] -
                  (5 * y[8] + 2 * 10 * y[9] + 3 * 10 * y[10] + 4 * 5 * y[11] + 5 * y[12]) -
                  y[15]) * y[8] - y[9])
        dy[10] = Timescale * (K3 *
                  (y[7] -
                   (5 * y[8] + 2 * 10 * y[9] + 3 * 10 * y[10] + 4 * 5 * y[11] + 5 * y[12]) -
                   y[15]) * y[9] - y[10])
        dy[11] = Timescale * (K4 *
                  (y[7] -
                   (5 * y[8] + 2 * 10 * y[9] + 3 * 10 * y[10] + 4 * 5 * y[11] + 5 * y[12]) -
                   y[15]) * y[10] - y[11])
        dy[12] = Timescale * (K5 *
                  (y[7] -
                   (5 * y[8] + 2 * 10 * y[9] + 3 * 10 * y[10] + 4 * 5 * y[11] + 5 * y[12]) -
                   y[15]) * y[11] - y[12])
        dy[13] = k_tgf - kd_tgf * (y[13] - y[15]) - kd_tgfR * y[15]
        dy[14] = k_OT + k_TGF * (y[13] - y[15]) - kd_TGF * y[14]
        dy[15] = Timescale * (TGF_flg +
                  KTGF *
                  (y[7] -
                   (5 * y[8] + 2 * 10 * y[9] + 3 * 10 * y[10] + 4 * 5 * y[11] + 5 * y[12]) -
                   y[15]) * (y[13] - y[15]) - y[15])
        dy[16] = GE * (k_ecad0 + k_ecad1 / (((y[2] / J_ecad1))^2 + 1) +
                  k_ecad2 / (((y[6] / J_ecad2))^2 + 1) - kd_ecad * y[16])
        dy[17] = k_ncad0 + k_ncad1 * (((y[2] / J_ncad1))^2) / (((y[2] / J_ncad1))^2 + 1) +
                 k_ncad2 * (((y[6] / J_ncad2))^2) /
                 (((y[6] / J_ncad2)^2 + 1) * (1 + y[19] / J_ncad3)) - kd_ncad * y[17]
        dy[18] = k0O + kO / (1 + ((y[6] / J_O))^nzo) - kdO * y[18]
        dy[19] = kOp * y[18] - kd_Op * y[19]
    end

    σ1 = function (dσ, y, p, t)
        dσ[1] = noiseLevel * 1.5y[1]
        dσ[18] = noiseLevel * 6y[18]
    end

    σ2 = function (dσ, y, p, t)
        dσ[1] = 0.02y[1]
        dσ[16] = 0.02y[16]
        dσ[18] = 0.2y[18]
        dσ[17] = 0.02y[17]
    end

    if largeFluctuations
        σ = σ1
    else
        σ = σ2
    end

    u0 = [0.128483, 1.256853, 0.0030203, 0.0027977, 0.0101511, 0.0422942, 0.2391346,
        0.0008014, 0.0001464, 2.67e-05, 4.8e-6, 9e-7, 0.0619917, 1.2444292, 0.0486676,
        199.9383546, 137.4267984, 1.5180203, 1.5180203] #Fig 9B
    if useBigs
        u0 = big.(u0)
    end
    #u0 =  [0.1701;1.6758;0.0027;0.0025;0.0141;0.0811;0.1642;0.0009;0.0001;0.0000;0.0000;0.0000;0.0697;1.2586;0.0478;194.2496;140.0758;1.5407;1.5407] #Fig 9A
    SDEProblem(f, σ, u0, (0.0, 500.0))
end

stiff_quad_f_ito(u, p, t) = -(p[1] + (p[2]^2) * u) * (1 - u^2)
stiff_quad_f_strat(u, p, t) = -p[1] * (1 - u^2)
stiff_quad_g(u, p, t) = p[2] * (1 - u^2)

function stiff_quad_f_ito_analytic(u0, p, t, W)
    α = p[1]
    β = p[2]
    exp_tmp = exp(-2 * α * t + 2 * β * W)
    tmp = 1 + u0
    (tmp * exp_tmp + u0 - 1) / (tmp * exp_tmp - u0 + 1)
end

ff_stiff_quad_ito = SDEFunction(stiff_quad_f_ito, stiff_quad_g,
    analytic = stiff_quad_f_ito_analytic)

function stiff_quad_f_strat_analytic(u0, p, t, W)
    α = p[1]
    β = p[2]
    exp_tmp = exp(-2 * α * t + 2 * β * W)
    tmp = 1 + u0
    (tmp * exp_tmp + u0 - 1) / (tmp * exp_tmp - u0 + 1)
end

ff_stiff_quad_strat = SDEFunction(stiff_quad_f_strat, stiff_quad_g,
    analytic = stiff_quad_f_strat_analytic)

@doc doc"""
The composite Euler method for stiff stochastic
differential equations

Kevin Burrage, Tianhai Tian

And

S-ROCK: CHEBYSHEV METHODS FOR STIFF STOCHASTIC
DIFFERENTIAL EQUATIONS

ASSYR ABDULLE AND STEPHANE CIRILLI

Stiffness of Euler is determined by α+β²<1
Higher α or β is stiff, with α being deterministic stiffness and
β being noise stiffness (and grows by square).
"""
prob_sde_stiffquadito = SDEProblem(ff_stiff_quad_ito, 0.5, (0.0, 3.0),
    (1.0, 1.0))

@doc doc"""
The composite Euler method for stiff stochastic
differential equations

Kevin Burrage, Tianhai Tian

And

S-ROCK: CHEBYSHEV METHODS FOR STIFF STOCHASTIC
DIFFERENTIAL EQUATIONS

ASSYR ABDULLE AND STEPHANE CIRILLI

Stiffness of Euler is determined by α+β²<1
Higher α or β is stiff, with α being deterministic stiffness and
β being noise stiffness (and grows by square).
"""
prob_sde_stiffquadstrat = SDEProblem(ff_stiff_quad_strat, 0.5, (0.0, 3.0),
    (1.0, 1.0))

@doc doc"""
Stochastic Heat Equation with scalar multiplicative noise

S-ROCK: CHEBYSHEV METHODS FOR STIFF STOCHASTIC
DIFFERENTIAL EQUATIONS

ASSYR ABDULLE AND STEPHANE CIRILLI

Raising D or k increases stiffness
"""
function generate_stiff_stoch_heat(D = 1, k = 1; N = 100)
    A = full(Tridiagonal([1.0 for i in 1:(N - 1)], [-2.0 for i in 1:N],
        [1.0 for i in 1:(N - 1)]))
    dx = 1 / N
    A = D / (dx^2) * A
    function f(du, u, p, t)
        mul!(du, A, u)
    end
    #=
    function f(::Type{Val{:analytic}},u0,p,t,W)
        expm(A*t+W*I)*u0
    end
    =#
    function g(du, u, p, t)
        @. du = k * u
    end
    SDEProblem(f, g, ones(N), (0.0, 3.0),
        noise = WienerProcess(0.0, 0.0, 0.0, rswm = RSWM(adaptivealg = :RSwM3)))
end

bistable_f(du, u, p, t) = du[1] = p[1] + p[2] * u[1]^4 / (u[1]^4 + 11.9^4) - p[3] * u[1]
function bistable_g(du, u, p, t)
    du[1, 1] = 0.1 * sqrt(p[1] + p[2] * u[1]^4 / (u[1]^4 + 11.9^4))
    du[1, 2] = -0.1 * sqrt(p[3] * u[1])
end
p = (5.0, 18.0, 1.0)
"""
Bistable chemical reaction network with a semi-stable lower state.
"""
prob_sde_bistable = SDEProblem(bistable_f, bistable_g, [3.0], (0.0, 300.0), p,
    noise_rate_prototype = zeros(1, 2))

function bruss_f(du, u, p, t)
    du[1] = p[1] + p[2] * u[1] * u[1] * u[2] - p[3] * u[1] - p[4] * u[1]
    du[2] = p[3] * u[1] - p[2] * u[1] * u[1] * u[2]
end
function bruss_g(du, u, p, t)
    du[1, 1] = 0.15 * sqrt(p[1])
    du[1, 2] = 0.15 * sqrt(p[2] * u[1] * u[1] * u[2])
    du[1, 3] = -0.15 * sqrt(p[3] * u[1])
    du[1, 4] = -0.15 * sqrt(p[4] * u[1])
    du[2, 1] = 0
    du[2, 2] = -0.15 * sqrt(p[2] * u[1] * u[1] * u[2])
    du[2, 3] = 0.15 * sqrt(p[3] * 2.5 * u[1])
    du[2, 4] = 0
end
p = (1.0, 1.0, 2.5, 1.0)
"""
Stochastic Brusselator
"""
prob_sde_bruss = SDEProblem(bruss_f, bruss_g, [3.0, 2.0], (0.0, 100.0), p,
    noise_rate_prototype = zeros(2, 4))

network = @reaction_network begin
    @parameters p1=0.01 p2=3.0 p3=3.0 p4=4.5 p5=2.0 p6=15.0 p7=20.0 p8=0.005 p9=0.01 p10=0.05
    p1, (X, Y, Z) --> 0
    hill(X, p2, 100.0, -4), 0 --> Y
    hill(Y, p3, 100.0, -4), 0 --> Z
    hill(Z, p4, 100.0, -4), 0 --> X
    hill(X, p5, 100.0, 6), 0 --> R
    hill(Y, p6, 100.0, 4) * 0.002, R --> 0
    p7, 0 --> S
    R * p8, S --> SP
    p9, SP + SP --> SP2
    p10, SP2 --> 0
end

"""
An oscillatory chemical reaction system
"""
prob_sde_oscilreact = SDEProblem(network, [200.0, 60.0, 120.0, 100.0, 50.0, 50.0, 50.0],
    (0.0, 4000.0), eval_module = @__MODULE__)

end # module

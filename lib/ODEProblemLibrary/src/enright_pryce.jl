using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit.Symbolics
using Symbolics: unwrap
using DiffEqBase, StaticArrays, LinearAlgebra

@variables y(t)[1:10]
y = collect(y)

# Stiff
sa1sys = let
    sa1eqs = [D(y[1]) ~ -0.5 * y[1], D(y[2]) ~ -y[2], D(y[3]) ~ -100*y[3], D(y[4]) ~ -90*y[4]]
    ODESystem(sa1eqs, t, name = :sa1)
end

sa1prob = ODEProblem{false}(structural_simplify(sa1sys), y[1:4] .=> 1.0, (0, 20.0), dt = 1e-2)

sa2sys = let
    sa2eqs = [D(y[1]) ~ -1800 * y[1] + 900 * y[2]]
    for i in 2:8
        push!(sa2eqs, D(y[i]) ~ y[i - 1] - 2 * y[i] + y[i + 1])
    end
    push!(sa2eqs, D(y[9]) ~ 1000 * y[8] - 2000 * y[9] + 1000)
    ODESystem(sa2eqs, t, name = :sa2)
end

sa2prob = ODEProblem{false}(structural_simplify(sa2sys), y[1:9] .=> 0.0, (0, 120.0), dt = 5e-4)

sa3sys = let
    sa3eqs = [
              D(y[1]) ~ -1e4 * y[1] + 100 * y[2] - 10 * y[3] + y[4],
              D(y[2]) ~ -1e3 * y[2] + 10 * y[3] - 10 * y[4],
              D(y[3]) ~ -y[3] + 10 * y[4],
              D(y[4]) ~ -0.1 * y[4],
             ]
    ODESystem(sa3eqs, t, name = :sa3)
end

sa3prob = ODEProblem{false}(structural_simplify(sa3sys), y[1:4] .=> 1.0, (0, 20.0), dt = 1e-5)

sa4sys = let
    sa4eqs = [D(y[i]) ~ -i^5 * y[i] for i in 1:10]
    ODESystem(sa4eqs, t, name = :sa4)
end

sa4prob = ODEProblem{false}(structural_simplify(sa4sys), y[1:10] .=> 1.0, (0, 1.0), dt = 1e-5)

#const SA_PROBLEMS = [sa1prob, sa2prob, sa3prob, sa4prob]

sb1sys = let
    sb1eqs = [D(y[1]) ~ -y[1] + y[2],
              D(y[2]) ~ -100y[1] - y[2],
              D(y[3]) ~ -100y[3] + y[4],
              D(y[4]) ~ -10_000y[3] - 100y[4]]

    ODESystem(sb1eqs, t, name = :sb1)
end

sb1prob = ODEProblem{false}(structural_simplify(sb1sys), [y[[1, 3]] .=> 1.0; y[[2, 4]] .=> 0.0;], (0, 20.0), dt = 7e-3)


@parameters α
sb2sys = let
    sb2eqs = [D(y[1]) ~ -10y[1] + α * y[2],
              D(y[2]) ~ -α * y[1] - 10 * y[2],
              D(y[3]) ~ -4y[3],
              D(y[4]) ~ -y[4],
              D(y[5]) ~ -0.5y[5],
              D(y[6]) ~ -0.1y[6]]

    ODESystem(sb2eqs, t, name = :sb2)
end

sb2prob = ODEProblem{false}(structural_simplify(sb2sys), y .=> 1.0, (0, 20.0), [α => 3], dt = 1e-2)
sb3prob = ODEProblem{false}(structural_simplify(sb2sys), y .=> 1.0, (0, 20.0), [α => 8], dt = 1e-2)
sb4prob = ODEProblem{false}(structural_simplify(sb2sys), y .=> 1.0, (0, 20.0), [α => 25], dt = 1e-2)
sb5prob = ODEProblem{false}(structural_simplify(sb2sys), y .=> 1.0, (0, 20.0), [α => 100], dt = 1e-2)


sc1sys = let
    sc1eqs = [
              D(y[1]) ~ -y[1] + y[2]^2 + y[3]^2 + y[4]^2,
              D(y[2]) ~ -10y[2] + 10*(y[3]^2 + y[4]^2),
              D(y[3]) ~ -40y[3] + 40 * y[4]^2,
              D(y[4]) ~ -100y[4] + 2]

    ODESystem(sc1eqs, t, name = :sc1)
end

sc1prob = ODEProblem{false}(structural_simplify(sc1sys), y .=> 1.0, (0, 20.0), dt = 1e-2)


@parameters β
sc2sys = let
    sc2eqs = [D(y[1]) ~ -y[1] + 2,
              D(y[2]) ~ -10y[2] + β * y[1]^2,
              D(y[3]) ~ -40y[3] + 4β * (y[1]^2 + y[2]^2),
              D(y[4]) ~ -100y[4] + 10β * (y[1]^2 + y[2]^2 + y[3]^2)]

    ODESystem(sc2eqs, t, name = :sc2)
end

sc2prob = ODEProblem{false}(structural_simplify(sc2sys), y .=> 1.0, (0, 20.0), [β => 0.1], dt = 1e-2)
sc3prob = ODEProblem{false}(structural_simplify(sc2sys), y .=> 1.0, (0, 20.0), [β => 1.0], dt = 1e-2)
sc4prob = ODEProblem{false}(structural_simplify(sc2sys), y .=> 1.0, (0, 20.0), [β => 10.0], dt = 1e-2)
sc5prob = ODEProblem{false}(structural_simplify(sc2sys), y .=> 1.0, (0, 20.0), [β => 20.0], dt = 1e-2)

sd1sys = let
    sd1eqs = [D(y[1]) ~ 0.2 * (y[2] - y[1]),
              D(y[2]) ~ 10y[1] - (60 - 0.125y[3])*y[2] + 0.125y[3],
              D(y[3]) ~ 1]

    ODESystem(sd1eqs, t, name = :sd1)
end

sd1prob = ODEProblem{false}(structural_simplify(sd1sys), y .=> 0.0, (0, 400.0), [β => 0.1], dt = 1.7e-2)

sd2sys = let
    sd2eqs = [D(y[1]) ~ -0.04y[1] + 0.01 * (y[2] * y[3]),
              D(y[2]) ~ 400y[1] - 100 * (y[2] * y[3]) - 3000 * y[2]^2,
              D(y[3]) ~ 30y[2]^2]

    ODESystem(sd2eqs, t, name = :sd2)
end

sd2prob = ODEProblem{false}(structural_simplify(sd2sys), [y[1] => 1.0; y[2:3] .=> 0.0], (0, 40.0), dt = 1e-5)

sd3sys = let
    sd3eqs = [D(y[1]) ~ y[3] - 100 * (y[1] * y[2]),
              D(y[2]) ~ y[3] + 2y[4] - 100 * (y[1] * y[2]) - 2e4 * y[2]^2,
              D(y[3]) ~ -y[3] + 100 * (y[1] * y[2]),
              D(y[4]) ~ -y[4] + 1e4*y[2]^2]

    ODESystem(sd3eqs, t, name = :sd3)
end

sd3prob = ODEProblem{false}(structural_simplify(sd3sys), [y[1:2] .=> 1; y[3:4] .=> 0.0], (0, 20.0), dt = 2.5e-5)

sd4sys = let
    sd4eqs = [D(y[1]) ~ -0.013y[1] - 1000 * (y[1] * y[3]),
              D(y[2]) ~ -2500*(y[2] * y[3]),
              D(y[3]) ~ -0.013y[1] - 1000 * (y[1] * y[3]) - 2500 * (y[2] * y[3])]

    ODESystem(sd4eqs, t, name = :sd4)
end

sd4prob = ODEProblem{false}(structural_simplify(sd4sys), [y[1:2] .=> 1; y[3] => 0.0], (0, 50.0), dt = 2.9e-4)

sd5sys = let
    sd5eqs = [D(y[1]) ~ 0.01 - (1 + (y[1] + 1000) * (y[1] + 1)) * (0.01 + y[1] + y[2]),
              D(y[2]) ~ 0.01 - (1 + y[2]^2) * (0.01 + y[1] + y[2])]

    ODESystem(sd5eqs, t, name = :sd5)
end

sd5prob = ODEProblem{false}(structural_simplify(sd5sys), y[1:2] .=> 0.0, (0, 100.0), dt = 1e-4)

sd6sys = let
    sd6eqs = [D(y[1]) ~ -y[1] + 10^8 * y[3] * (1 - y[1]),
              D(y[2]) ~ -10y[2] + 3e7 * y[3] * (1 - y[2]),
              D(y[3]) ~ -(-y[1] + 10^8 * y[3] * (1 - y[1])) - (-10y[2] + 3e7 * y[3] * (1 - y[2])),
             ]

    ODESystem(sd6eqs, t, name = :sd6)
end

sd6prob = ODEProblem{false}(structural_simplify(sd6sys), [y[1] => 1.0; y[2:3] .=> 0.0], (0, 1.0), dt = 3.3e-8)

se1sys = let
    Γ = 100
    se1eqs = [D(y[1]) ~ y[2],
              D(y[2]) ~ y[3],
              D(y[3]) ~ y[4],
              D(y[4]) ~ (y[1]^2 - sin(y[1]) - Γ^4) * y[1] + (y[2] * y[3] / (y[1]^2 + 1) - 4 * Γ^3) * y[2] + (1 - 6 * Γ^2) * y[3] + (10 * exp(-y[4]^2) - 4Γ) * y[4] + 1
             ]

    ODESystem(se1eqs, t, name = :se1)
end

se1prob = ODEProblem{false}(structural_simplify(se1sys), y .=> 0.0, (0, 1.0), dt = 6.8e-3)

se2sys = let
    se2eqs = [D(y[1]) ~ y[2],
              D(y[2]) ~ 5 * (1 - y[1]^2) * y[2] - y[1]
             ]

    ODESystem(se2eqs, t, name = :se2)
end

se2prob = ODEProblem{false}(structural_simplify(se2sys), [y[1] => 2.0, y[2] => 0.0], (0, 1.0), dt = 1e-3)

se3sys = let
    se3eqs = [D(y[1]) ~ -(55 + y[3]) * y[1] + 65 * y[2],
              D(y[2]) ~ 0.0785 * (y[1] - y[2]),
              D(y[3]) ~ 0.1 * y[1]
             ]

    ODESystem(se3eqs, t, name = :se3)
end

se3prob = ODEProblem{false}(structural_simplify(se3sys), [y[1:2] .=> 1.0; y[3] => 0.0], (0, 500.0), dt = 0.02)

se4sys = let y = y
    y = y[1:4]
    U = ones(4, 4)
    U[diagind(U)] .= -1
    U ./= 2
    Z = U * y
    G = U' * [(Z[1]^2 - Z[2]^2) / 2, Z[1] * Z[2], Z[3]^2, Z[4]^2]
    A = [-10 -10 0 0; 10 -10 0 0; 0 0 1000 0; 0 0 0 0.01]
    se4eqs = D.(y) .~ -(U' * A * Z) + G

    ODESystem(se4eqs, t, name = :se4)
end

se4prob = ODEProblem{false}(structural_simplify(se4sys), [y[1] => 0.0; y[2] => -2.0; y[3:4] .=> -1.0], (0, 1000.0), dt = 1e-3)

se5sys = let
    se5eqs = [
              D(y[1]) ~ -7.89e-10 * y[1] - 1.1e7 * y[1] * y[3],
              D(y[2]) ~ 7.89e-10 * y[1] - 1.13e9 * y[2] * y[3],
              D(y[3]) ~ 7.89e-10 * y[1] - 1.1e7 * y[1] * y[3] + 1.13e3 * y[4] - 1.13e9 * y[2] * y[3],
              D(y[4]) ~ 1.1e7 * y[1] * y[3] - 1.13e3 * y[4],
             ]

    ODESystem(se5eqs, t, name = :se5)
end

se5prob = ODEProblem{false}(structural_simplify(se5sys), [y[1] => 1.76e-3; y[2:4] .=> 0.0], (0, 1000.0), dt = 1e-3)

sf1sys = let
    k = exp(20.7 - 1500 / y[1])
    sf1eqs = [
              D(y[1]) ~ 1.3 * (y[3] - y[1]) + 10400 * k * y[2],
              D(y[2]) ~ 1880 * [y[4] - y[2] * (1 + k)],
              D(y[3]) ~ 1752 - 269 * y[3] + 267 * y[1],
              D(y[4]) ~ 0.1 + 320 * y[2] - 321 * y[4],
             ]

    ODESystem(sf1eqs, t, name = :sf1)
end

sf1prob = ODEProblem{false}(structural_simplify(sf1sys), [y[1] => 761.0; y[2] => 0.0; y[3] => 600.0; y[4] => 0.1], (0, 1000.0), dt = 1e-4)

sf2sys = let
    sf2eqs = [
              D(y[1]) ~ -y[1] - y[1] * y[2] + 294 * y[2],
              D(y[2]) ~ y[1] * (1 - y[2]) / 98 - 3 * y[2]
             ]

    ODESystem(sf2eqs, t, name = :sf2)
end

sf2prob = ODEProblem{false}(structural_simplify(sf2sys), [y[1] => 1.0; y[2] => 0.0], (0, 240.0), dt = 1e-2)

sf3sys = let
    sf3eqs = [
              D(y[1]) ~ -1e7 * y[2] * y[1] + 10 * y[3],
              D(y[2]) ~ -1e7 * y[2] * y[1] - 1e7 * y[2] * y[5] + 10 * y[3] + 10 * y[4],
              D(y[3]) ~ 1e7 * y[2] * y[1] - 1.001e4 * y[3] + 1e-3 * y[4],
              D(y[4]) ~ 1e4 * y[3] - 10.001 * y[4] + 1e7 * y[2] * y[5],
              D(y[5]) ~ 10 * y[4] - 1e7 * y[2] * y[5],
             ]

    ODESystem(sf3eqs, t, name = :sf3)
end

sf3prob = ODEProblem{false}(structural_simplify(sf3sys), [y[1] => 4e-6; y[2] => 1e-6; y[3:5] .=> 0.0], (0, 100.0), dt = 1e-6)

sf4sys = let
    sf4eqs = [
              D(y[1]) ~ 77.27 * (y[2] - y[1] * y[2] + y[1] - 8.375e-6 * y[1]^2),
              D(y[2]) ~ (-y[2] - y[1] * y[2] + y[3]) / 77.27,
              D(y[3]) ~ 0.161 * (y[1] - y[3]),
             ]

    ODESystem(sf4eqs, t, name = :sf4)
end

sf4prob = ODEProblem{false}(structural_simplify(sf4sys), [y[1] => 4.0; y[2] => 1.1; y[3] => 4.0], (0, 300.0), dt = 1e-3)

sf5sys = let
    sf5eqs = [
              D(y[1]) ~ 1e11 * (-3 * y[1] * y[2] + 0.0012 * y[4] - 9 * y[1] * y[3]),
              D(y[2]) ~ -3e11 * y[1] * y[2] + 2e7 * y[4],
              D(y[3]) ~ 1e11 * (-9 * y[1] * y[3] + 0.001 * y[4]),
              D(y[4]) ~ 1e11 * (3 * y[1] * y[2] - 0.0012 * y[4] + 9 * y[1] * y[3]),
             ]

    ODESystem(sf5eqs, t, name = :sf5)
end

sf5prob = ODEProblem{false}(structural_simplify(sf5sys), [y[1] => 3.365e-7; y[2] => 8.261e-3; y[3] => 1.642e-3; y[4] => 9.38e-6], (0, 100.0), dt = 1e-7)

# Non-stiff
na1sys = let y = y[1]
    na1eqs = D(y) ~ -y

    ODESystem(na1eqs, t, name = :na1)
end

na1prob = ODEProblem{false}(structural_simplify(na1sys), [y[1] => 1], (0, 20.0))

na2sys = let y = y[1]
    na2eqs = D(y) ~ -y^3/2

    ODESystem(na2eqs, t, name = :na2)
end

na2prob = ODEProblem{false}(structural_simplify(na2sys), [y[1] => 1], (0, 20.0))

na3sys = let y = y[1]
    na3eqs = D(y) ~ y * cos(t)

    ODESystem(na3eqs, t, name = :na3)
end

na3prob = ODEProblem{false}(structural_simplify(na3sys), [y[1] => 1], (0, 20.0))

na4sys = let y = y[1]
    na4eqs = D(y) ~ y/4 * (1 - y/20)

    ODESystem(na4eqs, t, name = :na4)
end

na4prob = ODEProblem{false}(structural_simplify(na4sys), [y[1] => 1], (0, 20.0))

na5sys = let y = y[1]
    na5eqs = D(y) ~ (y - t) / (y + t)

    ODESystem(na5eqs, t, name = :na5)
end

na5prob = ODEProblem{false}(structural_simplify(na5sys), [y[1] => 4], (0, 20.0))

nb1sys = let
    nb1eqs = [
              D(y[1]) ~ 2 * (y[1] - y[1] * y[2]),
              D(y[2]) ~ -(y[2] - y[1] * y[2]),
             ]

    ODESystem(nb1eqs, t, name = :nb1)
end

nb1prob = ODEProblem{false}(structural_simplify(nb1sys), [y[1] => 1.0, y[2] => 3], (0, 20.0))

nb2sys = let
    nb2eqs = [
              D(y[1]) ~ -y[1] + y[2],
              D(y[2]) ~ y[1] - 2y[2] + y[3],
              D(y[3]) ~ y[2] - y[3],
             ]

    ODESystem(nb2eqs, t, name = :nb2)
end

nb2prob = ODEProblem{false}(structural_simplify(nb2sys), [y[1] => 2.0, y[2] => 0.0, y[3] => 1.0], (0, 20.0))

nb3sys = let
    nb3eqs = [
              D(y[1]) ~ -y[1],
              D(y[2]) ~ y[1] - y[2]^2,
              D(y[3]) ~ y[2]^2,
             ]

    ODESystem(nb3eqs, t, name = :nb3)
end

nb3prob = ODEProblem{false}(structural_simplify(nb3sys), [y[1] => 1.0; y[2:3] .=> 0.0], (0, 20.0))

nb4sys = let
    r = sqrt(y[1]^2 + y[2]^2)
    nb4eqs = [
              D(y[1]) ~ -y[2] - (y[1] * y[3]) / r,
              D(y[2]) ~ y[1] - (y[2] * y[3]) / r,
              D(y[3]) ~ y[1] / r,
             ]

    ODESystem(nb4eqs, t, name = :nb4)
end

nb4prob = ODEProblem{false}(structural_simplify(nb4sys), [y[1] => 3.0; y[2:3] .=> 0.0], (0, 20.0))

nb5sys = let
    nb5eqs = [
              D(y[1]) ~ y[2] * y[3],
              D(y[2]) ~ -y[1] * y[3],
              D(y[3]) ~ -0.51 * y[1] * y[2],
             ]

    ODESystem(nb5eqs, t, name = :nb5)
end

nb5prob = ODEProblem{false}(structural_simplify(nb5sys), [y[1] => 0.0; y[2:3] .=> 1.0], (0, 20.0))

nc1sys = let y = y
    n = 10
    y = y[1:n]
    A = Bidiagonal(fill(-1, n), fill(1, n-1), :L)
    nc1eqs = D.(y) .~ A * y
    ODESystem(nc1eqs, t, name = :nc1)
end

nc1prob = ODEProblem{false}(structural_simplify(nc1sys), [y[1] => 1.0; y[2:10] .=> 0.0], (0, 20.0))

nc2sys = let y = y
    n = 10
    y = y[1:n]
    A = Bidiagonal([-1:-1:-n+1; 0], collect(1:n-1), :L)
    nc2eqs = D.(y) .~ A * y
    ODESystem(nc2eqs, t, name = :nc2)
end

nc2prob = ODEProblem{false}(structural_simplify(nc2sys), [y[1] => 1.0; y[2:10] .=> 0.0], (0, 20.0))

nc3sys = let y = y
    n = 10
    y = y[1:n]
    A = SymTridiagonal(fill(-2, n), fill(1, n - 1))
    nc3eqs = D.(y) .~ A * y
    ODESystem(nc3eqs, t, name = :nc3)
end

nc3prob = ODEProblem{false}(structural_simplify(nc3sys), [y[1] => 1.0; y[2:10] .=> 0.0], (0, 20.0))

@variables y(t)[1:51]
y = collect(y)
nc4sys = let y = y
    n = 51
    y = y[1:n]
    A = SymTridiagonal(fill(-2, n), fill(1, n - 1))
    nc4eqs = D.(y) .~ A * y
    ODESystem(nc4eqs, t, name = :nc4)
end

nc4prob = ODEProblem{false}(structural_simplify(nc4sys), [y[1] => 1.0; y[2:51] .=> 0.0], (0, 20.0))

@variables y(t)[1:3, 1:5]
y = collect(y)
nc5sys = let
    k_2 = 2.95912208286
    m_0 = 1.00000597682
    ms = [0.000954786104043
          0.000285583733151
          0.0000437273164546
          0.0000517759138449
          0.00000277777777778]

    r = [sqrt(sum(i->y[i, j]^2, 1:3)) for j in 1:5]
    d = [sqrt(sum(i->(y[i, k] - y[i, j])^2, 1:3)) for k in 1:5, j in 1:5]
    ssum(i, j) = sum(1:5) do k
        k == j && return 0
        ms[k] * (y[i, k] - y[i, j]) / d[j, k]^3
    end
    nc5eqs = [D(D(y[i, j])) ~ k_2 * (-(m_0 + ms[j]) * y[i, j])/r[j]^3 + ssum(i, j) for i in 1:3, j in 1:5]
    structural_simplify(ODESystem(nc5eqs, t, name = :nc5))
end

ys = [3.42947415189, 3.35386959711, 1.35494901715,
      6.64145542550, 5.97156957878, 2.18231499728,
      11.2630437207, 14.6952576794, 6.27960525067,
      -30.1552268759, 1.65699966404, 1.43785752721,
      -21.1238353380, 28.4465098142, 15.3882659679,]
ys′ = [-0.557160570446, 0.505696783289, 0.230578543901,
       -0.415570776342, 0.365682722812, 0.169143213293,
       -0.325325669158, 0.189706021964, 0.087726532278,
       -0.0240476254170, -0.287659532608, -0.117219543175,
       -0.176860753121, -0.216393453025, -0.0148647893090,]
y0 = y .=> reshape(ys, 3, 5)
y0′ = D.(y) .=> reshape(ys′, 3, 5)
# The orginal paper has t_f = 20, but 1000 looks way better
nc5prob = ODEProblem{false}(nc5sys, [y0; y0′], (0, 20.0))

@variables y(t)[1:4]
y = collect(y)
@parameters ε
nd1sys = let
    r = sqrt(y[1]^2 + y[2]^2)^3
    nd1eqs = [D(y[1]) ~ y[3],
              D(y[2]) ~ y[4],
              D(y[3]) ~ (-y[1]) / r,
              D(y[4]) ~ (-y[2]) / r,
             ]
    ODESystem(nd1eqs, t, name = :nd1)
end

function make_ds(nd1sys, e)
    y = collect(@nonamespace nd1sys.y)
    y0 = [y[1] => 1-e; y[2:3] .=> 0.0; y[4] => sqrt((1 + e) / (1 - e))]
    ODEProblem{false}(structural_simplify(nd1sys), y0, (0, 20.0), [ε => e])
end
nd1prob = make_ds(nd1sys, 0.1)
nd2prob = make_ds(nd1sys, 0.3)
nd3prob = make_ds(nd1sys, 0.5)
nd4prob = make_ds(nd1sys, 0.7)
nd5prob = make_ds(nd1sys, 0.9)

ne1sys = let
    ne1eqs = [D(y[1]) ~ y[2],
              D(y[2]) ~ -(y[2] / (t + 1) + (1 - 0.25 / (t + 1)^2) * y[1]),
             ]
    ODESystem(ne1eqs, t, name = :ne1)
end

y0 = [y[1] => 0.6713967071418030; y[2] => 0.09540051444747446]
ne1prob = ODEProblem{false}(structural_simplify(ne1sys), y0, (0, 20.0))

ne2sys = let
    ne2eqs = [D(y[1]) ~ y[2],
              D(y[2]) ~ (1 - y[1]^2) * y[2] - y[1],
             ]
    ODESystem(ne2eqs, t, name = :ne2)
end

y0 = [y[1] => 2.0; y[2] => 0.0]
ne2prob = ODEProblem{false}(structural_simplify(ne2sys), y0, (0, 20.0))

ne3sys = let
    ne3eqs = [D(y[1]) ~ y[2],
              D(y[2]) ~ y[1]^3/6 - y[1] + 2 * sin(2.78535t),
             ]
    ODESystem(ne3eqs, t, name = :ne3)
end

ne3prob = ODEProblem{false}(structural_simplify(ne3sys), y[1:2] .=> 0, (0, 20.0))

ne4sys = let
    ne4eqs = [D(y[1]) ~ y[2],
              D(y[2]) ~ 0.032 - 0.4 * y[2]^2,
             ]
    ODESystem(ne4eqs, t, name = :ne4)
end

ne4prob = ODEProblem{false}(structural_simplify(ne4sys), [y[1] => 30.0, y[2] => 0.0], (0, 20.0))

ne5sys = let
    ne5eqs = [D(y[1]) ~ y[2],
              D(y[2]) ~ sqrt(1 + y[2]^2) / (25 - t),]
    ODESystem(ne5eqs, t, name = :ne5)
end

ne5prob = ODEProblem{false}(structural_simplify(ne5sys), y[1:2] .=> 0.0, (0, 20.0))

cond = term(iseven, term(floor, Int, unwrap(t), type = Int), type = Bool)

nf1sys = let
    a = 0.1
    b = 2a * y[2] - (pi^2 + a^2) * y[1]
    nf1eqs = [D(y[1]) ~ y[2],
              D(y[2]) ~ b + ifelse(cond, 1, -1)]
    ODESystem(nf1eqs, t, name = :nf1)
end

nf1prob = ODEProblem{false}(structural_simplify(nf1sys), y[1:2] .=> 0.0, (0, 20.0))

#=
nf2sys = let
    nf2eqs = [D(y[1]) ~ 55 - ifelse(cond, 2y[1]/2, y[1]/2)]
    ODESystem(nf2eqs, t, name = :nf2)
end

nf2prob = ODEProblem{false}(structural_simplify(nf2sys), [y[1] .=> 110.0], (0, 20.0))
=#

nf3sys = let
    nf3eqs = [D(y[1]) ~ y[2],
              D(y[2]) ~ 0.01 * y[2] * (1 - y[1]^2) - y[1] - abs(sin(pi * t))]
    ODESystem(nf3eqs, t, name = :nf3)
end

nf3prob = ODEProblem{false}(structural_simplify(nf3sys), y[1:2] .=> 0.0, (0, 20.0))

nf4sys = let
    nf4eqs = [D(y[1]) ~ term(myifelse, t <= 10, -2/21 - (120 * (t - 5)) / (1 + 4 * (t - 5)^2), -2y[1], type=Real)]
    ODESystem(nf4eqs, t, name = :nf4)
end

nf4prob = ODEProblem{false}(structural_simplify(nf4sys), [y[1] => 1.0], (0, 20.0))

nf5sys = let
    c = sum(i->cbrt(i)^4, 1:19)
    p = sum(i->cbrt(t - i)^4, 1:19)
    nf5eqs = [D(y[1]) ~ inv(c) * Symbolics.derivative(p, t) * y[1]]
    ODESystem(nf5eqs, t, name = :nf5)
end

nf5prob = ODEProblem{false}(structural_simplify(nf5sys), [y[1] => 1.0], (0, 20.0))

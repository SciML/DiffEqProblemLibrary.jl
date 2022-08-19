using ModelingToolkit
using DiffEqBase, StaticArrays

@variables t y(t)[1:10]
y = collect(y)
D = Differential(t)

sa1sys = let
    sa1eqs = [D(y[1]) ~ -0.5 * y[1], D(y[2]) ~ -y[2], D(y[3]) ~ -100*y[3], D(y[4]) ~ -90*y[4]]
    ODESystem(sa1eqs, t, name = :sa1)
end

sa1prob = ODEProblem{false}(sa1sys, y[1:4] .=> 1.0, (0, 20.0), dt = 1e-2)

sa2sys = let
    sa2eqs = [D(y[1]) ~ -1800 * y[1] + 900 * y[2]]
    for i in 2:8
        push!(sa2eqs, D(y[i]) ~ y[i - 1] - 2 * y[i] + y[i + 1])
    end
    push!(sa2eqs, D(y[9]) ~ 1000 * y[8] - 2000 * y[9] + 1000)
    ODESystem(sa2eqs, t, name = :sa2)
end

sa2prob = ODEProblem{false}(sa2sys, y[1:9] .=> 0.0, (0, 120.0), dt = 5e-4)

sa3sys = let
    sa3eqs = [
              D(y[1]) ~ -1e4 * y[1] + 100 * y[2] - 10 * y[3] + y[4],
              D(y[2]) ~ -1e3 * y[2] + 10 * y[3] - 10 * y[4],
              D(y[3]) ~ -y[3] + 10 * y[4],
              D(y[4]) ~ -0.1 * y[4],
             ]
    ODESystem(sa3eqs, t, name = :sa3)
end

sa3prob = ODEProblem{false}(sa3sys, y[1:4] .=> 1.0, (0, 20.0), dt = 1e-5)

sa4sys = let
    sa4eqs = [D(y[i]) ~ -i^5 * y[i] for i in 1:10]
    ODESystem(sa4eqs, t, name = :sa4)
end

sa4prob = ODEProblem{false}(sa4sys, y[1:10] .=> 1.0, (0, 1.0), dt = 1e-5)

#const SA_PROBLEMS = [sa1prob, sa2prob, sa3prob, sa4prob]

sb1sys = let
    sb1eqs = [D(y[1]) ~ -y[1] + y[2],
              D(y[2]) ~ -100y[1] - y[2],
              D(y[3]) ~ -100y[3] + y[4],
              D(y[4]) ~ -10_000y[3] - 100y[4]]

    ODESystem(sb1eqs, t, name = :sb1)
end

sb1prob = ODEProblem{false}(sb1sys, [y[[1, 3]] .=> 1.0; y[[2, 4]] .=> 0.0;], (0, 20.0), dt = 7e-3)


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

sb2prob = ODEProblem{false}(sb2sys, y .=> 1.0, (0, 20.0), [α => 3], dt = 1e-2)
sb3prob = ODEProblem{false}(sb2sys, y .=> 1.0, (0, 20.0), [α => 8], dt = 1e-2)
sb4prob = ODEProblem{false}(sb2sys, y .=> 1.0, (0, 20.0), [α => 25], dt = 1e-2)
sb5prob = ODEProblem{false}(sb2sys, y .=> 1.0, (0, 20.0), [α => 100], dt = 1e-2)


sc1sys = let
    sc1eqs = [
              D(y[1]) ~ -y[1] + y[2]^2 + y[3]^2 + y[4]^2,
              D(y[2]) ~ -10y[2] + 10*(y[3]^2 + y[4]^2),
              D(y[3]) ~ -40y[3] + 40 * y[4]^2,
              D(y[4]) ~ -100y[4] + 2]

    ODESystem(sc1eqs, t, name = :sc1)
end

sc1prob = ODEProblem{false}(sc1sys, y .=> 1.0, (0, 20.0), dt = 1e-2)


@parameters β
sc2sys = let
    sc2eqs = [D(y[1]) ~ -y[1] + 2,
              D(y[2]) ~ -10y[2] + β * y[1]^2,
              D(y[3]) ~ -40y[3] + 4β * (y[1]^2 + y[2]^2),
              D(y[4]) ~ -100y[4] + 10β * (y[1]^2 + y[2]^2 + y[3]^2)]

    ODESystem(sc2eqs, t, name = :sc2)
end

sc2prob = ODEProblem{false}(sc2sys, y .=> 1.0, (0, 20.0), [β => 0.1], dt = 1e-2, cse = true)
sc3prob = ODEProblem{false}(sc2sys, y .=> 1.0, (0, 20.0), [β => 1.0], dt = 1e-2, cse = true)
sc3prob = ODEProblem{false}(sc2sys, y .=> 1.0, (0, 20.0), [β => 10.0], dt = 1e-2, cse = true)
sc3prob = ODEProblem{false}(sc2sys, y .=> 1.0, (0, 20.0), [β => 20.0], dt = 1e-2, cse = true)

sd1sys = let
    sd1eqs = [D(y[1]) ~ 0.2 * (y[2] - y[1]),
              D(y[2]) ~ 10y[1] - (60 - 0.125y[3])*y[2] + 0.125y[3],
              D(y[3]) ~ 1]

    ODESystem(sd1eqs, t, name = :sd1)
end

sd1prob = ODEProblem{false}(sd1sys, y .=> 0.0, (0, 400.0), [β => 0.1], dt = 1.7e-2)

sd2sys = let
    sd2eqs = [D(y[1]) ~ -0.04y[1] + 0.01 * (y[2] * y[3]),
              D(y[2]) ~ 400y[1] - 100 * (y[2] * y[3]) - 3000 * y[2]^2,
              D(y[3]) ~ 30y[2]^2]

    ODESystem(sd2eqs, t, name = :sd2)
end

sd2prob = ODEProblem{false}(sd2sys, [y[1] => 1.0; y[2:3] .=> 0.0], (0, 40.0), dt = 1e-5, cse = true)

sd3sys = let
    sd3eqs = [D(y[1]) ~ y[3] - 100 * (y[1] * y[2]),
              D(y[2]) ~ y[3] + 2y[4] - 100 * (y[1] * y[2]) - 2e4 * y[2]^2,
              D(y[3]) ~ -y[3] + 100 * (y[1] * y[2]),
              D(y[4]) ~ -y[4] + 1e4*y[2]^2]

    ODESystem(sd3eqs, t, name = :sd3)
end

sd3prob = ODEProblem{false}(sd3sys, [y[1:2] .=> 1; y[3:4] .=> 0.0], (0, 20.0), dt = 2.5e-5, cse = true)

sd4sys = let
    sd4eqs = [D(y[1]) ~ -0.013y[1] - 1000 * (y[1] * y[3]),
              D(y[2]) ~ -2500*(y[2] * y[3]),
              D(y[3]) ~ -0.013y[1] - 1000 * (y[1] * y[3]) - 2500 * (y[2] * y[3])]

    ODESystem(sd4eqs, t, name = :sd4)
end

sd4prob = ODEProblem{false}(sd4sys, [y[1:2] .=> 1; y[3] => 0.0], (0, 50.0), dt = 2.9e-4, cse = true)

sd5sys = let
    sd5eqs = [D(y[1]) ~ 0.01 - (1 + (y[1] + 1000) * (y[1] + 1)) * (0.01 + y[1] + y[2]),
              D(y[2]) ~ 0.01 - (1 + y[2]^2) * (0.01 + y[1] + y[2])]

    ODESystem(sd5eqs, t, name = :sd5)
end

sd5prob = ODEProblem{false}(sd5sys, y[1:2] .=> 0.0, (0, 100.0), dt = 1e-4, cse = true)

sd6sys = let
    sd6eqs = [D(y[1]) ~ -y[1] + 10^8 * y[3] * (1 - y[1]),
              D(y[2]) ~ -10y[2] + 3e7 * y[3] * (1 - y[2]),
              D(y[3]) ~ -(-y[1] + 10^8 * y[3] * (1 - y[1])) - (-10y[2] + 3e7 * y[3] * (1 - y[2])),
             ]

    ODESystem(sd6eqs, t, name = :sd6)
end

sd6prob = ODEProblem{false}(sd6sys, [y[1] => 1.0; y[2:3] .=> 0.0], (0, 1.0), dt = 3.3e-8, cse = true)

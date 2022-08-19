using ModelingToolkit
using DiffEqBase, StaticArrays

@variables t y(t)[1:10]
D = Differential(t)

sa1sys = let
    sa1eqs = [D(y[1]) ~ -0.5 * y[1], D(y[2]) ~ -y[2], D(y[3]) ~ -100*y[3], D(y[4]) ~ -90*y[4]]
    ODESystem(sa1eqs, t, name = :sa1)
end

sa1prob = ODEProblem{false}(sa1sys, collect(y[1:4]) .=> 1.0, (0, 20.0), dt = 1e-2)

sa2sys = let
    sa2eqs = [D(y[1]) ~ -1800 * y[1] + 900 * y[2]]
    for i in 2:8
        push!(sa2eqs, D(y[i]) ~ y[i - 1] - 2 * y[i] + y[i + 1])
    end
    push!(sa2eqs, D(y[9]) ~ 1000 * y[8] - 2000 * y[9] + 1000)
    ODESystem(sa2eqs, t, name = :sa2)
end

sa2prob = ODEProblem{false}(sa2sys, collect(y[1:9]) .=> 0.0, (0, 120.0), dt = 5e-4)

sa3sys = let
    sa3eqs = [
              D(y[1]) ~ -1e4 * y[1] + 100 * y[2] - 10 * y[3] + y[4],
              D(y[2]) ~ -1e3 * y[2] + 10 * y[3] - 10 * y[4],
              D(y[3]) ~ -y[3] + 10 * y[4],
              D(y[4]) ~ -0.1 * y[4],
             ]
    ODESystem(sa3eqs, t, name = :sa3)
end

sa3prob = ODEProblem{false}(sa3sys, collect(y[1:4]) .=> 1.0, (0, 20.0), dt = 1e-5)

sa4sys = let
    sa3eqs = [D(y[i]) ~ -i^5 * y[i] for i in 1:10]
    ODESystem(sa3eqs, t, name = :sa3)
end

sa4prob = ODEProblem{false}(sa4sys, collect(y[1:10]) .=> 1.0, (0, 1.0), dt = 1e-5)

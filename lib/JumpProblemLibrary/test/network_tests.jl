using JumpProblemLibrary, Test

networks = [
    "dnarepressor" => prob_jump_dnarepressor,
    "constproduct" => prob_jump_constproduct,
    "nonlinrxs" => prob_jump_nonlinrxs,
    "osc_mixed_jumptypes" => prob_jump_osc_mixed_jumptypes,
    "multistate" => prob_jump_multistate,
    "twentygenes" => prob_jump_twentygenes,
    "dnadimer_repressor" => prob_jump_dnadimer_repressor,
    "diffnetwork N = 6" => prob_jump_diffnetwork.prob_data["network_func"](6, 2.0),
]

propensities(jpn, u) = jpn.prob_data["jumps"](u, jpn.rates, 0.0)

@testset "$name: network data" for (name, jpn) in networks
    d = jpn.prob_data
    nspecies = length(d["species"])
    nma = length(d["reactant_stoich"])
    njumps = nma + length(d["constant_rate_jumps"])
    @test length(jpn.u0) == nspecies
    @test length(jpn.rates) == nma == length(d["net_stoich"])
    @test jpn.discrete_prob.u0 == jpn.u0
    @test jpn.discrete_prob.tspan == (0.0, jpn.tstop)
    @test length(d["nu"]) == length(d["dep_graph"]) == length(d["jumptovars_map"]) == njumps
    @test length(d["vartojumps_map"]) == nspecies

    all_net_stoich = vcat(d["net_stoich"], last.(d["constant_rate_jumps"]))
    for (nu, ns) in zip(d["nu"], all_net_stoich)
        @test length(nu) == nspecies
        @test nu == [sum((c for (s, c) in ns if s == i); init = 0) for i in 1:nspecies]
    end
    @test d["jumptovars_map"] == [findall(!iszero, nu) for nu in d["nu"]]

    # The dependency graphs are checked against the propensity functions themselves:
    # perturbing species s changes only propensities in vartojumps_map[s], and firing
    # jump j changes only those in dep_graph[j]. Taking the union over several states
    # guards against a perturbation that leaves a propensity unchanged by coincidence.
    states = [convert.(eltype(jpn.u0), 50 .+ mod.(37 .* (1:nspecies) .+ 11m^2, 101)) for m in 1:5]
    for s in 1:nspecies
        changed = map(states) do u
            up = copy(u)
            up[s] += 1
            findall(propensities(jpn, up) .!= propensities(jpn, u))
        end
        @test all(c -> c ⊆ d["vartojumps_map"][s], changed)
        @test sort!(union(changed...)) == d["vartojumps_map"][s]
    end
    for j in 1:njumps
        changed = map(states) do u
            findall(propensities(jpn, u .+ d["nu"][j]) .!= propensities(jpn, u))
        end
        @test all(c -> c ⊆ d["dep_graph"][j], changed)
        @test sort!(union([j], changed...)) == d["dep_graph"][j]
    end
    @test all(u -> length(propensities(jpn, u)) == njumps && all(>(0), propensities(jpn, u)), states)
end

@testset "Propensities" begin
    k = prob_jump_dnarepressor.rates
    u = [1, 5, 7, 2]
    @test propensities(prob_jump_dnarepressor, u) ≈
        [k[1] * 1, k[2] * 5, k[3] * 5, k[4] * 7, k[5] * 1 * 7, k[6] * 2]

    @test propensities(prob_jump_constproduct, [3]) ≈ [1000.0, 10.0 * 3]
    @test prob_jump_constproduct.prob_data["expected_mean_at_t"](1.0) ≈ 100 * (1 - exp(-10))

    k = prob_jump_nonlinrxs.rates
    @test propensities(prob_jump_nonlinrxs, [200, 100, 150]) ≈ [
        k[1] * 200 * 199 / 2, k[2] * 100, k[3] * 200 * 100, k[4] * 150,
        k[5] * 150 * 149 * 148 / 6,
    ]

    hill(x, v, K, n) = v * x^n / (K^n + x^n)
    X, Y, Z, R, S, SP, SP2 = prob_jump_osc_mixed_jumptypes.u0
    @test propensities(prob_jump_osc_mixed_jumptypes, prob_jump_osc_mixed_jumptypes.u0) ≈ [
        0.01 * X, 0.01 * Y, 0.01 * Z, 20.0, 0.01 * SP * (SP - 1) / 2, 0.05 * SP2,
        hill(X, 3.0, 100.0, -4), hill(Y, 3.0, 100.0, -4), hill(Z, 4.5, 100.0, -4),
        hill(X, 2.0, 100.0, 6), hill(Y, 15.0, 100.0, 4) * 0.002 * R, 0.005 * R * S,
    ]
    @test propensities(prob_jump_osc_mixed_jumptypes, [0.0, 0.0, 0.0, R, S, SP, SP2])[7:10] ==
        [3.0, 3.0, 4.5, 0.0]

    k = prob_jump_dnadimer_repressor.rates
    G, M, P, P2, P2G = 1000, 10, 20, 5, 3
    @test propensities(prob_jump_dnadimer_repressor, [G, M, P, P2, P2G]) ≈ [
        k[1] * G, k[2] * M, k[3] * M, k[4] * P, k[5] * P * (P - 1) / 2, k[6] * P2,
        k[7] * P2 * G, k[8] * P2G,
    ]

    kon, koff, kAon, kAoff, kAp, kAdp = 0.01, 0.1, 0.01, 0.1, 0.01, 0.1
    S = collect(1:9) .+ 10
    @test propensities(prob_jump_multistate, S) ≈ [
        kon * S[1] * S[2], kAon * S[1] * S[3], kon * S[2] * S[5], koff * S[4],
        kAon * S[3] * S[4], kAoff * S[5], koff * S[6], kAoff * S[6], kAp * S[6],
        koff * S[7], kAoff * S[7], kAdp * S[7], kon * S[2] * S[8], kAon * S[1] * S[9],
        kAon * S[4] * S[9], kAoff * S[8], kAdp * S[8], kAdp * S[9],
    ]
    @test prob_jump_multistate.u0 == [5360, 1160, 5360, 0, 0, 0, 0, 0, 0]

    d = prob_jump_twentygenes.prob_data
    @test length(d["species"]) == 80
    @test prob_jump_twentygenes.u0 == [fill(1, 20); fill(0, 60)]
    sp = Dict(s => i for (i, s) in enumerate(d["species"]))
    u = collect(1:80)
    a = propensities(prob_jump_twentygenes, u)
    @test length(a) == 100
    for i in 1:10, (g, kg) in ((2i - 1, 10.0), (2i, 5.0))
        @test kg * u[sp[Symbol(:G, g)]] in a
    end
    for i in 1:10
        cross = 0.0001 * u[sp[Symbol(:G, 2i)]] * u[sp[Symbol(:P, 2i - 1)]]
        @test a[10i - 1] ≈ cross
        @test d["net_stoich"][10i - 1] ==
            sort([sp[Symbol(:G, 2i)] => -1, sp[Symbol(:P, 2i - 1)] => -1, sp[Symbol(:G_ind, 2i)] => 1])
    end

    diff3 = prob_jump_diffnetwork.prob_data["network_func"](3, 2.0)
    @test diff3.u0 == prob_jump_diffnetwork.u0(3) == [10, 10, 10]
    @test propensities(diff3, [1, 2, 3]) == 2.0 .* [1, 2, 2, 3]
    @test diff3.prob_data["nu"] == [[-1, 1, 0], [1, -1, 0], [0, -1, 1], [0, 1, -1]]
end

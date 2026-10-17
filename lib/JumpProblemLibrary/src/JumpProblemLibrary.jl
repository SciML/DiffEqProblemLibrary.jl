module JumpProblemLibrary

using SciMLBase: DiscreteProblem

# Jump Example Problems
export prob_jump_dnarepressor, prob_jump_constproduct, prob_jump_nonlinrxs,
    # examples mixing mass action and constant rate jumps
    prob_jump_osc_mixed_jumptypes,
    # examples used in published benchmarks / comparisons
    prob_jump_multistate, prob_jump_twentygenes, prob_jump_dnadimer_repressor,
    # examples approximating diffusion by continuous time random walks
    prob_jump_diffnetwork

"""
    JumpProblemNetwork

Container for the inputs needed to construct a `JumpProblem`. It does not store a
`JumpProblem` because that requires choosing an aggregator.

Fields:

  - `rates`: rate constants of the mass action reactions, one per reaction.
  - `tstop`: final simulation time.
  - `u0`: initial species counts, ordered as `prob_data["species"]`.
  - `discrete_prob`: `DiscreteProblem(u0, (0.0, tstop), rates)`.
  - `prob_data`: a `Dict` holding the reaction network and any reference results.

Every network stores these `prob_data` entries, with jumps ordered mass action
reactions first, then constant rate jumps:

  - `"species"`: species names, in state-vector order.
  - `"reactant_stoich"`, `"net_stoich"`: stoichiometry of the mass action reactions as
    vectors of `species index => coefficient` pairs (the `MassActionJump` format).
  - `"constant_rate_jumps"`: vector of `(rate, net_stoich)` tuples for the reactions
    whose propensity is not of mass action form, with `rate(u, p, t)` the propensity.
  - `"jumps"`: `jumps(u, p, t)` returns the vector of all propensities, with
    mass action propensities `k * ∏ binomial(u[s], c)`.
  - `"nu"`: dense state change vector of each jump.
  - `"dep_graph"`: `dep_graph[j]` lists the jumps whose propensity changes when jump
    `j` fires, including `j`.
  - `"vartojumps_map"`: `vartojumps_map[s]` lists the jumps whose propensity
    depends on species `s`.
  - `"jumptovars_map"`: `jumptovars_map[j]` lists the species changed by jump `j`.

Example, using JumpProcesses:

```julia
using JumpProcesses
jpn = prob_jump_dnarepressor
d = jpn.prob_data
majumps = MassActionJump(jpn.rates, d["reactant_stoich"], d["net_stoich"])
jprob = JumpProblem(jpn.discrete_prob, NRM(), majumps; dep_graph = d["dep_graph"])
```

`prob_jump_diffnetwork` instead stores functions of the lattice size, see its docstring.
"""
struct JumpProblemNetwork
    rates::Any
    tstop::Any
    u0::Any
    discrete_prob::Any
    prob_data::Any
end

"""
Convert reactions given as `reactants => products`, with each side a vector of
`species => coefficient` pairs, into `MassActionJump` reactant and net stoichiometry.
"""
function stoichiometry(species, reactions)
    idx = Dict(s => i for (i, s) in enumerate(species))
    toidx(side) = sort!(Pair{Int, Int}[idx[s] => c for (s, c) in side])
    reactant_stoich = [toidx(r) for (r, _) in reactions]
    net_stoich = map(reactions) do (r, p)
        net = Dict{Int, Int}()
        for (s, c) in r
            net[idx[s]] = get(net, idx[s], 0) - c
        end
        for (s, c) in p
            net[idx[s]] = get(net, idx[s], 0) + c
        end
        sort!(Pair{Int, Int}[s => c for (s, c) in net if c != 0])
    end
    return reactant_stoich, net_stoich
end

function massaction_propensity(u, k, reactants)
    a = float(k)
    for (s, c) in reactants
        a *= prod(u[s] - i for i in 0:(c - 1)) / factorial(c)
    end
    return a
end

"""
Assemble the `prob_data` entries of a network. Each constant
rate jump is a `(rate, net_stoich, rate_species)` tuple, with `rate_species` the
indices of the species its propensity depends on.
"""
function network_data(
        species, reactant_stoich, net_stoich;
        constant_rate_jumps = Tuple{Function, Vector{Pair{Int, Int}}, Vector{Int}}[]
    )
    nspecies = length(species)
    crj_rates = [f for (f, _, _) in constant_rate_jumps]
    all_net_stoich = vcat(net_stoich, [ns for (_, ns, _) in constant_rate_jumps])
    rate_species = vcat(
        [first.(rs) for rs in reactant_stoich], [deps for (_, _, deps) in constant_rate_jumps]
    )

    function jumps(u, p, t)
        props = [massaction_propensity(u, p[j], rs) for (j, rs) in enumerate(reactant_stoich)]
        return vcat(props, [f(u, p, t) for f in crj_rates])
    end

    nu = map(all_net_stoich) do ns
        v = zeros(Int, nspecies)
        for (s, c) in ns
            v[s] = c
        end
        v
    end

    vartojumps_map = [Int[] for _ in 1:nspecies]
    for (j, specs) in enumerate(rate_species), s in specs
        push!(vartojumps_map[s], j)
    end
    jumptovars_map = [first.(ns) for ns in all_net_stoich]
    dep_graph = [
        sort!(union([j], (vartojumps_map[s] for s in jumptovars_map[j])...))
            for j in eachindex(all_net_stoich)
    ]

    return Dict{String, Any}(
        "species" => species,
        "reactant_stoich" => reactant_stoich, "net_stoich" => net_stoich,
        "constant_rate_jumps" => [(f, ns) for (f, ns, _) in constant_rate_jumps],
        "jumps" => jumps, "nu" => nu, "dep_graph" => dep_graph,
        "vartojumps_map" => vartojumps_map, "jumptovars_map" => jumptovars_map
    )
end

function build_network(rates, tstop, u0, species, reactions, extra_data = Dict(); kwargs...)
    data = network_data(species, stoichiometry(species, reactions)...; kwargs...)
    merge!(data, extra_data)
    prob = DiscreteProblem(u0, (0.0, tstop), rates)
    return JumpProblemNetwork(rates, tstop, u0, prob, data)
end

const ∅ = Pair{Symbol, Int}[]

dna_species = [:DNA, :mRNA, :P, :DNAR]
dna_reactions = [
    [:DNA => 1] => [:mRNA => 1, :DNA => 1],
    [:mRNA => 1] => [:mRNA => 1, :P => 1],
    [:mRNA => 1] => ∅,
    [:P => 1] => ∅,
    [:DNA => 1, :P => 1] => [:DNAR => 1],
    [:DNAR => 1] => [:DNA => 1, :P => 1],
]
rates = [0.5, (20 * log(2.0) / 120.0), (log(2.0) / 120.0), (log(2.0) / 600.0), 0.025, 1.0]
"""
DNA negative feedback autoregulatory model. Protein acts as repressor.

Species `[DNA, mRNA, P, DNAR]`, reactions:

 1. `k1, DNA --> mRNA + DNA`
 2. `k2, mRNA --> mRNA + P`
 3. `k3, mRNA --> 0`
 4. `k4, P --> 0`
 5. `k5, DNA + P --> DNAR`
 6. `k6, DNAR --> DNA + P`
"""
prob_jump_dnarepressor = build_network(
    rates, 1000.0, [1, 0, 0, 0], dna_species, dna_reactions,
    Dict("num_sims_for_mean" => 8000, "expected_mean" => 5.92655375e+2)
)

bd_rates = [1000.0, 10.0]
expected_avg = t -> bd_rates[1] / bd_rates[2] .* (1.0 .- exp.(-bd_rates[2] * t))
"""
Simple birth-death process with constant production and degradation.

Species `[A]`, reactions:

 1. `k1, 0 --> A`
 2. `k2, A --> 0`
"""
prob_jump_constproduct = build_network(
    bd_rates, 1.0, [0], [:A], [∅ => [:A => 1], [:A => 1] => ∅],
    Dict("num_sims_for_mean" => 16000, "expected_mean_at_t" => expected_avg)
)

nonlin_reactions = [
    [:A => 2] => [:B => 1],
    [:B => 1] => [:A => 2],
    [:A => 1, :B => 1] => [:C => 1],
    [:C => 1] => [:A => 1, :B => 1],
    [:C => 3] => [:A => 3],
]
"""
Example with a mix of nonlinear reactions, including third order.

Species `[A, B, C]`, reactions:

 1. `k1, 2A --> B`
 2. `k2, B --> 2A`
 3. `k3, A + B --> C`
 4. `k4, C --> A + B`
 5. `k5, 3C --> 3A`
"""
prob_jump_nonlinrxs = build_network(
    [1.0, 2.0, 0.5, 0.75, 0.25], 0.01, [200, 100, 150], [:A, :B, :C], nonlin_reactions,
    Dict("num_sims_for_mean" => 32000, "expected_mean" => 84.876015624999994)
)

"""
    hill(x, v, K, n)

Hill function `v * x^n / (K^n + x^n)`, written in a form that is also finite at `x = 0`.
"""
hill(x, v, K, n) = v / (1 + (K / x)^n)

oscil_species = [:X, :Y, :Z, :R, :S, :SP, :SP2]
oscil_reactions = [
    [:X => 1] => ∅,
    [:Y => 1] => ∅,
    [:Z => 1] => ∅,
    ∅ => [:S => 1],
    [:SP => 2] => [:SP2 => 1],
    [:SP2 => 1] => ∅,
]
oscil_crjs = [
    ((u, p, t) -> hill(u[1], 3.0, 100.0, -4), [2 => 1], [1]),
    ((u, p, t) -> hill(u[2], 3.0, 100.0, -4), [3 => 1], [2]),
    ((u, p, t) -> hill(u[3], 4.5, 100.0, -4), [1 => 1], [3]),
    ((u, p, t) -> hill(u[1], 2.0, 100.0, 6), [4 => 1], [1]),
    ((u, p, t) -> hill(u[2], 15.0, 100.0, 4) * 0.002 * u[4], [4 => -1], [2, 4]),
    ((u, p, t) -> 0.005 * u[4] * u[5], [5 => -1, 6 => 1], [4, 5]),
]
"""
Oscillatory system, uses a mixture of jump types.

Species `[X, Y, Z, R, S, SP, SP2]`. Mass action reactions:

 1. `0.01, X --> 0`
 2. `0.01, Y --> 0`
 3. `0.01, Z --> 0`
 4. `20, 0 --> S`
 5. `0.01, SP + SP --> SP2`
 6. `0.05, SP2 --> 0`

Constant rate jumps, with `hill(x, v, K, n) = v * x^n / (K^n + x^n)`:

 7. `hill(X, 3.0, 100.0, -4), 0 --> Y`
 8. `hill(Y, 3.0, 100.0, -4), 0 --> Z`
 9. `hill(Z, 4.5, 100.0, -4), 0 --> X`
 10. `hill(X, 2.0, 100.0, 6), 0 --> R`
 11. `hill(Y, 15.0, 100.0, 4) * 0.002, R --> 0`
 12. `R * 0.005, S --> SP`

Species counts are floats since the Hill functions are evaluated on them.
"""
prob_jump_osc_mixed_jumptypes = build_network(
    [0.01, 0.01, 0.01, 20.0, 0.01, 0.05], 4000.0,
    [200.0, 60.0, 120.0, 100.0, 50.0, 50.0, 50.0], oscil_species, oscil_reactions;
    constant_rate_jumps = oscil_crjs
)

specs_sym_to_name = Dict(
    :S1 => "R(a,l)",
    :S2 => "L(r)",
    :S3 => "A(Y~U,r)",
    :S4 => "L(r!1).R(a,l!1)",
    :S5 => "A(Y~U,r!1).R(a!1,l)",
    :S6 => "A(Y~U,r!1).L(r!2).R(a!1,l!2)",
    :S7 => "A(Y~P,r!1).L(r!2).R(a!1,l!2)",
    :S8 => "A(Y~P,r!1).R(a!1,l)",
    :S9 => "A(Y~P,r)"
)
rsi = Dict(
    :R0 => 1, :L0 => 2, :A0 => 3, :kon => 4, :koff => 5,
    :kAon => 6, :kAoff => 7, :kAp => 8, :kAdp => 9
)
params = (5360, 1160, 5360, 0.01, 0.1, 0.01, 0.1, 0.01, 0.1)
multistate_reactions = [
    :kon => [:S1 => 1, :S2 => 1] => [:S4 => 1],
    :kAon => [:S1 => 1, :S3 => 1] => [:S5 => 1],
    :kon => [:S2 => 1, :S5 => 1] => [:S6 => 1],
    :koff => [:S4 => 1] => [:S1 => 1, :S2 => 1],
    :kAon => [:S3 => 1, :S4 => 1] => [:S6 => 1],
    :kAoff => [:S5 => 1] => [:S1 => 1, :S3 => 1],
    :koff => [:S6 => 1] => [:S2 => 1, :S5 => 1],
    :kAoff => [:S6 => 1] => [:S3 => 1, :S4 => 1],
    :kAp => [:S6 => 1] => [:S7 => 1],
    :koff => [:S7 => 1] => [:S2 => 1, :S8 => 1],
    :kAoff => [:S7 => 1] => [:S4 => 1, :S9 => 1],
    :kAdp => [:S7 => 1] => [:S6 => 1],
    :kon => [:S2 => 1, :S8 => 1] => [:S7 => 1],
    :kAon => [:S1 => 1, :S9 => 1] => [:S8 => 1],
    :kAon => [:S4 => 1, :S9 => 1] => [:S7 => 1],
    :kAoff => [:S8 => 1] => [:S1 => 1, :S9 => 1],
    :kAdp => [:S8 => 1] => [:S5 => 1],
    :kAdp => [:S9 => 1] => [:S3 => 1],
]
"""
Multistate model from Gupta and Mendes,
"An Overview of Network-Based and -Free Approaches for Stochastic Simulation of Biochemical Systems",
Computation 2018, 6, 9; doi:10.3390/computation6010009
Translated from supplementary data file: Models/Multi-state/fixed_multistate.xml

Species `[S1, ..., S9]`, reactions:

 1. `kon, S1 + S2 --> S4`
 2. `kAon, S1 + S3 --> S5`
 3. `kon, S2 + S5 --> S6`
 4. `koff, S4 --> S1 + S2`
 5. `kAon, S3 + S4 --> S6`
 6. `kAoff, S5 --> S1 + S3`
 7. `koff, S6 --> S2 + S5`
 8. `kAoff, S6 --> S3 + S4`
 9. `kAp, S6 --> S7`
 10. `koff, S7 --> S2 + S8`
 11. `kAoff, S7 --> S4 + S9`
 12. `kAdp, S7 --> S6`
 13. `kon, S2 + S8 --> S7`
 14. `kAon, S1 + S9 --> S8`
 15. `kAon, S4 + S9 --> S7`
 16. `kAoff, S8 --> S1 + S9`
 17. `kAdp, S8 --> S5`
 18. `kAdp, S9 --> S3`

`prob_data["params"]` holds the initial counts and rate constants, indexed by
`prob_data["rates_sym_to_idx"]`.
"""
prob_jump_multistate = build_network(
    [params[rsi[k]] for (k, _) in multistate_reactions], 100.0,
    [params[1], params[2], params[3], 0, 0, 0, 0, 0, 0], Symbol.(:S, 1:9),
    last.(multistate_reactions),
    Dict(
        "specs_to_sym_name" => specs_sym_to_name,
        "rates_sym_to_idx" => rsi, "params" => params
    )
)

function twentygenes_reactions(N)
    G(i) = Symbol(:G, i)
    M(i) = Symbol(:M, i)
    P(i) = Symbol(:P, i)
    G_ind(i) = Symbol(:G_ind, i)
    rates = Float64[]
    reactions = Pair{Vector{Pair{Symbol, Int}}, Vector{Pair{Symbol, Int}}}[]
    for i in 1:N
        for (g, k) in ((2i - 1, 10.0), (2i, 5.0))
            append!(rates, [k, k, 1.0, 1.0])
            push!(reactions, [G(g) => 1] => [G(g) => 1, M(g) => 1])
            push!(reactions, [M(g) => 1] => [M(g) => 1, P(g) => 1])
            push!(reactions, [M(g) => 1] => ∅)
            push!(reactions, [P(g) => 1] => ∅)
        end
        append!(rates, [0.0001, 100.0])
        push!(reactions, [G(2i) => 1, P(2i - 1) => 1] => [G_ind(2i) => 1])
        push!(reactions, [G_ind(2i) => 1] => [G_ind(2i) => 1, M(2i) => 1])
    end
    species = [f(i) for f in (G, M, P, G_ind) for i in 1:(2N)]
    return rates, species, reactions
end
rates, species, reactions = twentygenes_reactions(10)
"""
Twenty-gene model from McCollum et al,
"The sorting direct method for stochastic simulation of biochemical systems with varying reaction execution behavior"
Comp. Bio. and Chem., 30, pg. 39-49 (2006).

Species `[G1, ..., G20, M1, ..., M20, P1, ..., P20, G_ind1, ..., G_ind20]`. For each
gene `g` with rate `k = 10` for odd and `k = 5` for even `g`:

  - `k, G_g --> G_g + M_g`
  - `k, M_g --> M_g + P_g`
  - `1, M_g --> 0`
  - `1, P_g --> 0`

and for each gene pair `i = 1, ..., 10`:

  - `0.0001, G_(2i) + P_(2i-1) --> G_ind_(2i)`
  - `100, G_ind_(2i) --> G_ind_(2i) + M_(2i)`
"""
prob_jump_twentygenes = build_network(
    rates, 2000.0, [i <= 20 ? 1 : 0 for i in eachindex(species)], species, reactions
)

dimer_species = [:G, :M, :P, :P2, :P2G]
dimer_reactions = [
    [:G => 1] => [:G => 1, :M => 1],
    [:M => 1] => [:M => 1, :P => 1],
    [:M => 1] => ∅,
    [:P => 1] => ∅,
    [:P => 2] => [:P2 => 1],
    [:P2 => 1] => [:P => 2],
    [:P2 => 1, :G => 1] => [:P2G => 1],
    [:P2G => 1] => [:P2 => 1, :G => 1],
]
"""
Negative feedback autoregulatory gene expression model. Dimer is the repressor.
Taken from Marchetti, Priami and Thanh,
"Simulation Algorithms for Computational Systems Biology",
Springer (2017).

Species `[G, M, P, P2, P2G]`, reactions:

 1. `c1, G --> G + M`
 2. `c2, M --> M + P`
 3. `c3, M --> 0`
 4. `c4, P --> 0`
 5. `c5, 2P --> P2`
 6. `c6, P2 --> 2P`
 7. `c7, P2 + G --> P2G`
 8. `c8, P2G --> P2 + G`
"""
prob_jump_dnadimer_repressor = build_network(
    [0.09, 0.05, 0.001, 0.0009, 0.00001, 0.0005, 0.005, 0.9], 4000.0,
    [1000, 0, 0, 0, 0], dimer_species, dimer_reactions,
    Dict("specs_names" => string.(dimer_species))
)

function getDiffNetwork(N, K = 1.0)
    species = Symbol.(:X, 1:N)
    reactions = Pair{Vector{Pair{Symbol, Int}}, Vector{Pair{Symbol, Int}}}[]
    for i in 1:(N - 1)
        push!(reactions, [species[i] => 1] => [species[i + 1] => 1])
        push!(reactions, [species[i + 1] => 1] => [species[i] => 1])
    end
    rates = fill(float(K), length(reactions))
    return build_network(rates, 10.0, getDiffu0(N), species, reactions)
end
getDiffu0(N) = fill(10, N)
"""
Continuous time random walk (i.e. diffusion approximation) example on a 1D lattice
of `N` sites with hopping rate `K` between neighboring sites:

  - `K, X_i --> X_(i+1)`
  - `K, X_(i+1) --> X_i`

for `i = 1, ..., N - 1`, starting with 10 particles per site. As the network depends
on `N`, `prob_data["network_func"](N, K)` returns the full `JumpProblemNetwork` for a
lattice of `N` sites, and `u0(N)` returns its initial condition.
"""
prob_jump_diffnetwork = JumpProblemNetwork(
    [1.0], 10.0, getDiffu0, nothing, Dict("network_func" => getDiffNetwork)
)

end # module

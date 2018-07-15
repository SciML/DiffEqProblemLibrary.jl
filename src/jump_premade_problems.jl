using DiffEqBase, DiffEqBiological
# Jump Example Problems
export prob_jump_dnarepressor, prob_jump_constproduct, prob_jump_nonlinrxs,
# examples mixing mass action and constant rate jumps
       prob_jump_osc_mixed_jumptypes,
# examples used in published benchmarks / comparisions
       prob_jump_multistate, prob_jump_twentygenes, prob_jump_dnadimer_repressor,
# examples approximating diffusion by continuous time random walks
       prob_jump_diffnetwork

"""
    General structure to hold JumpProblem info. Needed since
    the JumpProblem constructor requires the algorithm, so we
    don't create the JumpProblem here.
"""
struct JumpProblemNetwork
    network         # DiffEqBiological network
    rates           # vector of rate constants or nothing
    tstop           # time to end simulation
    u0              # initial values
    discrete_prob   # initialized discrete problem
    prob_data       # additional problem data, stored as a Dict
end

dna_rs = @reaction_network ptype begin
    k1, DNA --> mRNA + DNA
    k2, mRNA --> mRNA + P
    k3, mRNA --> 0
    k4, P --> 0
    k5, DNA + P --> DNAR
    k6, DNAR --> DNA + P
end k1 k2 k3 k4 k5 k6
rates = [.5, (20*log(2.)/120.), (log(2.)/120.), (log(2.)/600.), .025, 1.]
tf = 1000.0
u0 = [1,0,0,0]
prob = DiscreteProblem(u0, (0.0, tf), rates)
Nsims = 8000
expected_avg = 5.926553750000000e+02
prob_data = Dict("num_sims_for_mean" => Nsims, "expected_mean" => expected_avg)
"""
    DNA negative feedback autoregulatory model. Protein acts as repressor.
"""
prob_jump_dnarepressor = JumpProblemNetwork(dna_rs, rates, tf, u0, prob, prob_data)

bd_rs = @reaction_network pdtype begin
    k1, 0 --> A
    k2, A --> 0
end k1 k2
rates = [1000., 10.]
tf = 1.0
u0 = [0]
prob = DiscreteProblem(u0, (0., tf), rates)
Nsims = 16000
expected_avg = t -> rates[1] / rates[2] .* ( 1. - exp.(-rates[2] * t))
prob_data = Dict("num_sims_for_mean" => Nsims, "expected_mean_at_t" => expected_avg)
"""
    Simple birth-death process with constant production and degradation.
"""
prob_jump_constproduct = JumpProblemNetwork(bd_rs, rates, tf, u0, prob, prob_data)

nonlin_rs = @reaction_network dtype begin
    k1, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
end k1 k2 k3 k4 k5
rates = [1., 2., .5, .75, .25]
tf = .01
u0 = [200, 100, 150]
prob = DiscreteProblem(u0, (0., tf), rates)
Nsims = 32000
expected_avg = 84.876015624999994
prob_data = Dict("num_sims_for_mean" => Nsims, "expected_mean" => expected_avg)
"""
    Example with a mix of nonlinear reactions, including third order
"""
prob_jump_nonlinrxs = JumpProblemNetwork(nonlin_rs, rates, tf, u0, prob, prob_data)


oscil_rs = @reaction_network rnoscType  begin
    0.01, (X,Y,Z) --> 0
    hill(X,3.,100.,-4), 0 --> Y
    hill(Y,3.,100.,-4), 0 --> Z
    hill(Z,4.5,100.,-4), 0 --> X
    hill(X,2.,100.,6), 0 --> R
    hill(Y,15.,100.,4)*0.002, R --> 0
    20, 0 --> S
    R*0.005, S --> SP
    0.01, SP + SP --> SP2
    0.05, SP2 --> 0
end
u0 = [200.,60.,120.,100.,50.,50.,50.]  # Hill equations force use of floats!
tf = 4000.
prob = DiscreteProblem(u0, (0.,tf))
"""
    Oscillatory system, uses a mixture of jump types.
"""
prob_jump_osc_mixed_jumptypes = JumpProblemNetwork(oscil_rs, nothing, tf, u0, prob, nothing)


specs_sym_to_name = Dict(
    :S1 => "R(a,l)",
    :S2 => "L(r)",
    :S3 => "A(Y~U,r)",
    :S4 => "L(r!1).R(a,l!1)",
    :S5 => "A(Y~U,r!1).R(a!1,l)",
    :S6 => "A(Y~U,r!1).L(r!2).R(a!1,l!2)",
    :S7 => "A(Y~P,r!1).L(r!2).R(a!1,l!2)",
    :S8 => "A(Y~P,r!1).R(a!1,l)",
    :S9 => "A(Y~P,r)")
rates_sym_to_idx = Dict(
    :R0 => 1,   :L0 => 2,    :A0 => 3,  :kon => 4, :koff => 5,
    :kAon => 6, :kAoff => 7, :kAp => 8, :kAdp => 9)
params = [5360, 1160, 5360, 0.01, 0.1, 0.01, 0.1, 0.01, 0.1]
rs = @reaction_network msType begin
    kon, S1 + S2 --> S4
    kAon, S1 + S3 --> S5
    kon, S2 + S5 --> S6
    koff, S4 --> S1 + S2
    kAon, S3 + S4 --> S6
    kAoff, S5 --> S1 + S3
    koff, S6 --> S2 + S5
    kAoff, S6 --> S3 + S4
    kAp, S6 --> S7
    koff, S7 --> S2 + S8
    kAoff, S7 --> S4 + S9
    kAdp, S7 --> S6
    kon, S2 + S8 --> S7
    kAon, S1 + S9 --> S8
    kAon, S4 + S9 --> S7
    kAoff, S8 --> S1 + S9
    kAdp, S8 --> S5
    kAdp, S9 --> S3
end kon kAon koff kAoff kAp kAdp
rsi   = rates_sym_to_idx
rates = params[[rsi[:kon], rsi[:kAon], rsi[:koff], rsi[:kAoff], rsi[:kAp], rsi[:kAdp]]]
u0    = zeros(Int,9)
u0[ something(findfirst(isequal(:S1), rs.syms),0)] = params[1]
u0[ something(findfirst(isequal(:S2), rs.syms),0)] = params[2]
u0[ something(findfirst(isequal(:S3), rs.syms),0)] = params[3]
tf    = 100.
prob = DiscreteProblem(u0, (0., tf), rates)
"""
    Multistate model from Gupta and Mendes,
    "An Overview of Network-Based and -Free Approaches for Stochastic Simulation of Biochemical Systems",
    Computation 2018, 6, 9; doi:10.3390/computation6010009
    Translated from supplementary data file: Models/Multi-state/fixed_multistate.xml
"""
prob_jump_multistate = JumpProblemNetwork(rs, rates, tf, u0, prob,
                Dict("specs_to_sym_name" => specs_sym_to_name, "rates_sym_to_idx" => rates_sym_to_idx, "params" => params))


# generate the network
N = 10  # number of genes
function construct_genenetwork(N)
  genenetwork = "@reaction_network twentgtype begin\n"
  for i in 1:N
      genenetwork *= "\t 10.0, G$(2*i-1) --> G$(2*i-1) + M$(2*i-1)\n"
      genenetwork *= "\t 10.0, M$(2*i-1) --> M$(2*i-1) + P$(2*i-1)\n"
      genenetwork *= "\t 1.0,  M$(2*i-1) --> 0\n"
      genenetwork *= "\t 1.0,  P$(2*i-1) --> 0\n"

      genenetwork *= "\t 5.0, G$(2*i) --> G$(2*i) + M$(2*i)\n"
      genenetwork *= "\t 5.0, M$(2*i) --> M$(2*i) + P$(2*i)\n"
      genenetwork *= "\t 1.0,  M$(2*i) --> 0\n"
      genenetwork *= "\t 1.0,  P$(2*i) --> 0\n"

      genenetwork *= "\t 0.0001, G$(2*i) + P$(2*i-1) --> G$(2*i)_ind \n"
      genenetwork *= "\t 100., G$(2*i)_ind --> G$(2*i)_ind + M$(2*i)\n"
  end
  genenetwork *= "end"
end
genenetwork = construct_genenetwork(N)
rs = eval( Meta.parse(genenetwork) )
u0 = zeros(Int, length(rs.syms))
for i = 1:(2*N)
    u0[something(findfirst(isequal(Symbol("G$(i)")),rs.syms),0)] = 1
end
tf = 2000.0
prob = DiscreteProblem(u0, (0.0, tf))
"""
    Twenty-gene model from McCollum et al,
    "The sorting direct method for stochastic simulation of biochemical systems with varying reaction execution behavior"
    Comp. Bio. and Chem., 30, pg. 39-49 (2006).
"""
prob_jump_twentygenes = JumpProblemNetwork(rs, nothing, tf, u0, prob, nothing)


rn = @reaction_network gnrdtype begin
    c1, G --> G + M
    c2, M --> M + P
    c3, M --> 0
    c4, P --> 0
    c5, 2P --> P2
    c6, P2 --> 2P
    c7, P2 + G --> P2G
    c8, P2G --> P2 + G
end c1 c2 c3 c4 c5 c6 c7 c8
rnpar = [.09, .05, .001, .0009, .00001, .0005, .005, .9]
varlabels = ["G", "M", "P", "P2","P2G"]
u0 = [1000, 0, 0, 0,0]
tf = 4000.
prob = DiscreteProblem(u0, (0.0, tf), rnpar)
"""
    Negative feedback autoregulatory gene expression model. Dimer is the repressor.
    Taken from Marchetti, Priami and Thanh,
    "Simulation Algorithms for Comptuational Systems Biology",
    Springer (2017).
"""
prob_jump_dnadimer_repressor = JumpProblemNetwork(rn, rnpar, tf, u0, prob,
                                                Dict("specs_names" => varlabels))


# diffusion model
function getDiffNetwork(N)
    diffnetwork = "@reaction_network dpldifftype begin\n"
    for i in 1:(N-1)
        diffnetwork *= "\t K, X$(i) --> X$(i+1)\n"
        diffnetwork *= "\t K, X$(i+1) --> X$(i)\n"
    end
    diffnetwork *= "end K"
    rs = eval( Meta.parse(diffnetwork) )
    rs
end
params = (1.,)
function getDiffu0(N)
    10*ones(Int64, N)
end
tf = 10.
"""
    Continuous time random walk (i.e. diffusion approximation) example.
    Here the network in the JumpProblemNetwork is a function that returns a
    network given the number of lattice sites.
    u0 is a similar function that returns the initial condition vector.
"""
prob_jump_diffnetwork = JumpProblemNetwork(getDiffNetwork, params, tf, getDiffu0, nothing, nothing)

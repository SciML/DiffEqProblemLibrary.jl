module JumpProblemLibrary

using DiffEqBase

import RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

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
    rates::Any           # vector of rate constants or jumps
    tstop::Any           # time to end simulation
    u0::Any              # initial values
    discrete_prob::Any   # initialized discrete problem
    prob_data::Any       # additional problem data, stored as a Dict
end

# DNA repressor model - direct jump implementation
# Species: DNA=1, mRNA=2, P=3, DNAR=4
dna_rates = [0.5, (20 * log(2.0) / 120.0), (log(2.0) / 120.0), (log(2.0) / 600.0), 0.025, 1.0]

# Define jump processes for DNA repressor model
function dna_repressor_jumps(u, p, t)
    DNA, mRNA, P, DNAR = u
    k1, k2, k3, k4, k5, k6 = p
    [
        k1 * DNA,           # k1: DNA --> mRNA + DNA (mRNA production)
        k2 * mRNA,          # k2: mRNA --> mRNA + P (protein production) 
        k3 * mRNA,          # k3: mRNA --> 0 (mRNA degradation)
        k4 * P,             # k4: P --> 0 (protein degradation)
        k5 * DNA * P,       # k5: DNA + P --> DNAR (repressor binding)
        k6 * DNAR           # k6: DNAR --> DNA + P (repressor unbinding)
    ]
end

# State change vectors for each reaction
dna_nu = [
    [0, 1, 0, 0],    # k1: mRNA production
    [0, 0, 1, 0],    # k2: protein production  
    [0, -1, 0, 0],   # k3: mRNA degradation
    [0, 0, -1, 0],   # k4: protein degradation
    [-1, 0, -1, 1],  # k5: repressor binding
    [1, 0, 1, -1]    # k6: repressor unbinding
]

tf = 1000.0
u0 = [1, 0, 0, 0]  # [DNA, mRNA, P, DNAR]
prob = DiscreteProblem(u0, (0.0, tf), dna_rates)
Nsims = 8000
expected_avg = 5.926553750000000e+02
prob_data = Dict("num_sims_for_mean" => Nsims, "expected_mean" => expected_avg, 
                 "jumps" => dna_repressor_jumps, "nu" => dna_nu)
"""
    DNA negative feedback autoregulatory model. Protein acts as repressor.
    
    Reactions:
    1. DNA --> mRNA + DNA (mRNA production, rate k1=0.5)
    2. mRNA --> mRNA + P (protein production, rate k2≈0.115)  
    3. mRNA --> 0 (mRNA degradation, rate k3≈0.0058)
    4. P --> 0 (protein degradation, rate k4≈0.0012)
    5. DNA + P --> DNAR (repressor binding, rate k5=0.025)
    6. DNAR --> DNA + P (repressor unbinding, rate k6=1.0)
"""
prob_jump_dnarepressor = JumpProblemNetwork(dna_rates, tf, u0, prob, prob_data)

# Birth-death process - direct jump implementation  
# Species: A=1
bd_rates = [1000.0, 10.0]

function bd_jumps(u, p, t)
    A = u[1] 
    k1, k2 = p
    [
        k1,      # k1: 0 --> A (birth)
        k2 * A   # k2: A --> 0 (death)
    ]
end

bd_nu = [
    [1],   # k1: birth
    [-1]   # k2: death  
]

tf = 1.0
u0 = [0]  # [A]
prob = DiscreteProblem(u0, (0.0, tf), bd_rates)
Nsims = 16000
expected_avg = t -> bd_rates[1] / bd_rates[2] .* (1.0 .- exp.(-bd_rates[2] * t))
prob_data = Dict("num_sims_for_mean" => Nsims, "expected_mean_at_t" => expected_avg,
                 "jumps" => bd_jumps, "nu" => bd_nu)
"""
    Simple birth-death process with constant production and degradation.
    
    Reactions:
    1. 0 --> A (birth, rate k1=1000.0)
    2. A --> 0 (death, rate k2=10.0)
"""
prob_jump_constproduct = JumpProblemNetwork(bd_rates, tf, u0, prob, prob_data)

# Nonlinear reactions - direct jump implementation
# Species: A=1, B=2, C=3
nonlin_rates = [1.0, 2.0, 0.5, 0.75, 0.25]

function nonlin_jumps(u, p, t)
    A, B, C = u
    k1, k2, k3, k4, k5 = p
    [
        k1 * A * (A - 1) / 2,  # k1: 2A --> B (binomial for 2nd order)
        k2 * B,                # k2: B --> 2A  
        k3 * A * B,            # k3: A + B --> C
        k4 * C,                # k4: C --> A + B
        k5 * C * (C - 1) * (C - 2) / 6  # k5: 3C --> 3A (binomial for 3rd order)
    ]
end

nonlin_nu = [
    [-2, 1, 0],   # k1: 2A --> B
    [2, -1, 0],   # k2: B --> 2A  
    [-1, -1, 1],  # k3: A + B --> C
    [1, 1, -1],   # k4: C --> A + B
    [3, 0, -3]    # k5: 3C --> 3A
]

tf = 0.01
u0 = [200, 100, 150]  # [A, B, C]
prob = DiscreteProblem(u0, (0.0, tf), nonlin_rates)
Nsims = 32000
expected_avg = 84.876015624999994
prob_data = Dict("num_sims_for_mean" => Nsims, "expected_mean" => expected_avg,
                 "jumps" => nonlin_jumps, "nu" => nonlin_nu)
"""
    Example with a mix of nonlinear reactions, including third order
    
    Reactions:
    1. 2A --> B (dimerization, rate k1=1.0)
    2. B --> 2A (dissociation, rate k2=2.0)  
    3. A + B --> C (complex formation, rate k3=0.5)
    4. C --> A + B (complex dissociation, rate k4=0.75)
    5. 3C --> 3A (third order conversion, rate k5=0.25)
"""
prob_jump_nonlinrxs = JumpProblemNetwork(nonlin_rates, tf, u0, prob, prob_data)

# Oscillatory system with hill functions - direct jump implementation  
# Species: X=1, Y=2, Z=3, R=4, S=5, SP=6, SP2=7

# Hill function helper
hill_jump(x, k, n, h) = h > 0 ? (x/k)^n / (1 + (x/k)^n) : 1 / (1 + (x/k)^abs(n))

function oscil_jumps(u, p, t)
    X, Y, Z, R, S, SP, SP2 = u
    [
        0.01,                                          # 1: X --> 0 (degradation)
        0.01,                                          # 2: Y --> 0 (degradation)  
        0.01,                                          # 3: Z --> 0 (degradation)
        hill_jump(X, 3.0, 100.0, -4),                # 4: 0 --> Y (X-regulated)
        hill_jump(Y, 3.0, 100.0, -4),                # 5: 0 --> Z (Y-regulated)
        hill_jump(Z, 4.5, 100.0, -4),                # 6: 0 --> X (Z-regulated)  
        hill_jump(X, 2.0, 100.0, 6),                 # 7: 0 --> R (X-activated)
        hill_jump(Y, 15.0, 100.0, 4) * 0.002,        # 8: R --> 0 (Y-regulated)
        20.0,                                         # 9: 0 --> S (constant production)
        R * 0.005,                                    # 10: S --> SP (R-catalyzed)
        0.01 * SP * (SP - 1) / 2,                     # 11: SP + SP --> SP2 (dimerization)
        0.05                                          # 12: SP2 --> 0 (degradation)
    ]
end

oscil_nu = [
    [-1, 0, 0, 0, 0, 0, 0],   # 1: X --> 0
    [0, -1, 0, 0, 0, 0, 0],   # 2: Y --> 0 
    [0, 0, -1, 0, 0, 0, 0],   # 3: Z --> 0
    [0, 1, 0, 0, 0, 0, 0],    # 4: 0 --> Y
    [0, 0, 1, 0, 0, 0, 0],    # 5: 0 --> Z
    [1, 0, 0, 0, 0, 0, 0],    # 6: 0 --> X
    [0, 0, 0, 1, 0, 0, 0],    # 7: 0 --> R
    [0, 0, 0, -1, 0, 0, 0],   # 8: R --> 0
    [0, 0, 0, 0, 1, 0, 0],    # 9: 0 --> S
    [0, 0, 0, 0, -1, 1, 0],   # 10: S --> SP
    [0, 0, 0, 0, 0, -2, 1],   # 11: SP + SP --> SP2
    [0, 0, 0, 0, 0, 0, -1]    # 12: SP2 --> 0
]

u0 = [200.0, 60.0, 120.0, 100.0, 50.0, 50.0, 50.0]  # [X, Y, Z, R, S, SP, SP2]
tf = 4000.0
prob = DiscreteProblem(u0, (0.0, tf), Float64[])
prob_data = Dict("jumps" => oscil_jumps, "nu" => oscil_nu)
"""
    Oscillatory system, uses a mixture of jump types.
    
    Reactions include hill function regulation:
    - X, Y, Z degradation (rate 0.01 each)
    - Hill-regulated production: X→Y, Y→Z, Z→X  
    - X-activated R production
    - Y-regulated R degradation
    - Constant S production (rate 20)
    - R-catalyzed S→SP conversion
    - SP dimerization to SP2
    - SP2 degradation
"""
prob_jump_osc_mixed_jumptypes = JumpProblemNetwork(Float64[], tf, u0, prob, prob_data)

# Multistate model - direct jump implementation
# Species: S1=1, S2=2, S3=3, S4=4, S5=5, S6=6, S7=7, S8=8, S9=9
specs_sym_to_name = Dict(:S1 => "R(a,l)",
    :S2 => "L(r)",
    :S3 => "A(Y~U,r)",
    :S4 => "L(r!1).R(a,l!1)",
    :S5 => "A(Y~U,r!1).R(a!1,l)",
    :S6 => "A(Y~U,r!1).L(r!2).R(a!1,l!2)",
    :S7 => "A(Y~P,r!1).L(r!2).R(a!1,l!2)",
    :S8 => "A(Y~P,r!1).R(a!1,l)",
    :S9 => "A(Y~P,r)")
rsi = Dict(:R0 => 1, :L0 => 2, :A0 => 3, :kon => 4, :koff => 5,
    :kAon => 6, :kAoff => 7, :kAp => 8, :kAdp => 9)
params = (5360, 1160, 5360, 0.01, 0.1, 0.01, 0.1, 0.01, 0.1)

# Rate constants
multistate_rates = [0.01, 0.01, 0.1, 0.1, 0.01, 0.1]  # [kon, kAon, koff, kAoff, kAp, kAdp]

function multistate_jumps(u, p, t)
    S1, S2, S3, S4, S5, S6, S7, S8, S9 = u
    kon, kAon, koff, kAoff, kAp, kAdp = p
    [
        kon * S1 * S2,      # 1: S1 + S2 --> S4
        kAon * S1 * S3,     # 2: S1 + S3 --> S5  
        kon * S2 * S5,      # 3: S2 + S5 --> S6
        koff * S4,          # 4: S4 --> S1 + S2
        kAon * S3 * S4,     # 5: S3 + S4 --> S6
        kAoff * S5,         # 6: S5 --> S1 + S3
        koff * S6,          # 7: S6 --> S2 + S5
        kAoff * S6,         # 8: S6 --> S3 + S4
        kAp * S6,           # 9: S6 --> S7
        koff * S7,          # 10: S7 --> S2 + S8
        kAoff * S7,         # 11: S7 --> S4 + S9
        kAdp * S7,          # 12: S7 --> S6
        kon * S2 * S8,      # 13: S2 + S8 --> S7
        kAon * S1 * S9,     # 14: S1 + S9 --> S8
        kAon * S4 * S9,     # 15: S4 + S9 --> S7
        kAoff * S8,         # 16: S8 --> S1 + S9
        kAdp * S8,          # 17: S8 --> S5
        kAdp * S9           # 18: S9 --> S3
    ]
end

# State change vectors for each reaction
multistate_nu = [
    [-1, -1, 0, 1, 0, 0, 0, 0, 0],    # 1: S1 + S2 --> S4
    [-1, 0, -1, 0, 1, 0, 0, 0, 0],    # 2: S1 + S3 --> S5  
    [0, -1, 0, 0, -1, 1, 0, 0, 0],    # 3: S2 + S5 --> S6
    [1, 1, 0, -1, 0, 0, 0, 0, 0],     # 4: S4 --> S1 + S2
    [0, 0, -1, -1, 0, 1, 0, 0, 0],    # 5: S3 + S4 --> S6
    [1, 0, 1, 0, -1, 0, 0, 0, 0],     # 6: S5 --> S1 + S3
    [0, 1, 0, 0, 1, -1, 0, 0, 0],     # 7: S6 --> S2 + S5
    [0, 0, 1, 1, 0, -1, 0, 0, 0],     # 8: S6 --> S3 + S4
    [0, 0, 0, 0, 0, -1, 1, 0, 0],     # 9: S6 --> S7
    [0, 1, 0, 0, 0, 0, -1, 1, 0],     # 10: S7 --> S2 + S8
    [0, 0, 0, 1, 0, 0, -1, 0, 1],     # 11: S7 --> S4 + S9
    [0, 0, 0, 0, 0, 1, -1, 0, 0],     # 12: S7 --> S6
    [0, -1, 0, 0, 0, 0, 1, -1, 0],    # 13: S2 + S8 --> S7
    [-1, 0, 0, 0, 0, 0, 0, 1, -1],    # 14: S1 + S9 --> S8
    [0, 0, 0, -1, 0, 0, 1, 0, -1],    # 15: S4 + S9 --> S7
    [1, 0, 0, 0, 0, 0, 0, -1, 1],     # 16: S8 --> S1 + S9
    [0, 0, 0, 0, 1, 0, 0, -1, 0],     # 17: S8 --> S5
    [0, 0, 1, 0, 0, 0, 0, 0, -1]      # 18: S9 --> S3
]

u0 = [5360, 1160, 5360, 0, 0, 0, 0, 0, 0]  # [S1, S2, S3, S4, S5, S6, S7, S8, S9]
tf = 100.0
prob = DiscreteProblem(u0, (0.0, tf), multistate_rates)
"""
    Multistate model from Gupta and Mendes,
    "An Overview of Network-Based and -Free Approaches for Stochastic Simulation of Biochemical Systems",
    Computation 2018, 6, 9; doi:10.3390/computation6010009
    Translated from supplementary data file: Models/Multi-state/fixed_multistate.xml
    
    Reactions:
    1. S1 + S2 --> S4 (kon * S1 * S2)
    2. S1 + S3 --> S5 (kAon * S1 * S3)  
    3. S2 + S5 --> S6 (kon * S2 * S5)
    4. S4 --> S1 + S2 (koff * S4)
    5. S3 + S4 --> S6 (kAon * S3 * S4)
    6. S5 --> S1 + S3 (kAoff * S5)
    7. S6 --> S2 + S5 (koff * S6)
    8. S6 --> S3 + S4 (kAoff * S6)
    9. S6 --> S7 (kAp * S6)
    10. S7 --> S2 + S8 (koff * S7)
    11. S7 --> S4 + S9 (kAoff * S7)
    12. S7 --> S6 (kAdp * S7)
    13. S2 + S8 --> S7 (kon * S2 * S8)
    14. S1 + S9 --> S8 (kAon * S1 * S9)
    15. S4 + S9 --> S7 (kAon * S4 * S9)
    16. S8 --> S1 + S9 (kAoff * S8)
    17. S8 --> S5 (kAdp * S8)
    18. S9 --> S3 (kAdp * S9)
"""
prob_jump_multistate = JumpProblemNetwork(multistate_rates, tf, u0, prob,
    Dict("specs_to_sym_name" => specs_sym_to_name,
        "rates_sym_to_idx" => rsi, "params" => params,
        "jumps" => multistate_jumps, "nu" => multistate_nu))

# Twenty gene network - direct jump implementation
# Species: G[1:20], M[1:20], P[1:20], G_ind[1:20] (total 80 species)
N = 10  # number of genes

# No rate constants vector needed - rates are embedded in jump function
twentygenes_rates = Float64[]

function twentygenes_jumps(u, p, t)
    # u has 80 elements: G[1:20], M[1:20], P[1:20], G_ind[1:20]
    G = u[1:20]
    M = u[21:40] 
    P = u[41:60]
    G_ind = u[61:80]
    
    propensities = Float64[]
    
    for i in 1:N
        # Odd genes (2*i-1)
        odd_idx = 2*i - 1
        push!(propensities, 10.0 * G[odd_idx])          # G --> G + M (mRNA production)
        push!(propensities, 10.0 * M[odd_idx])          # M --> M + P (protein production)  
        push!(propensities, 1.0 * M[odd_idx])           # M --> 0 (mRNA degradation)
        push!(propensities, 1.0 * P[odd_idx])           # P --> 0 (protein degradation)
        
        # Even genes (2*i)
        even_idx = 2*i
        push!(propensities, 5.0 * G[even_idx])          # G --> G + M (mRNA production)
        push!(propensities, 5.0 * M[even_idx])          # M --> M + P (protein production)
        push!(propensities, 1.0 * M[even_idx])          # M --> 0 (mRNA degradation)  
        push!(propensities, 1.0 * P[even_idx])          # P --> 0 (protein degradation)
        
        # Cross-regulation reactions
        push!(propensities, 0.0001 * G[even_idx] * P[odd_idx])  # G + P --> G_ind
        push!(propensities, 100.0 * G_ind[even_idx])            # G_ind --> G_ind + M
    end
    
    propensities
end

function twentygenes_nu(N)
    # Build state change vectors for 80 species
    nu_matrix = []
    
    for i in 1:N
        odd_idx = 2*i - 1
        even_idx = 2*i
        
        # Odd gene reactions
        # G --> G + M (mRNA production) 
        delta = zeros(Int, 80)
        delta[20 + odd_idx] = 1  # M increases
        push!(nu_matrix, delta)
        
        # M --> M + P (protein production)
        delta = zeros(Int, 80) 
        delta[40 + odd_idx] = 1  # P increases
        push!(nu_matrix, delta)
        
        # M --> 0 (mRNA degradation)
        delta = zeros(Int, 80)
        delta[20 + odd_idx] = -1  # M decreases
        push!(nu_matrix, delta)
        
        # P --> 0 (protein degradation)
        delta = zeros(Int, 80)
        delta[40 + odd_idx] = -1  # P decreases  
        push!(nu_matrix, delta)
        
        # Even gene reactions
        # G --> G + M (mRNA production)
        delta = zeros(Int, 80)
        delta[20 + even_idx] = 1  # M increases
        push!(nu_matrix, delta)
        
        # M --> M + P (protein production)
        delta = zeros(Int, 80)
        delta[40 + even_idx] = 1  # P increases
        push!(nu_matrix, delta)
        
        # M --> 0 (mRNA degradation)
        delta = zeros(Int, 80)
        delta[20 + even_idx] = -1  # M decreases
        push!(nu_matrix, delta)
        
        # P --> 0 (protein degradation)
        delta = zeros(Int, 80)
        delta[40 + even_idx] = -1  # P decreases
        push!(nu_matrix, delta)
        
        # Cross-regulation reactions
        # G + P --> G_ind
        delta = zeros(Int, 80)
        delta[even_idx] = -1          # G decreases
        delta[40 + odd_idx] = -1      # P decreases  
        delta[60 + even_idx] = 1      # G_ind increases
        push!(nu_matrix, delta)
        
        # G_ind --> G_ind + M
        delta = zeros(Int, 80)
        delta[20 + even_idx] = 1      # M increases
        push!(nu_matrix, delta)
    end
    
    nu_matrix
end

twentygenes_nu_matrix = twentygenes_nu(N)

# Initial conditions: G[i] = 1 for all genes, everything else = 0
u0 = zeros(Int, 80)
u0[1:20] .= 1  # All G species start at 1

tf = 2000.0
prob = DiscreteProblem(u0, (0.0, tf), twentygenes_rates)

"""
    Twenty-gene model from McCollum et al,
    "The sorting direct method for stochastic simulation of biochemical systems with varying reaction execution behavior"
    Comp. Bio. and Chem., 30, pg. 39-49 (2006).
    
    For each gene pair i (i = 1 to 10):
    Odd genes (2*i-1):
    - G --> G + M (rate 10.0)
    - M --> M + P (rate 10.0)  
    - M --> 0 (rate 1.0)
    - P --> 0 (rate 1.0)
    
    Even genes (2*i):
    - G --> G + M (rate 5.0)
    - M --> M + P (rate 5.0)
    - M --> 0 (rate 1.0)
    - P --> 0 (rate 1.0)
    
    Cross-regulation:
    - G(2*i) + P(2*i-1) --> G_ind(2*i) (rate 0.0001)
    - G_ind(2*i) --> G_ind(2*i) + M(2*i) (rate 100.0)
"""
prob_jump_twentygenes = JumpProblemNetwork(twentygenes_rates, tf, u0, prob, 
    Dict("jumps" => twentygenes_jumps, "nu" => twentygenes_nu_matrix))

# DNA dimer repressor - direct jump implementation  
# Species: G=1, M=2, P=3, P2=4, P2G=5
dnadimer_rates = [0.09, 0.05, 0.001, 0.0009, 0.00001, 0.0005, 0.005, 0.9]

function dnadimer_jumps(u, p, t)
    G, M, P, P2, P2G = u
    c1, c2, c3, c4, c5, c6, c7, c8 = p
    [
        c1 * G,                    # c1: G --> G + M (mRNA production)
        c2 * M,                    # c2: M --> M + P (protein production)
        c3 * M,                    # c3: M --> 0 (mRNA degradation)  
        c4 * P,                    # c4: P --> 0 (protein degradation)
        c5 * P * (P - 1) / 2,      # c5: 2P --> P2 (dimerization)
        c6 * P2,                   # c6: P2 --> 2P (dimer dissociation)
        c7 * P2 * G,               # c7: P2 + G --> P2G (repressor binding)
        c8 * P2G                   # c8: P2G --> P2 + G (repressor unbinding)
    ]
end

# State change vectors for each reaction
dnadimer_nu = [
    [0, 1, 0, 0, 0],    # c1: G --> G + M (mRNA production)
    [0, 0, 1, 0, 0],    # c2: M --> M + P (protein production)
    [0, -1, 0, 0, 0],   # c3: M --> 0 (mRNA degradation)
    [0, 0, -1, 0, 0],   # c4: P --> 0 (protein degradation)
    [0, 0, -2, 1, 0],   # c5: 2P --> P2 (dimerization)
    [0, 0, 2, -1, 0],   # c6: P2 --> 2P (dimer dissociation)
    [-1, 0, 0, -1, 1],  # c7: P2 + G --> P2G (repressor binding) 
    [1, 0, 0, 1, -1]    # c8: P2G --> P2 + G (repressor unbinding)
]

varlabels = ["G", "M", "P", "P2", "P2G"]
u0 = [1000, 0, 0, 0, 0]  # [G, M, P, P2, P2G]
tf = 4000.0
prob = DiscreteProblem(u0, (0.0, tf), dnadimer_rates)
"""
    Negative feedback autoregulatory gene expression model. Dimer is the repressor.
    Taken from Marchetti, Priami and Thanh,
    "Simulation Algorithms for Comptuational Systems Biology",
    Springer (2017).
    
    Reactions:
    1. G --> G + M (mRNA production, rate c1=0.09)
    2. M --> M + P (protein production, rate c2=0.05)
    3. M --> 0 (mRNA degradation, rate c3=0.001)
    4. P --> 0 (protein degradation, rate c4=0.0009)
    5. 2P --> P2 (dimerization, rate c5=0.00001)
    6. P2 --> 2P (dimer dissociation, rate c6=0.0005)
    7. P2 + G --> P2G (repressor binding, rate c7=0.005)
    8. P2G --> P2 + G (repressor unbinding, rate c8=0.9)
"""
prob_jump_dnadimer_repressor = JumpProblemNetwork(dnadimer_rates, tf, u0, prob,
    Dict("specs_names" => varlabels, "jumps" => dnadimer_jumps, "nu" => dnadimer_nu))

# Diffusion network - direct jump implementation
# Creates a 1D diffusion lattice with N sites and nearest-neighbor hopping

function getDiffJumps(N)
    function diffusion_jumps(u, p, t)
        K = p[1]  # diffusion rate
        propensities = Float64[]
        
        # For each adjacent pair of sites, add forward and backward jumps
        for i in 1:(N-1)
            push!(propensities, K * u[i])      # X[i] --> X[i+1]
            push!(propensities, K * u[i+1])    # X[i+1] --> X[i]  
        end
        
        propensities
    end
    diffusion_jumps
end

function getDiffNu(N)
    nu_matrix = []
    
    # For each adjacent pair of sites, create state change vectors
    for i in 1:(N-1)
        # X[i] --> X[i+1]
        delta = zeros(Int, N)
        delta[i] = -1      # X[i] decreases
        delta[i+1] = 1     # X[i+1] increases
        push!(nu_matrix, delta)
        
        # X[i+1] --> X[i]  
        delta = zeros(Int, N)
        delta[i+1] = -1    # X[i+1] decreases
        delta[i] = 1       # X[i] increases
        push!(nu_matrix, delta)
    end
    
    nu_matrix
end

function getDiffNetwork(N)
    # Return jump functions and nu matrix for N-site diffusion lattice
    jumps = getDiffJumps(N)
    nu = getDiffNu(N)
    rates = [1.0]  # K = 1.0
    return (jumps, nu, rates)
end

params = [1.0]  # K = 1.0
function getDiffu0(diffnetwork, N)
    10 * ones(Int64, N)  # 10 particles per site initially
end

tf = 10.0
"""
    Continuous time random walk (i.e. diffusion approximation) example.
    Creates a 1D lattice with N sites where particles can hop between
    adjacent sites with rate K.
    
    For each adjacent pair of sites i and i+1:
    - X[i] --> X[i+1] (rate K * X[i])
    - X[i+1] --> X[i] (rate K * X[i+1])
    
    The network function returns jump functions and state change vectors
    for a lattice of given size N.
"""
prob_jump_diffnetwork = JumpProblemNetwork(params, tf, getDiffu0, nothing,
    Dict("network_func" => getDiffNetwork))

end # module

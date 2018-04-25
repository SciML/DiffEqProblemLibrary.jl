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

"""
    DNA negative feedback autoregulatory model. Protein acts as repressor.    
"""
rs = @reaction_network ptype begin
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
prob_jump_dnarepressor = JumpProblemNetwork(rs, rates, tf, u0, prob, prob_data)


"""
    Simple constant production with degradation example
"""
rs = @reaction_network pdtype begin
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
prob_jump_constproduct = JumpProblemNetwork(rs, rates, tf, u0, prob, prob_data)


"""
    Example with a mix of nonlinear reactions, including third order
"""
rs = @reaction_network dtype begin
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
prob_jump_nonlinrxs = JumpProblemNetwork(rs, rates, tf, u0, prob, prob_data)


"""
    Oscillatory system, uses a mixture of jump types.
"""
rs = @reaction_network rnType  begin
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
prob_jump_osc_mixed_jumptypes = JumpProblemNetwork(rs, nothing, tf, u0, prob, nothing)
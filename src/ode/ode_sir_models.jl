function SIR(du,u,p,t)
    S = u[1] ; I = u[2] ; R = u[3];
    β = p[1] ; γ = p[2] ; N = p[3]
    du[1] = -β*S*I/N
    du[2] = β*S*I/N - γ*I/N
    du[3] = γ*I/N
end
p = [0.25,0.1,1000]
u₀ = [999,1,0]
tspan = (0,200)
prob_sir = ODEProblem(SIR,u₀,tspan,p)

function SIRD(du,u,p,t)
     S = u[1] ; I = u[2] ; R = u[3] ; D = u[4]
     β = p[1] ; γ = p[2] ; δ = p[3] ; N = p[4]

     du[1] = -β*S*I/N
     du[2] = β*S*I/N - γ*I - δ*I
     du[3] = γ*I
     du[4] = δ*I
end
p = [0.25,0.1,0.01,1000]
u₀ = [999,1,0,0]
tspan = (0,200)
prob_sird = ODEProblem(SIRD,u₀,tspan,p)

function SEIR(du,u,p,t)
    S = u[1] ; E = u[2] ; I = u[3] ; R = u[4]
    μ = p[1] ; β = p[2] ; α = p[3] ; γ = p[4] ; N = p[5]

    du[1] = μ*N - μ*S -β/N*S*I
    du[2] = β*I*S/N - (μ+α)*E
    du[3] = α*E - γ*I - μ*I
    du[4] = γ*I - μ*R
end
p = [0.01,0.25,0.01,0.1,1000]
u₀ = [900,0,100,0]
tspan = (0,200)
prob_seir = ODEProblem(SEIR,u₀,tspan,p)

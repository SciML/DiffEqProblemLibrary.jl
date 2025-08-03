#!/usr/bin/env julia

# Test script to verify SDE problem conversion works correctly
using Pkg
Pkg.activate(".")

println("Testing SDE problem conversion...")

# Add StochasticDiffEq if not present
try
    using StochasticDiffEq
    println("✓ StochasticDiffEq.jl loaded successfully")
catch
    println("Installing StochasticDiffEq.jl...")
    Pkg.add("StochasticDiffEq")
    using StochasticDiffEq
    println("✓ StochasticDiffEq.jl installed and loaded")
end

# Load our converted SDEProblemLibrary
include("lib/SDEProblemLibrary/src/SDEProblemLibrary.jl")
using .SDEProblemLibrary

println("\n=== Testing converted SDE problems ===")

# Test the converted oscillatory reaction problem
println("\n1. Testing prob_sde_oscilreact...")
try
    prob = prob_sde_oscilreact
    println("   ✓ Problem loaded: $(typeof(prob))")
    println("   ✓ Initial condition: $(prob.u0)")
    println("   ✓ Time span: $(prob.tspan)")
    println("   ✓ Parameters: $(prob.p)")
    
    # Test a short simulation
    println("   ✓ Running short simulation...")
    sol = solve(prob, EM(), dt=0.1, save_everystep=false, save_end=true)
    println("   ✓ Simulation completed successfully")
    println("   ✓ Final state: $(sol.u[end])")
    
catch e
    println("   ✗ Error with prob_sde_oscilreact: $e")
end

# Test a few other SDE problems to ensure they still work
println("\n2. Testing prob_sde_linear...")
try
    prob = prob_sde_linear
    sol = solve(prob, EM(), dt=0.01, save_everystep=false)
    println("   ✓ prob_sde_linear works correctly")
catch e
    println("   ✗ Error with prob_sde_linear: $e")
end

println("\n3. Testing prob_sde_lorenz...")
try
    prob = prob_sde_lorenz
    sol = solve(prob, EM(), dt=0.01, save_everystep=false, save_end=true)
    println("   ✓ prob_sde_lorenz works correctly")
    println("   ✓ Final state: $(sol.u[end])")
catch e
    println("   ✗ Error with prob_sde_lorenz: $e")
end

println("\n4. Testing prob_sde_bruss...")
try
    prob = prob_sde_bruss
    sol = solve(prob, EM(), dt=0.01, save_everystep=false, save_end=true)
    println("   ✓ prob_sde_bruss works correctly")
    println("   ✓ Final state: $(sol.u[end])")
catch e
    println("   ✗ Error with prob_sde_bruss: $e")
end

println("\n=== SDE testing completed ===")
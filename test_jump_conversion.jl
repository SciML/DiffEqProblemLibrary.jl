#!/usr/bin/env julia

# Test script to verify Jump problem conversion works correctly
using Pkg
Pkg.activate(".")

println("Testing Jump problem conversion...")

# Add JumpProcesses if not present
try
    using JumpProcesses
    println("✓ JumpProcesses.jl loaded successfully")
catch
    println("Installing JumpProcesses.jl...")
    Pkg.add("JumpProcesses")
    using JumpProcesses
    println("✓ JumpProcesses.jl installed and loaded")
end

# Load our converted JumpProblemLibrary
include("lib/JumpProblemLibrary/src/JumpProblemLibrary.jl")
using .JumpProblemLibrary

println("\n=== Testing converted Jump problems ===")

# Test the converted jump problems
println("\n1. Testing prob_jump_dnarepressor...")
try
    prob_network = prob_jump_dnarepressor
    println("   ✓ Network loaded: $(typeof(prob_network))")
    println("   ✓ Initial condition: $(prob_network.u0)")
    println("   ✓ Time span: (0.0, $(prob_network.tstop))")
    println("   ✓ Rates: $(prob_network.rates)")
    
    # Create jumps using the stored functions
    if haskey(prob_network.prob_data, "jumps")
        jumps_func = prob_network.prob_data["jumps"]
        nu = prob_network.prob_data["nu"]
        
        # Create MassActionJumps manually for testing
        println("   ✓ Testing jump rate calculations...")
        rates = jumps_func(prob_network.u0, prob_network.rates, 0.0)
        println("   ✓ Initial jump rates: $(rates)")
        
        # Test basic functionality
        if all(r >= 0 for r in rates)
            println("   ✓ All jump rates are non-negative")
        else
            println("   ⚠ Some jump rates are negative: $(rates)")
        end
    else
        println("   ⚠ No jump functions found in prob_data")
    end
    
catch e
    println("   ✗ Error with prob_jump_dnarepressor: $e")
end

println("\n2. Testing prob_jump_constproduct...")
try
    prob_network = prob_jump_constproduct
    println("   ✓ Network loaded successfully")
    
    if haskey(prob_network.prob_data, "jumps")
        jumps_func = prob_network.prob_data["jumps"] 
        rates = jumps_func(prob_network.u0, prob_network.rates, 0.0)
        println("   ✓ Jump rates: $(rates)")
        
        # Test expected mean function
        if haskey(prob_network.prob_data, "expected_mean_at_t")
            expected_func = prob_network.prob_data["expected_mean_at_t"]
            mean_at_1 = expected_func(1.0)
            println("   ✓ Expected mean at t=1: $(mean_at_1)")
        end
    end
    
catch e
    println("   ✗ Error with prob_jump_constproduct: $e")
end

println("\n3. Testing prob_jump_nonlinrxs...")
try
    prob_network = prob_jump_nonlinrxs
    println("   ✓ Network loaded successfully")
    
    if haskey(prob_network.prob_data, "jumps")
        jumps_func = prob_network.prob_data["jumps"]
        rates = jumps_func(prob_network.u0, prob_network.rates, 0.0)
        println("   ✓ Jump rates: $(rates)")
        println("   ✓ Initial condition: $(prob_network.u0)")
    end
    
catch e
    println("   ✗ Error with prob_jump_nonlinrxs: $e")
end

println("\n4. Testing prob_jump_osc_mixed_jumptypes...")
try
    prob_network = prob_jump_osc_mixed_jumptypes
    println("   ✓ Network loaded successfully")
    
    if haskey(prob_network.prob_data, "jumps")
        jumps_func = prob_network.prob_data["jumps"]
        rates = jumps_func(prob_network.u0, Float64[], 0.0)  # No parameters for this one
        println("   ✓ Jump rates: $(rates)")
        println("   ✓ Number of species: $(length(prob_network.u0))")
        println("   ✓ Number of reactions: $(length(rates))")
    end
    
catch e
    println("   ✗ Error with prob_jump_osc_mixed_jumptypes: $e")
end

println("\n=== Jump testing completed ===")
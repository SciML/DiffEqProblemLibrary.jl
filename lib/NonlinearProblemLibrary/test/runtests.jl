# The test is simply that all of the examples build!
using NonlinearProblemLibrary
using Test

# Check that there are no undefined exports, stale dependencies, etc.
# Ambiguity checks are disabled since tests fail due to ambiguities
# in dependencies
using Aqua
Aqua.test_all(NonlinearProblemLibrary; ambiguities = false)

# Allocation tests - ensure key functions don't allocate
if get(ENV, "GROUP", "all") == "all" || get(ENV, "GROUP", "all") == "nopre"
    @testset "Allocation Tests" begin
        using AllocCheck

        # Test p23_f! (Chandrasekhar function) - the fixed function
        @testset "p23_f! zero allocations" begin
            x = ones(10)
            out = zeros(10)
            # Warmup
            NonlinearProblemLibrary.p23_f!(out, x)
            # Test
            allocs = @allocated NonlinearProblemLibrary.p23_f!(out, x)
            @test allocs == 0
        end

        # Test other key functions
        @testset "p1_f! zero allocations" begin
            x = ones(10)
            x[1] = -1.2
            out = zeros(10)
            NonlinearProblemLibrary.p1_f!(out, x)
            allocs = @allocated NonlinearProblemLibrary.p1_f!(out, x)
            @test allocs == 0
        end

        @testset "p2_f! zero allocations" begin
            x = [3.0, -1.0, 0.0, 1.0]
            out = zeros(4)
            NonlinearProblemLibrary.p2_f!(out, x)
            allocs = @allocated NonlinearProblemLibrary.p2_f!(out, x)
            @test allocs == 0
        end
    end
end

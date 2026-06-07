using NonlinearProblemLibrary
using Test

const TEST_GROUP = get(ENV, "DIFFEQPROBLEMLIBRARY_TEST_GROUP", "All")

if TEST_GROUP == "Core" || TEST_GROUP == "All"
    @time @testset "Load Tests" begin
        @test NonlinearProblemLibrary isa Module
    end

    # Allocation tests - ensure key functions don't allocate.
    @time @testset "Allocation Tests" begin
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

# Quality assurance: no undefined exports, stale dependencies, etc.
# Ambiguity checks are disabled since tests fail due to ambiguities in dependencies.
if TEST_GROUP == "QA" || TEST_GROUP == "All"
    using Aqua
    @time @testset "Aqua" begin
        Aqua.test_all(NonlinearProblemLibrary; ambiguities = false)
    end
end

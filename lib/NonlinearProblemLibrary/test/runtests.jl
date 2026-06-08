using Pkg
using NonlinearProblemLibrary
using SafeTestsets, Test

const TEST_GROUP = get(ENV, "DIFFEQPROBLEMLIBRARY_TEST_GROUP", "All")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

if TEST_GROUP == "Core" || TEST_GROUP == "All"
    @time @testset "Load Tests" begin
        @test NonlinearProblemLibrary isa Module
    end
end

# Quality assurance: allocation checks (AllocCheck), then no undefined exports,
# stale dependencies, etc. (Aqua). These run in the isolated test/qa environment
# so the heavy tooling stays out of the main test env. Ambiguity checks are
# disabled since tests fail due to ambiguities in dependencies.
if TEST_GROUP == "QA" || TEST_GROUP == "All"
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end

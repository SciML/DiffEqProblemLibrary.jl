using Pkg
using DAEProblemLibrary
using SafeTestsets, Test

const TEST_GROUP = get(ENV, "DIFFEQPROBLEMLIBRARY_TEST_GROUP", "All")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

# The Core test is simply that all of the examples build (load the module).
if TEST_GROUP == "Core" || TEST_GROUP == "All"
    @time @testset "Load Tests" begin
        @test DAEProblemLibrary isa Module
    end
end

# Quality assurance: no undefined exports, stale dependencies, etc.
# Ambiguity checks are disabled since tests fail due to ambiguities in dependencies.
if TEST_GROUP == "QA" || TEST_GROUP == "All"
    activate_qa_env()
    @time @safetestset "Aqua" include("qa/qa.jl")
end

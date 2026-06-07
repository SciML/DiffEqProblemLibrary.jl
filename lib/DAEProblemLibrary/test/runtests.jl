using DAEProblemLibrary
using Test

const TEST_GROUP = get(ENV, "DIFFEQPROBLEMLIBRARY_TEST_GROUP", "ALL")

# The Core test is simply that all of the examples build (load the module).
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @testset "Load Tests" begin
        @test DAEProblemLibrary isa Module
    end
end

# Quality assurance: no undefined exports, stale dependencies, etc.
# Ambiguity checks are disabled since tests fail due to ambiguities in dependencies.
if TEST_GROUP == "QA" || TEST_GROUP == "ALL"
    using Aqua
    @time @testset "Aqua" begin
        Aqua.test_all(DAEProblemLibrary; ambiguities = false)
    end
end

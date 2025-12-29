# The test is simply that all of the examples build!
using DiffEqProblemLibrary
using ExplicitImports
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(DiffEqProblemLibrary) === nothing
    @test check_no_stale_explicit_imports(DiffEqProblemLibrary) === nothing
end

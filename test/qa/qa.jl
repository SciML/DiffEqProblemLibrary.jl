# The root umbrella package's own test: it re-exports each sublibrary,
# so the test is that the whole thing builds with no implicit or stale
# explicit imports.
using DiffEqProblemLibrary
using ExplicitImports
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(DiffEqProblemLibrary) === nothing
    @test check_no_stale_explicit_imports(DiffEqProblemLibrary) === nothing
end

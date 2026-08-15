using SciMLTesting, DiffEqProblemLibrary, Test

# extras-only harness deps; Aqua sees them after Pkg.test extras-merge
run_qa(
    DiffEqProblemLibrary;
    aqua_kwargs = (; stale_deps = (; ignore = [:SafeTestsets, :SciMLTesting])),
)

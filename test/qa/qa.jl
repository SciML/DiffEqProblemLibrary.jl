using SciMLTesting, DiffEqProblemLibrary, Test

# The root umbrella package re-exports each sublibrary; it owns no methods of its own,
# so the reported method ambiguities all come from the imported dependency tree
# (DiffEqBase / SciMLBase / Catalyst), matching the `ambiguities = false` each sublib uses.
run_qa(DiffEqProblemLibrary; explicit_imports = true, aqua_kwargs = (; ambiguities = false))

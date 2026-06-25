using SciMLTesting, SDEProblemLibrary, Test

# `RuntimeGeneratedFunctions.init(@__MODULE__)` qualified-accesses the non-public
# `RuntimeGeneratedFunctions.init` (the package's documented init boilerplate), so it
# is ignored for the public-API qualified-access check.
run_qa(
    SDEProblemLibrary; explicit_imports = true,
    aqua_kwargs = (; ambiguities = false, persistent_tasks = false),
    ei_kwargs = (; all_qualified_accesses_are_public = (; ignore = (:init,)))
)
